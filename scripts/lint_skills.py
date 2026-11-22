#!/usr/bin/env python3
"""Lint Claude Code skills under .claude/skills.

Checks (all run by default; --offline skips the nf-core registry check):
  frontmatter   name/description present, name == directory, name charset
  size          SKILL.md <= 500 lines (warn > 400); supporting .md <= 600
  rules         sub-skills carry the shared operating-rules block
  fences        code fences have a language; bash blocks have no comment/blank
                line inside a backslash continuation
  shellcheck    every bash block parses (shellcheck -S error), if installed
  deprecated    forbidden tool/flag patterns inside code blocks
  nfcore        every nf-core/<name> exists in https://nf-co.re/pipelines.json
  paths         every `.claude/skills/...` path mentioned exists
  settings      .claude/settings.json has no legacy "skills" key

Exit 1 on any error; warnings never fail.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / ".claude" / "skills"
NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
RULES_MARK = "Follow the eight rules in the `genomics` router skill"
DEPRECATED = [
    (re.compile(r"\bmacs2 callpeak"), "use macs3"),
    (re.compile(r"\blumpyexpress\b"), "use smoove"),
    (re.compile(r"\bgemini (load|query|de_novo|comp_hets)\b"), "use slivar"),
    (re.compile(r"library\(copynumber\)"), "package removed from Bioconductor"),
    (re.compile(r"nf-core/(copynumber|fastqc)\b"), "pipeline does not exist"),
    (re.compile(r"samtools view [^\n|]*-bS\b"), "-S is obsolete"),
    (re.compile(r"\bINFO/FILTER\b"), "FILTER is not an INFO field"),
    (re.compile(r"--plugin ClinVar\b"), "no such VEP plugin; use --custom"),
    (re.compile(r"--plugin gnomADc\b"), "gnomADc is coverage, not constraint"),
    (re.compile(r"^\s*HMMRATAC\b", re.M), "use macs3 hmmratac"),
    (re.compile(r"deepvariant:1\.[0-9]\.[0-9]"), "pin deepvariant >= 1.10.0"),
    (re.compile(r"\bmosdepth\b[^\n]*\s-b\s[^\n]*--by\b"), "-b and --by are the same flag"),
]

errors: list[str] = []
warnings: list[str] = []


def err(msg: str) -> None:
    errors.append(msg)


def warn(msg: str) -> None:
    warnings.append(msg)


def parse_frontmatter(text: str, path: Path) -> dict[str, str]:
    if not text.startswith("---\n"):
        err(f"{path}: missing YAML frontmatter")
        return {}
    end = text.find("\n---", 4)
    if end < 0:
        err(f"{path}: unterminated frontmatter")
        return {}
    fm: dict[str, str] = {}
    key = None
    for line in text[4:end].splitlines():
        if not line.strip():
            continue
        m = re.match(r"^([A-Za-z_][\w-]*):\s*(.*)$", line)
        if m:
            key, val = m.group(1), m.group(2).strip()
            fm[key] = val.strip('"').strip("'")
        elif line.lstrip().startswith("-") and key:
            fm[key] = (fm.get(key, "") + " " + line.strip()).strip()
        elif key and line.startswith((" ", "\t")):
            fm[key] = (fm.get(key, "") + " " + line.strip()).strip()
        else:
            err(f"{path}: unparseable frontmatter line: {line!r}")
    return fm


def code_blocks(text: str, path: Path):
    """Yield (lang, body, start_line) for fenced blocks."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^(\s*)(```+|~~~+)\s*([\w+-]*)(?:\s[^`]*)?$", lines[i])
        if not m:
            i += 1
            continue
        indent, fence, lang = m.group(1), m.group(2), m.group(3)
        start = i + 1
        j = i + 1
        while j < len(lines) and not re.match(rf"^{re.escape(indent)}{re.escape(fence)}\s*$", lines[j]):
            j += 1
        body = "\n".join(l[len(indent):] if l.startswith(indent) else l for l in lines[start:j])
        if not lang:
            err(f"{path}:{start}: code fence without language")
        yield lang, body, start
        i = j + 1


def is_command(row: str) -> bool:
    """A line that can fail: not blank, not a comment, not a plain variable assignment."""
    s = row.strip()
    if not s or s.startswith("#"):
        return False
    return not re.match(r"^[A-Za-z_][A-Za-z0-9_]*=(?!.*\$\()", s)


def check_bash_block(body: str, path: Path, start: int, shellcheck: str | None) -> None:
    rows = body.splitlines()
    for k, row in enumerate(rows[:-1]):
        if row.rstrip().endswith("\\"):
            nxt = rows[k + 1]
            if not nxt.strip():
                err(f"{path}:{start + k + 1}: blank line after backslash continuation")
            elif nxt.lstrip().startswith("#"):
                err(f"{path}:{start + k + 1}: comment inside backslash continuation")
            elif nxt.strip() == "\\":
                err(f"{path}:{start + k + 1}: bare backslash line")
    # Rule 4: a block with two or more top-level command lines must start with set -euo pipefail.
    # Exempt: blocks that open with "# alternatives" (independent one-liners, run one) and sbatch headers.
    logical = re.sub(r"\\\n\s*", " ", body).splitlines()
    cmds = [r for r in logical if is_command(r)]
    alternatives = bool(rows) and rows[0].strip().lower().startswith("# alternatives")
    if len(cmds) >= 2 and not alternatives and not any(r.strip().startswith("set -euo pipefail") for r in cmds[:2]):
        if not any(r.strip().startswith("#SBATCH") for r in logical[:3]):
            err(f"{path}:{start}: multi-command bash block without `set -euo pipefail` as its first command (or open it with `# alternatives: run one`)")
    if shellcheck:
        with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
            fh.write("#!/usr/bin/env bash\n" + body + "\n")
            tmp = fh.name
        res = subprocess.run(
            [shellcheck, "-S", "error", "-e", "SC2154,SC1091,SC2148,SC1090", "-f", "gcc", tmp],
            capture_output=True, text=True,
        )
        Path(tmp).unlink(missing_ok=True)
        for line in res.stdout.splitlines():
            m = re.match(r"^[^:]+:(\d+):(\d+): error: (.*)$", line)
            if m:
                err(f"{path}:{start + int(m.group(1)) - 1}: shellcheck: {m.group(3)}")


def main(argv: list[str]) -> int:
    offline = "--offline" in argv
    shellcheck = shutil.which("shellcheck")
    if not shellcheck:
        warn("shellcheck not installed; bash blocks are not parsed")

    if not SKILLS.is_dir():
        err(f"{SKILLS} missing")
        return report()

    settings = ROOT / ".claude" / "settings.json"
    if settings.exists():
        try:
            if "skills" in json.loads(settings.read_text()):
                err(f"{settings}: legacy 'skills' key is not a Claude Code setting; remove it")
        except json.JSONDecodeError as e:
            err(f"{settings}: invalid JSON: {e}")

    skill_dirs = sorted(p for p in SKILLS.iterdir() if p.is_dir())
    stray = [p for p in SKILLS.iterdir() if p.is_file() and p.suffix == ".md"]
    for p in stray:
        err(f"{p}: flat skill file is not loaded by Claude Code; move to {p.stem}/SKILL.md")

    nfcore_mentions: dict[str, list[str]] = {}
    all_md: list[Path] = []
    for d in skill_dirs:
        skill = d / "SKILL.md"
        if not skill.exists():
            err(f"{d}: no SKILL.md")
            continue
        text = skill.read_text()
        fm = parse_frontmatter(text, skill)
        name = fm.get("name", "")
        if not name:
            err(f"{skill}: frontmatter lacks name")
        elif not NAME_RE.match(name):
            err(f"{skill}: invalid name {name!r}")
        elif name != d.name:
            err(f"{skill}: name {name!r} != directory {d.name!r}")
        desc = fm.get("description", "")
        if not desc:
            err(f"{skill}: frontmatter lacks description")
        elif len(desc) > 1024:
            err(f"{skill}: description > 1024 chars")
        elif "Use when" not in desc and "use when" not in desc:
            warn(f"{skill}: description has no 'Use when' trigger clause")
        for key in ("allowed-tools", "disallowed-tools", "model", "context", "agent", "hooks", "effort"):
            if key in fm:
                err(f"{skill}: frontmatter key {key!r} is not permitted in this repo (see CLAUDE.md)")
        if "CLAUDE_PROJECT_DIR" in text:
            err(f"{skill}: references ${{CLAUDE_PROJECT_DIR}}; resolve scripts via ${{CLAUDE_SKILL_DIR}} instead")
        inj = [l for l in text.splitlines() if l.startswith("!`")]
        for l in inj:
            if not re.match(r'^!`f="\$\{CLAUDE_SKILL_DIR\}(/\.\./genomics)?/scripts/check_env\.sh"; \[ -f "\$f" \] && bash "\$f" [A-Za-z0-9._+ -]+ \|\| true`$', l):
                err(f"{skill}: dynamic-context line must be the guarded check_env.sh form: {l[:80]}")
        n = text.count("\n") + 1
        if n > 500:
            err(f"{skill}: {n} lines > 500")
        elif n > 400:
            warn(f"{skill}: {n} lines > 400 (move material to a supporting file)")
        if d.name != "genomics" and RULES_MARK not in text:
            err(f"{skill}: missing shared operating-rules block")
        for sup in sorted(d.glob("*.md")):
            if sup.name == "SKILL.md":
                continue
            m = sup.read_text().count("\n") + 1
            if m > 600:
                err(f"{sup}: {m} lines > 600")
            if sup.name not in text:
                warn(f"{skill}: supporting file {sup.name} is never referenced")
        all_md.extend(sorted(d.glob("*.md")))
        for sh in d.glob("scripts/*.sh"):
            if shellcheck:
                res = subprocess.run([shellcheck, "-S", "warning", str(sh)], capture_output=True, text=True)
                if res.returncode:
                    err(f"{sh}: shellcheck\n{res.stdout.strip()}")

    for md in all_md:
        text = md.read_text()
        for lang, body, start in code_blocks(text, md):
            if lang in ("bash", "sh", "shell"):
                check_bash_block(body, md, start, shellcheck)
            # Join backslash continuations so `samtools view \<nl> -bS` and `mosdepth -b x \<nl> --by` are caught.
            joined = re.sub(r"\\\n\s*", " ", body)
            for rx, why in DEPRECATED:
                for m in rx.finditer(joined):
                    line = start + joined[: m.start()].count("\n")  # approximate when continuations precede the match
                    err(f"{md}:{line}: deprecated `{m.group(0).strip()}` ({why})")
        for m in re.finditer(r"nf-core/([a-z0-9]+)", text):
            ls = text.rfind("\n", 0, m.start()) + 1
            le = text.find("\n", m.end())
            line_text = text[ls:le if le >= 0 else len(text)]
            if re.search(r"(does|do) not exist|not found|There is no|no such|Neither ", line_text):
                continue  # prose warning the reader that the pipeline is fictional
            nfcore_mentions.setdefault(m.group(1), []).append(f"{md}:{text[:m.start()].count(chr(10)) + 1}")
        for m in re.finditer(r"`?(\.claude/skills/[\w./-]+?)`?(?=[\s`,;:)\]])", text):
            ref = m.group(1).rstrip(".")
            if not (ROOT / ref).exists():
                err(f"{md}:{text[:m.start()].count(chr(10)) + 1}: path does not exist: {ref}")

    if nfcore_mentions and not offline:
        try:
            with urllib.request.urlopen("https://nf-co.re/pipelines.json", timeout=30) as r:
                registry = {w["name"] for w in json.load(r)["remote_workflows"]}
            registry |= {"tools", "modules", "configs", "website"}
            for name, where in sorted(nfcore_mentions.items()):
                if name not in registry:
                    err(f"nf-core/{name} not in registry ({', '.join(where[:3])})")
        except Exception as e:  # network failure must not mask other findings
            warn(f"nf-core registry unreachable ({e}); skipped")

    return report()


def report() -> int:
    for w in warnings:
        print(f"WARN  {w}")
    for e in errors:
        print(f"ERROR {e}")
    print(f"{len(errors)} error(s), {len(warnings)} warning(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
