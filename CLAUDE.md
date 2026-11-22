# claude-skill-genomic-pipelines

Claude Code skills for genomics. This file is for agents editing **this repository**; the skills themselves are for agents doing genomics work.

## What is here

- `.claude/skills/<name>/SKILL.md` — one skill per directory, Claude Code format (YAML frontmatter + Markdown). `genomics` is the router; the twelve `genomics-*` directories are assay/task skills.
- `.claude/skills/<name>/*.md` — supporting files loaded on demand (`reference.md`, `bulk.md`, `atac.md`, …). Every supporting file must be referenced from its `SKILL.md`.
- `.claude/skills/genomics/scripts/` — `peek.sh`, `ref_sanity.sh`, `check_env.sh`; the only code the skills execute. Every `SKILL.md` runs `check_env.sh` through dynamic context injection, resolved from the skill's own directory and guarded: `` !`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" <tools…> || true` `` (the router uses `${CLAUDE_SKILL_DIR}/scripts/…`). Never reference scripts through `${CLAUDE_PROJECT_DIR}`: a personal install would then execute whatever the opened project ships.
- `scripts/lint_skills.py`, `scripts/test_scripts.sh`, `.github/workflows/lint.yml` — the only checks. No build, no package manager, no runtime dependencies beyond bash, python3, shellcheck.

## Rules for editing skills

1. Keep the contract: frontmatter `name` equals the directory name (`[a-z0-9-]`), `description` says what it does **and** when to use it, `paths:` lists only distinctive globs (generic `*.bam`/`*.vcf*` belong to `genomics` only). Never add `allowed-tools`, `model`, `context`, or `hooks` to frontmatter.
2. Section order in every sub-skill: title + role, `## Operating rules` (verbatim block referencing the router's eight rules), `## Environment`, `## Inspect first`, `## Decision table`, `## Workflow`, `## QC thresholds`, `## Pitfalls`, `## nf-core`, `## Validation`, `## Additional resources`.
3. `SKILL.md` ≤ 400 lines (CI fails at 500); supporting files ≤ 600. Move long command dumps out of `SKILL.md`.
4. Bash blocks must parse: no comment or blank line inside a `\` continuation, every fence has a language, `${UPPER_CASE}` placeholders declared at the top of the block, multi-step blocks start with `set -euo pipefail`, aligner pipes end with `samtools quickcheck`.
5. Versions and URLs: cite only what you verified against a primary source today; otherwise write `<check current: URL>`. Never invent database dates. Update `## Verified against` in the README and `CHANGELOG.md` when you bump pins.
6. One owner per topic (no duplication): raw-read QC → `genomics-qc`; Nextflow/Snakemake/SLURM mechanics → `genomics-workflow`; reference/database downloads and SRA/ENA/GEO retrieval → `genomics/reference.md`; truth sets and hap.py/Truvari → `genomics-benchmark`; somatic SNV/indel → `genomics-somatic`; long-read callers → `genomics-longread`. Link by path (`.claude/skills/genomics-qc/SKILL.md`) instead of repeating.
7. Deprecated tools are lint errors inside code blocks: MACS2, lumpyexpress, GEMINI, R `copynumber`, `nf-core/copynumber`, `nf-core/fastqc`, `samtools view -bS`, `INFO/FILTER`, VEP `--plugin ClinVar`, `--plugin gnomADc`, standalone HMMRATAC, DeepVariant < 1.10. Mention them only in prose as "legacy".
8. Clinical content (ACMG, CNV scoring, somatic tiers) carries the verbatim "Research use only" block.
9. Scripts stay dependency-free bash, `shellcheck -S warning` clean, never write outside `mktemp -d`, and `check_env.sh` always exits 0 (it runs at skill load).

## Verify before you finish

```bash
python3 scripts/lint_skills.py     # 0 errors required; --offline skips the nf-core registry check
bash scripts/test_scripts.sh       # "scripts OK"
```

Both run in CI on push/PR and monthly (nf-core renames, dead links). Do not add other CI steps that need network secrets; the repo has none and should stay that way.

## Adding a skill

1. `mkdir .claude/skills/genomics-<topic>` and write `SKILL.md` following rules 1–8; copy the operating-rules block from any sibling.
2. Add a routing row in `.claude/skills/genomics/SKILL.md` and a row in the README `Skills` table and `Layout` tree.
3. Add a `CHANGELOG.md` entry; run the two checks above.

## Do not

- Commit data files, reference genomes, or anything under `results/`, `work/`, `.nextflow*`, `.snakemake/`.
- Add `settings.json` skill registration; Claude Code discovers `.claude/skills/*/SKILL.md` on its own.
- Weaken the operating rules in `genomics/SKILL.md` to make a workflow shorter; they exist because the consumer is an autonomous agent handling controlled-access data.
