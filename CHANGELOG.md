# Changelog

## 2.0.1 — 2026-11-22

Security and adversarial review of 2.0.0 (two independent read-only review passes; all findings addressed).

### Security
- Scripts are resolved from `${CLAUDE_SKILL_DIR}` instead of `${CLAUDE_PROJECT_DIR}`, and the dynamic-context line is guarded (`[ -f "$f" ] && … || true`): a personal install no longer executes whatever the opened project ships, and a missing script no longer injects error text. Lint rejects `CLAUDE_PROJECT_DIR` and any non-guarded injection line.
- `check_env.sh`: refuses relative `PATH` entries, time-limits every probe (10 s), truncates output to one line, rejects tool names with shell metacharacters.
- `ref_sanity.sh`: never writes beside the reference (indexes into its temp dir); contig style is a majority vote over primary contigs under `LC_ALL=C` (was: first sorted name, which misclassified Broad references with `HLA-*` contigs); alt/HLA/decoy-only gaps are `WARN`, missing tools are `SKIP` (exit 2) instead of a false `MISMATCH` from binary bytes; BED `track`/`browser` headers ignored.
- `peek.sh`: symlink-safe sizes, no arithmetic error on <4-line FASTQ, never dumps BCF bytes when bcftools is absent.
- SLURM array template: samplesheet fields passed as arguments (no shell interpolation of sample IDs/paths), sample ID whitelisted, per-job `mktemp -d` scratch so the `EXIT` trap cannot delete a shared `SLURM_TMPDIR`.
- CI: `permissions: contents: read`, actions pinned to commit SHAs, `persist-credentials: false`, lychee `--exclude-private --exclude-link-local --exclude-mail`; lint forbids `allowed-tools`/`model`/`context`/`hooks`/`effort` frontmatter.
- Router rule 6 narrowed to primary-data destruction (regenerating indexes/reports is allowed); rule 8 now covers workflow telemetry, web calculators and licensed APIs, and tells Claude to take only version strings from `<check current: URL>` pages. ChIP/CUT&RUN/methylation default to marking duplicates and excluding by flag instead of physical removal; Seqera `-with-tower`, Terra and the ClinGen CNV web calculator carry DUA caveats.

### Fixed
- ATAC chrM filter used a non-existent htslib identifier (`mrnm` → `mrname`); the step failed under `set -e`.
- RNA-seq strandedness block did not persist `KALLISTO_STRAND`/`RMATS_LIB` into `strand.env`; kallisto and rMATS blocks aborted with unbound variables.
- WGS: `REF_MISMATCH` validation grepped the wrong file (stderr is now captured to `.norm.log`); somalier sex check printed the wrong columns and lacked `--infer`; `verifybamid2` vs `VerifyBamID` binary name unified; WES/WGS interval handling made an explicit `MODE` switch; HaplotypeCaller scatter uses `xargs -P` with exit-code propagation instead of 24 unchecked background jobs.
- QC: somalier sex-mismatch awk compared sex to phenotype (`$6` → `$7`).
- Annotation: VEP field-count validation anchored (`grep -x`) so dbNSFP `REVEL`/`LOEUF` columns are not double-counted.
- Benchmark: hap.py summary column indices (`11-14`). Methylation: Bowtie2 `.bt2` vs `.bt2l` index check, GENCODE v50 and gene-name extraction in `dmr.md`. ChIP description referenced a non-existent `genomics-methyl` skill.
- Missing `set -euo pipefail` in six multi-command blocks; missing `mkdir -p` before redirects in four blocks. Lint now enforces both the pipefail rule and deprecated patterns across backslash-continued lines; fence parser tolerates info strings.
- Cross-skill threshold drift (duplicate rate, FREEMIX, WES mapping %, WES Ti/Tv) aligned to `genomics-qc` as the single owner; `genomics/reference.md` gained the sections sibling skills already pointed to (spike-in genomes, mappability/segdup tracks, PAR/TR/exclusion BEDs, dosage maps, CADD/dbNSFP/chain files, metagenomics databases, bisulfite controls); `nf-core/smrnaseq` row added to rnaseq; router resource table no longer duplicates wgs numbers.

### Added
- `CLAUDE.md` (contributor rules for agents editing this repo), `.gitignore` for data/work directories, README "Security model" section, expanded script smoke tests (locale, alt contigs, missing tools, symlinks, stub tools, hung-probe timeout).

## 2.0.0 — 2026-11-22

Complete rewrite. The 1.x files were never loaded by Claude Code (flat `.claude/skills/*.md` plus an invented `settings.json` `skills` key); this release moves to the supported `.claude/skills/<name>/SKILL.md` layout with YAML frontmatter.

### Structure
- One directory per skill, `SKILL.md` ≤ 400 lines, long material in supporting files (`reference.md`, `bulk.md`, `atac.md`, …).
- Shared operating rules in the `genomics` router (no data files in context, inspect-before-ask, reference sanity gate, `pipefail` + per-step validation, pin and record versions, estimate before spending, thresholds are defaults, controlled-access data).
- Shared scripts: `peek.sh` (bounded preview of FASTQ/BAM/CRAM/VCF/BED/GTF), `ref_sanity.sh` (contig name/length/build consistency gate), `check_env.sh` (tool/version/scheduler probe used via dynamic context injection).
- `paths:` frontmatter so skills activate on relevant files; `genomics` owns the generic globs.
- CI: `scripts/lint_skills.py` (frontmatter, size, bash-block syntax via shellcheck, deprecated tools, nf-core registry, internal paths), `scripts/test_scripts.sh`, lychee link check. Removed `settings.json`; added `LICENSE` (MIT, as the README already claimed).

### New skills
`genomics-longread`, `genomics-somatic`, `genomics-methylation`, `genomics-metagenomics`, `genomics-benchmark`, `genomics-workflow`.

### Fixed (commands that were broken or wrong)
- VEP command broken by comments inside backslash continuation; non-existent `--plugin ClinVar`; `--plugin gnomADc` misdescribed as constraint; filtering on `INFO/SYMBOL`/`INFO/gnomAD_AF` (fields live in CSQ → `bcftools +split-vep`); `INFO/FILTER`.
- ATAC-seq: MACS `-f BAMPE` combined with `--shift/--extsize` on already Tn5-shifted reads; gawk-only `and()` read shift that broke mate coordinates; `grep -v chrM` on SAM text; standalone HMMRATAC; `TSSEscore()` given a path.
- GATK gCNV/somatic CNV missing `PreprocessIntervals`; `FilterIntervals --minimum-mappability` without a mappability track.
- Mutect2 shown without orientation-bias model, contamination estimation or PoN construction.
- DeepVariant run on BQSR output; `featureCounts` without `-s`; Salmon index without decoys; `lfcShrink(coef=)` without `relevel()`; DiffBind `minMembers` trap.
- `mosdepth -b … --by 1`; `multiqc -m featureCounts`; `verifybamid2 --Ref`; `nf-core/fastqc` and `nf-core/copynumber` (do not exist); R `copynumber` (removed from Bioconductor); `wget` HTTPS glob; NCBI dbSNP contig names; sarek `--genome` inconsistency and `--wes` without `--intervals`.

### Updated
DeepVariant 1.10, VEP/Ensembl 116, GENCODE v50 / vM39, gnomAD v4.1 (730,947 exomes + 76,215 genomes), MACS3, GATK 4.7 with NVScoreVariants, JASPAR2024, GRCm39, ClinGen SVI calibrated PP3/BP4 thresholds (REVEL, CADD, SpliceAI, AlphaMissense), ACMG without PP5/BP6, smoove for LUMPY, slivar for GEMINI, ENCODE TSS-enrichment bands.

## 1.0.1 — 2026-01-22
- Fix reference file URLs to use verified public sources.

## 1.0.0 — 2026-01-22
- Initial release: genomics, wgs, rnaseq, chipseq, annotation, qc, cnv.
