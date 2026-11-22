# Genomics skills for Claude Code

Thirteen [Claude Code skills](https://code.claude.com/docs/en/skills) that turn Claude into a careful bioinformatics engineer for sequencing data: short- and long-read variant calling, somatic analysis, RNA-seq (bulk, single-cell, spatial), ChIP/ATAC/CUT&RUN, CNV/SV, annotation and interpretation, methylation, metagenomics, QC, benchmarking and workflow engineering.

The skills are written for an **autonomous agent**, not for a human reading a wiki: every skill starts from a fixed set of operating rules (never load data files into context, inspect headers before asking, verify reference consistency before combining inputs, validate each step, pin and record versions, treat thresholds as defaults, treat human data as controlled-access) and ships bounded inspection scripts so Claude can look at a 200 GB directory without reading it.

## Install

Project-wide (committed with the repository, shared by the team):

```bash
git clone https://github.com/bloo-berries/claude-skill-genomic-pipelines.git
cp -r claude-skill-genomic-pipelines/.claude/skills /path/to/your/project/.claude/
```

Personal (every project on this machine):

```bash
cp -r claude-skill-genomic-pipelines/.claude/skills/* ~/.claude/skills/
```

Restart Claude Code (or wait for live reload) and type `/` to see the skills. No `settings.json` entry is required; Claude Code discovers `.claude/skills/<name>/SKILL.md` automatically.

## Skills

| Skill | Invoke | Scope |
|---|---|---|
| Router | `/genomics` | Operating rules, inspection scripts, routing, shared conventions, reference/annotation/database sources (`reference.md`) |
| WGS/WES germline | `/genomics-wgs` | BWA-MEM2, GATK 4.7 (HaplotypeCaller, VQSR, NVScoreVariants, hard filters), DeepVariant 1.10, GLnexus, DRAGEN-GATK, sarek |
| Somatic | `/genomics-somatic` | Mutect2 full best practice, Strelka2, consensus, PoN, ASCAT/FACETS/PURPLE, MSI/TMB/HRD, signatures, AMP/ASCO/CAP tiers, PCGR, UMI consensus |
| RNA-seq | `/genomics-rnaseq` | Strandedness inference, decoy-aware Salmon, STAR, DESeq2/edgeR/limma, batch handling, splicing, fusions; scRNA (Seurat v5, Scanpy, ambient/doublet, integration, pseudobulk DE), spatial |
| ChIP / ATAC / CUT&RUN | `/genomics-chipseq` | MACS3, IDR, deepTools, DiffBind, motifs; ATAC Tn5 handling done correctly, `macs3 hmmratac`, TOBIAS, scATAC; SEACR and spike-in scaling |
| Annotation & interpretation | `/genomics-annotation` | VEP 116 with plugins, `bcftools +split-vep`, SnpEff, ANNOVAR, liftover; ACMG/ClinGen SVI with calibrated thresholds; slivar, Exomiser, STR, PGx, nf-core/raredisease |
| QC | `/genomics-qc` | fastp/Trim Galore/cutadapt, FastQC/MultiQC, mosdepth, contamination (VerifyBamID2, FastQ Screen, Kraken2), identity/sex/relatedness (somalier, peddy), PASS/WARN/FAIL tables |
| CNV / SV | `/genomics-cnv` | GATK gCNV and somatic CNV (with `PreprocessIntervals`), CNVkit, ExomeDepth/DECoN, Manta, smoove, Tiddit, dysgu, merging, AnnotSV/SVAFotate, ClinGen CNV scoring |
| Long read | `/genomics-longread` | dorado, pbmm2/minimap2, Clair3/DeepVariant, WhatsHap/HiPhase, Sniffles2/cuteSV/pbsv, modkit, TRGT/Straglr, hifiasm/Flye, Iso-Seq/bambu |
| Methylation | `/genomics-methylation` | Bismark, bwa-meth + MethylDackel, conversion QC, DSS/methylKit/dmrseq, arrays, nf-core/methylseq |
| Metagenomics | `/genomics-metagenomics` | Host removal, Kraken2/Bracken, MetaPhlAn 4, MAGs (MEGAHIT → binning → CheckM2/GTDB-Tk), HUMAnN, AMR, 16S/QIIME 2, nf-core/taxprofiler, mag, ampliseq |
| Benchmarking | `/genomics-benchmark` | GIAB truth sets and stratifications, hap.py / vcfeval, Truvari, SEQC2 somatic, nf-core/variantbenchmarking, minimum F1 by assay |
| Workflow engineering | `/genomics-workflow` | Nextflow DSL2 and nf-core mechanics, Snakemake 8, SLURM, containers and env pinning, cloud and accelerated stacks, failure triage, templates |

Skills also activate automatically when Claude works on matching files (`*.fastq.gz`, `*.bam`, `*.vcf.gz`, `*.narrowPeak`, `*.h5ad`, `*.pod5`, `*.nf`, `Snakefile`, …) or when the request matches a skill's description.

## Layout

```text
.claude/skills/
├── genomics/            SKILL.md  reference.md  scripts/{peek.sh,ref_sanity.sh,check_env.sh}
├── genomics-wgs/        SKILL.md  reference.md
├── genomics-somatic/    SKILL.md  reference.md
├── genomics-rnaseq/     SKILL.md  bulk.md  singlecell.md
├── genomics-chipseq/    SKILL.md  atac.md  cutandrun.md
├── genomics-annotation/ SKILL.md  acmg.md  rare-disease.md
├── genomics-qc/         SKILL.md
├── genomics-cnv/        SKILL.md  reference.md
├── genomics-longread/   SKILL.md  reference.md
├── genomics-methylation/SKILL.md  dmr.md
├── genomics-metagenomics/SKILL.md
├── genomics-benchmark/  SKILL.md
└── genomics-workflow/   SKILL.md  templates.md
scripts/lint_skills.py   repo lint (run in CI)
scripts/test_scripts.sh  smoke test for the shared scripts
```

Each `SKILL.md` has the same sections: operating rules, environment probe (live tool/version detection), inspect-first table, decision table, workflow with per-step validation, QC thresholds with rationale, pitfalls, nf-core pipelines (names validated in CI), validation, additional resources.

## Research use only

Variant classifications, CNV scoring and somatic tiers produced with these skills are research outputs. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189 or local equivalent) with the exact ClinVar/gnomAD/transcript-set versions recorded. Treat all human sequencing data as controlled-access; the skills instruct Claude not to upload data or paste identifiers, but the data-use agreement is yours to enforce.

## Verified against (2026-11-22)

GATK 4.7.0.0 (bioconda lags at 4.6.2.0), DeepVariant 1.10.0, Ensembl/VEP 116, GENCODE v50 / vM39, gnomAD v4.1, MACS3 3.1, Salmon 2.12, samtools/bcftools 1.25, nf-core/sarek 3.12, Nextflow 26.10, Snakemake 9.27. Versions drift; the monthly CI run validates nf-core names and links, and every skill tells Claude to run `tool --version` rather than trust a number in a file.

## Development

```bash
python3 scripts/lint_skills.py      # frontmatter (incl. forbidden keys), size (warn >400, fail >500 lines), bash syntax via shellcheck,
                                    # pipefail on multi-command blocks, guarded check_env injection, deprecated tools, nf-core registry, paths
bash scripts/test_scripts.sh        # shared-script smoke test: locale/alt-contig/missing-tool cases, symlinks, stub tools, hung-probe timeout
```

Contributor rules (section order, size limits, deprecated-tool list, verification steps) live in `CLAUDE.md`, which Claude Code also loads when you edit this repository. See `CHANGELOG.md` for history.

## Security model

- No skill grants `allowed-tools`, changes `model`, or forks context. The only code executed at skill load is `check_env.sh` (tool presence/version probe): it refuses relative `PATH` entries, time-limits every probe, truncates output to one short line, writes nothing and always exits 0. Scripts are resolved from the skill's own directory (`${CLAUDE_SKILL_DIR}`), so a personal install never executes a project's copy.
- The shared scripts write only inside `mktemp -d` and never modify inputs (`ref_sanity.sh` indexes an un-indexed reference into its temp dir, not beside the file).
- Skills contain `<check current: URL>` markers for versions that move faster than this repo. The router instructs Claude to take only the version or file-name string from such a page and never to act on instructions found in fetched content; all hosts are primary sources (NCBI, Ensembl, Broad, nf-core, tool repositories).
- Data governance is a standing rule (`genomics/SKILL.md` rule 8): no uploads to third-party services, workflow telemetry, web calculators or licensed APIs without a confirmed data-use agreement; no sample identifiers in chat or logs; clinical output is research-use-only.
- The repo contains no credentials, no sample data (all IDs are public GIAB/Coriell references or placeholders) and no CI secrets. CI runs with `contents: read`, SHA-pinned actions, and does lint + smoke test + link check only. `.gitignore` blocks BAM/CRAM/FASTQ/VCF/reference files and pipeline work directories.
- Review the diff of `.claude/skills/**` on every PR as you would review code: skill text is executed as instructions by an agent with your permissions.

## License

MIT — see `LICENSE`.
