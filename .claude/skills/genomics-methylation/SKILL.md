---
name: genomics-methylation
description: DNA methylation analysis from bisulfite/enzymatic conversion sequencing (WGBS, RRBS, EM-seq, PBAT, targeted) and brief EPIC array handling; alignment with Bismark or bwa-meth, M-bias-aware extraction, conversion-rate QC, DMR calling and nf-core/methylseq. Use when inputs are BS-converted FASTQs, Bismark/MethylDackel outputs (`*.cov.gz`, `*.bedGraph`, `CpG_report.txt`), IDATs, or the user says WGBS, RRBS, EM-seq, methylseq, DMR, 5mC. Not for native long-read 5mC/6mA calls (`genomics-longread`) or ChIP/ATAC (`genomics-chipseq`).
argument-hint: "[wgbs|rrbs|emseq|targeted|array|dmr] [path]"
paths:
  - "**/*.bismark.cov.gz"
  - "**/*CpG_report.txt*"
  - "**/*_CpG.bedGraph"
  - "**/*.idat"
---

# Genomics: DNA methylation

Expert in conversion-based methylation sequencing: library-type detection, three-letter alignment (Bismark, bwa-meth/MethylDackel), deduplication rules, M-bias trimming, conversion-efficiency QC, per-CpG coverage files, and group-level DMR/DML statistics. Outputs are commands, decision tables and validation checks, not tutorials.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" bismark bowtie2 bwameth.py bwa-mem2 MethylDackel trim_galore samtools nextflow Rscript metilene || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

Bismark ≥3.0 is a single Rust binary (`bismark align|dedup|extract|…`); the classic names (`bismark`, `deduplicate_bismark`, `bismark_methylation_extractor`, `bismark_genome_preparation`) remain supported aliases and are used below. Perl Bismark v0.25.x is legacy (`legacy_perl/`).

## Inspect first

| Check | Command | Infers |
| --- | --- | --- |
| Read composition | `peek.sh R1.fastq.gz 400 \| awk 'NR%4==2' \| grep -o '[ACGT]' \| sort \| uniq -c` | Converted library if C ≪ G on R1 (directional); C-poor **and** G-poor on both reads suggests PBAT/non-directional; balanced bases = not converted (wrong file or TAPS/5-Base). |
| RRBS signature | `peek.sh R1.fastq.gz 400 \| awk 'NR%4==2{print substr($1,1,3)}' \| sort \| uniq -c \| sort -rn \| head -3` | ≥50 % reads starting `CGG`/`TGG` = MspI RRBS → `--rrbs`, no dedup. |
| Read length / layout | `peek.sh R1.fastq.gz 4`; presence of `_R2`/`_2` file | SE vs PE; 2×150 WGBS vs 2×75 RRBS; EM-seq usually 2×150. |
| Existing BAM | `samtools view -H x.bam \| grep -E '^@(PG\|SQ)' \| head` | Bismark (`@PG ID:Bismark`, `XM:Z` tags) vs bwa-meth (`@PG ID:bwa`, `YD:Z` tag); contig naming (`chr1` vs `1`); assembly by `LN` of chr1 (GRCh38 248956422, GRCm39 195154279, mm10/GRCm38 195471971). |
| Spike-ins in reference | `grep -c -E '^>(lambda\|J02459\|pUC19\|L09137)' ref.fa` | Whether conversion controls can be measured directly. |
| Coverage files | `peek.sh x.cov.gz 5` | Bismark cov (`chr start end %meth countM countU`) vs MethylDackel bedGraph (`chr start end %meth countM countU` after 1-line track header) vs cytosine report (7 cols with strand+context). |
| Array | `ls *.idat \| head; ls *_Grn.idat \| wc -l` | EPIC v2 (~935k probes, 2023+) vs EPIC v1 vs 450k — determined from the manifest SeSAMe infers, not from filenames. |

Ask only if still ambiguous: library kit (directional vs PBAT/non-directional), whether lambda/pUC19 spike-in was added, target BED for capture panels, which groups to contrast, and whether the reference assembly must match a prior cohort (e.g. GRCm39 vs legacy mm10).

## Decision table

| Situation | Choose | Why |
| --- | --- | --- |
| WGBS/EM-seq, human/mouse, ≥30× | `--aligner bwameth` (+ `--use_mem2`) and MethylDackel | Faster, standard BAM, MAPQ usable; EM-seq is directional and bisulfite-equivalent downstream. |
| RRBS, non-directional, PBAT, single-cell, legacy comparability | Bismark (bowtie2) | Native handling of `--non_directional`, `--pbat`, RRBS no-dedup conventions, richer reports. |
| Directionality | `--non_directional` only for kits documented as such (e.g. some PBAT/post-bisulfite, Zymo Pico, scBS) | Default directional discards CTOT/CTOB reads; wrong setting halves mapping or inflates ambiguity. |
| Deduplication | WGBS/EM-seq: yes (`deduplicate_bismark` or `samtools markdup` for bwa-meth); RRBS/amplicon/capture: **no** | Restriction/amplicon fragments start at identical positions by design; dedup removes real coverage. |
| Read trimming | Trim Galore `--rrbs` for MspI RRBS; EM-seq/PBAT/Pico: `--clip_R1 10 --clip_R2 10` (then verify with M-bias) | Removes end-repair fill-in (RRBS 3' `CG`) and random-primer bias. Raw-read QC otherwise lives in `.claude/skills/genomics-qc/SKILL.md`. |
| Conversion QC | Lambda/pUC19 spike-in if present; else non-CpG (CHH) methylation on autosomes as proxy | Mammalian somatic CHH ≈ 0.3–1 %; higher = incomplete conversion (or ESC/brain biology). |
| DMR caller | 2–3 reps/group, WGBS: DSS (smoothing) or bsseq/dmrseq; RRBS/targeted, many samples/covariates: methylKit (logistic) or DSS with design; fast no-R: metilene | dmrseq controls FDR at region level; methylKit handles covariates/overdispersion; DSS tolerates low replication. |
| Array (EPIC v2/v1/450k) | SeSAMe (`openSesame`) for QC+normalisation; minfi/`minfiData`/`ChAMP` for legacy pipelines | SeSAMe handles EPIC v2 replicate probes and pOOBAH masking; keep array work brief — this skill is sequencing-first. |
| Native 5mC/5hmC on ONT/PacBio | Do not use this skill's pipeline; see `.claude/skills/genomics-longread/SKILL.md` (modkit, bedMethyl) | No conversion; methylation comes from basecaller mod tags. |
| cfDNA methylation / fragmentomics | Align as EM-seq (bwa-meth), keep duplicates marked not removed, retain fragment length (`TLEN`) | Tissue-of-origin deconvolution needs fragment-level calls; dedup+UMI handling kit-specific. |
| TAPS / Illumina 5-Base | nf-core/methylseq `--aligner bwamem --taps` (rastair) or Bismark `--taps`/`--illumina_5base` | Chemistry inverts conversion (5mC→T); bisulfite aligners mis-call. |

## Workflow

Variables used throughout: `REF` (FASTA with spike-ins appended), `SAMPLE`, `THREADS`.

1. Prepare the reference once (append lambda/pUC19 if the kit used them, then index for both aligners).

```bash
set -euo pipefail
REF=ref/GRCh38_plus_lambda_pUC19.fa
THREADS=16
samtools faidx "${REF}"
bismark_genome_preparation --bowtie2 --parallel 4 "$(dirname "${REF}")"
bwameth.py index-mem2 "${REF}"
```
Validate: `ls "$(dirname "${REF}")"/Bisulfite_Genome/CT_conversion/*.bt2 "$(dirname "${REF}")"/Bisulfite_Genome/CT_conversion/*.bt2l 2>/dev/null | head -n 1 && ls "${REF}".bwameth.c2t*` (bowtie2 writes `.bt2l` only for genomes > 4 Gbp; human/mouse get `.bt2`).

2. Trim by library type (adapter/quality defaults are Trim Galore's; this step only sets conversion-specific clipping).

```bash
set -euo pipefail
SAMPLE=S1
LIBTYPE=wgbs
case "${LIBTYPE}" in
  rrbs)  EXTRA="--rrbs" ;;
  emseq|pbat|pico) EXTRA="--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10" ;;
  *)     EXTRA="" ;;
esac
trim_galore --paired --cores 4 ${EXTRA} -o trimmed "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz"
```
Validate: `grep -E 'Reads written|Total written' trimmed/${SAMPLE}_R1.fastq.gz_trimming_report.txt`

3a. Align with Bismark (RRBS/non-directional/PBAT branch). Add `--non_directional` only when the kit is documented non-directional; `--pbat` for PBAT.

```bash
set -euo pipefail
SAMPLE=S1
GENOME_DIR=ref
THREADS=16
bismark --genome "${GENOME_DIR}" --parallel 4 -p 4 -1 "trimmed/${SAMPLE}_R1_val_1.fq.gz" -2 "trimmed/${SAMPLE}_R2_val_2.fq.gz" -o align --temp_dir tmp
samtools quickcheck "align/${SAMPLE}_R1_val_1_bismark_bt2_pe.bam"
```
Validate: `grep -E 'Mapping efficiency|C methylated in CHH context' align/${SAMPLE}_R1_val_1_bismark_bt2_PE_report.txt`

3b. Align with bwa-meth (WGBS/EM-seq branch).

```bash
set -euo pipefail
SAMPLE=S1
REF=ref/GRCh38_plus_lambda_pUC19.fa
THREADS=16
bwameth.py --threads "${THREADS}" --reference "${REF}" --read-group "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "trimmed/${SAMPLE}_R1_val_1.fq.gz" "trimmed/${SAMPLE}_R2_val_2.fq.gz" | samtools sort -@ 8 -m 2G -o "align/${SAMPLE}.bam" - && samtools quickcheck "align/${SAMPLE}.bam"
samtools index "align/${SAMPLE}.bam"
```
Validate: `samtools flagstat align/${SAMPLE}.bam | grep -E 'mapped \(|properly paired'`

4. Deduplicate — skip for RRBS, amplicon, hybrid-capture.

```bash
set -euo pipefail
SAMPLE=S1
deduplicate_bismark --paired --bam -o "${SAMPLE}" --output_dir dedup "align/${SAMPLE}_R1_val_1_bismark_bt2_pe.bam"
samtools sort -@ 8 -o "dedup/${SAMPLE}.sorted.bam" "dedup/${SAMPLE}.deduplicated.bam" && samtools index "dedup/${SAMPLE}.sorted.bam"
```
For bwa-meth BAMs: `samtools markdup -@ 8 align/${SAMPLE}.bam dedup/${SAMPLE}.bam` (after `samtools fixmate -m` on a name-sorted BAM) marks duplicates without removing them; downstream MethylDackel calls exclude them with `--ignoreFlags 1024` (any `samtools view` on the BAM: `-F 1024`). Validate: `grep 'Total number duplicated' dedup/${SAMPLE}.deduplication_report.txt` or `samtools flagstat dedup/${SAMPLE}.bam | grep duplicates`.

5. Measure M-bias, then extract with trimming derived from it (do not guess trim lengths).

```bash
set -euo pipefail
SAMPLE=S1
REF=ref/GRCh38_plus_lambda_pUC19.fa
mkdir -p qc
MethylDackel mbias -@ 8 --ignoreFlags 1024 "${REF}" "dedup/${SAMPLE}.bam" "qc/${SAMPLE}" > "qc/${SAMPLE}.mbias_suggest.txt" 2>&1
grep -- '--OT' "qc/${SAMPLE}.mbias_suggest.txt"
```
`mbias` prints a `Suggested inclusion options: --OT a,b,c,d --OB ...` line; pass exactly that to `extract`:

```bash
set -euo pipefail
SAMPLE=S1
REF=ref/GRCh38_plus_lambda_pUC19.fa
INCL="$(grep -o -- '--OT.*' "qc/${SAMPLE}.mbias_suggest.txt")"
MethylDackel extract -@ 8 --ignoreFlags 1024 --mergeContext --minDepth 1 ${INCL} "${REF}" "dedup/${SAMPLE}.bam" -o "meth/${SAMPLE}"
```
Bismark branch (M-bias comes from the extractor report; rerun with `--ignore_r2 N`/`--ignore_3prime N` if the plot shows end bias):

```bash
set -euo pipefail
SAMPLE=S1
GENOME_DIR=ref
bismark_methylation_extractor --paired-end --gzip --parallel 4 --bedGraph --cytosine_report --genome_folder "${GENOME_DIR}" --no_overlap -o meth "dedup/${SAMPLE}.deduplicated.bam"
```
Validate: `peek.sh meth/${SAMPLE}_CpG.bedGraph 3` (MethylDackel) or `zcat meth/${SAMPLE}.deduplicated.bismark.cov.gz | head -3`; `bismark2report` / `bismark2summary` for HTML.

6. Conversion-rate QC. Spike-in: methylation on lambda (should be ≈0 %) and pUC19 (≈100 %, EM-seq kits):

```bash
set -euo pipefail
SAMPLE=S1
for CTG in lambda pUC19; do
  awk -v c="${CTG}" '$1==c{m+=$5;u+=$6} END{if(m+u>0) printf "%s conv=%.3f%% (n=%d)\n", c, 100*u/(m+u), m+u; else print c, "no coverage"}' "meth/${SAMPLE}_CpG.bedGraph"
done
```
No spike-in: CHH proxy from Bismark report (`C methylated in CHH context`) or `MethylDackel extract --CHH --noCpG`. Validate: lambda conversion ≥99 %, pUC19 methylation ≥95 %, CHH ≤1.5 % (see thresholds).

7. Filter and prepare for DMR calling: drop CpGs overlapping known C>T/G>A SNPs (dbSNP or sample VCF), chrM, and sites outside the depth window; convert MethylDackel bedGraph to Bismark cov (`--methylKit` output option) if the caller needs it.

```bash
set -euo pipefail
SAMPLE=S1
SNPS=ref/dbsnp_CT_GA_only.bed.gz
tail -n +2 "meth/${SAMPLE}_CpG.bedGraph" | awk '$1!="chrM" && $1!="lambda" && $1!="pUC19"' | bedtools intersect -v -a - -b "${SNPS}" > "meth/${SAMPLE}.filtered.bedGraph"
```
Validate: `wc -l meth/${SAMPLE}_CpG.bedGraph meth/${SAMPLE}.filtered.bedGraph` (expect ≤5 % loss; more suggests contig mismatch — rerun `ref_sanity.sh`).

8. DMR calling (DSS shown; see `dmr.md` for methylKit, dmrseq, metilene and covariate designs).

```r
suppressPackageStartupMessages({library(DSS); library(bsseq)})
files  <- c(ctrl1="meth/C1.cov", ctrl2="meth/C2.cov", case1="meth/T1.cov", case2="meth/T2.cov")
read1  <- function(f) { d <- read.table(f); data.frame(chr=d$V1, pos=d$V2, N=d$V5+d$V6, X=d$V5) }
bs     <- makeBSseqData(lapply(files, read1), names(files))
dml    <- DMLtest(bs, group1=c("ctrl1","ctrl2"), group2=c("case1","case2"), smoothing=TRUE, ncores=8)
dmls   <- callDML(dml, p.threshold=0.001)
dmrs   <- callDMR(dml, delta=0.1, p.threshold=0.01, minlen=50, minCG=3, dis.merge=100)
write.table(dmrs, "dmr/DSS_dmrs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
sessionInfo()
```
Validate: `wc -l dmr/DSS_dmrs.tsv`; check `areaStat` sign distribution is not 100 % one-directional unless biology predicts it; plot 3 top DMRs with `showOneDMR`.

## QC thresholds

| Metric | Default | Rationale / when to deviate |
| --- | --- | --- |
| Mapping efficiency (Bismark unique) | WGBS ≥60 %, RRBS ≥50 %, EM-seq ≥70 % | Three-letter genome lowers uniqueness; <40 % → wrong directionality, untrimmed adapters, or wrong genome. |
| Conversion (lambda) | ≥99.0 % (EM-seq typically ≥99.5 %) | Below → under-converted library; every CpG inflates. Report but do not "correct" by subtraction without a stated model. |
| pUC19 methylated control | ≥95 % methylated | Low → over-conversion/enzymatic damage (EM-seq). |
| CHH methylation (proxy) | ≤1.5 % somatic mammalian tissue | ESCs, oocytes, brain and plants have real non-CpG methylation — do not treat as failure there. |
| Duplication (WGBS) | ≤20 % | Higher at low input/PCR cycles; for RRBS this metric is meaningless. |
| CpG depth for per-site calls | ≥10× per CpG (≥5× with smoothing/DSS) | Binomial noise; targeted/RRBS often use ≥10×, cap at ≤99.9th percentile to drop repeats. |
| Genome-wide CpG covered ≥10× | WGBS ≥70 % at 30× raw; RRBS 1–3 M CpGs | Lower → insufficient depth; RRBS count depends on size selection. |
| M-bias | Flat after trimming; drop positions deviating >±5 % from plateau | Untrimmed EM-seq/PBAT shows 5′ dips; `mbias` suggestions are per strand/read. |
| Insert size | WGBS mode 200–400 bp; cfDNA mode ~167 bp | Short inserts in WGBS imply degraded input; dovetailing requires `--no_overlap`. |
| DMR call | Δβ ≥0.1 (≥0.2 for RRBS/clinical), FDR ≤0.05, ≥3 CpGs, ≥50 bp | Smaller Δβ rarely validates by pyrosequencing; justify any relaxation by replicate count. |
| Array detection p / probe masking | SeSAMe pOOBAH p<0.05, drop samples with >10 % failed probes | Replace minfi detectionP for EPIC v2. |

## Pitfalls

| Symptom | Cause | Fix |
| --- | --- | --- |
| Mapping ≈ half expected; ambiguity high | Non-directional library aligned as directional (or vice versa) | Inspect base composition of R1 and R2; rerun with/without `--non_directional`/`--pbat`. |
| Duplication 60–90 % on RRBS | Dedup applied to restriction fragments | Never deduplicate RRBS/amplicon; use `samtools markdup` stats only as information. |
| Methylation ~0 % everywhere or ~50 % with balanced reads | Not a converted library (TAPS, 5-Base, or plain WGS) | Check composition; route to `--taps`/`--illumina_5base` or stop. |
| 5′ end of R2 shows 10–20 % lower methylation in M-bias | End-repair fill-in (EM-seq, PBAT, Pico) | `--clip_R2`/`--ignore_r2` by M-bias plot; re-extract. |
| Spurious hypomethylated CpGs cluster at polymorphic sites | C>T SNP read as unmethylated C | Filter CpGs overlapping dbSNP C/T and G/A SNVs or sample genotype; for mouse strains use strain-specific SNPs. |
| Coordinates off by one / duplicated CpGs on both strands | Mixing `--mergeContext` (destranded) and per-strand files | Pick one convention per project; Bismark `coverage2cytosine --merge_CpG`. |
| Mouse DMRs do not lift to published coordinates | Reference is GRCm39 but cohort used legacy mm10 (GRCm38) | Record assembly; liftover regions with UCSC chain or re-align; never mix. |
| `ref_sanity.sh` fails on cov files | Spike-in contigs or `chr` prefix differ between reference and annotation | Use the same FASTA (with spike-ins) for alignment, extraction and annotation. |
| DSS/bsseq runs out of memory | Loading all CpGs across many samples | Filter to ≥5× sites, run per chromosome, or use `HDF5Array` backend. |
| EM-seq reported conversion 97 % by CHH but lambda 99.6 % | Real non-CpG methylation or bias in CHH proxy | Trust spike-in; CHH is a proxy only. |
| CpG islands under-covered in Illumina 5-Base | Methylated CpG-dense reads fail to align (documented Bismark limitation) | Use `--bwamem3` engine (Bismark ≥3.3) and report the caveat. |

## nf-core

`nf-core/methylseq` 4.2.0 (exists; `--aligner bismark|bismark_hisat|bwameth|bwamem`; library flags `--rrbs`, `--em_seq`, `--pbat`, `--non_directional`, `--single_cell`, `--zymo`, `--accel`, `--taps`; `--use_mem2`; `--cytosine_report`; `--run_targeted_sequencing --target_regions_file x.bed --collecthsmetrics`).

Samplesheet header (4.2.0; `genome` column optional):

```csv
sample,fastq_1,fastq_2,genome
S1,S1_R1.fastq.gz,S1_R2.fastq.gz,
```

```bash
set -euo pipefail
nextflow run nf-core/methylseq -r 4.2.0 -profile singularity --input samplesheet.csv --outdir results --fasta ref/GRCh38_plus_lambda_pUC19.fa --aligner bwameth --use_mem2 --em_seq --cytosine_report
```
Validate: `ls results/multiqc/*/multiqc_report.html results/methyldackel/*/*.bedGraph | head`. Nextflow/profile/resume mechanics: `.claude/skills/genomics-workflow/SKILL.md`. Pass `--rrbs` (disables dedup, sets Trim Galore `--rrbs`) or `--pbat`/`--non_directional` exactly as decided in the table above.

## Validation

- Spike-in controls are the primary truth: lambda ≈0 %, pUC19 ≈100 %; report both per sample.
- Global CpG methylation distribution bimodal (peaks near 0 and 80–100 %) for somatic mammalian WGBS; unimodal mid-range = failed conversion or mixed context.
- Imprinted DMRs (e.g. `H19/IGF2`, `SNRPN`, `KCNQ1OT1`) ≈50 % in normal somatic tissue; X-inactivation: female X promoters ≈50 %.
- Sex check: chrY CpG coverage near zero in XX samples; mismatch → sample swap.
- Technical replicates / EPIC-vs-WGBS overlap: Pearson r ≥0.9 on CpGs with ≥10×.
- Public benchmarks: NA12878/HG002 WGBS and EM-seq (ENCODE, NEB) for pipeline comparison; see `.claude/skills/genomics-benchmark/SKILL.md` for truth-set mechanics.
- Record: aligner + version, Bismark suite version (`bismark -v`), reference FASTA checksum, spike-in contigs, trim/M-bias parameters, DMR caller version and design formula (`sessionInfo()`).

## Additional resources

- `dmr.md` (this dir): methylKit, dmrseq, metilene and SeSAMe/EPIC v2 command blocks, covariate designs, annotation of DMRs — load when the user needs a caller other than DSS or array processing.
- `.claude/skills/genomics-qc/SKILL.md`: FastQC/fastp/Trim Galore adapter and quality defaults — load before deciding non-conversion-specific trimming.
- `.claude/skills/genomics-longread/SKILL.md`: ONT/PacBio native modification calling (modkit, bedMethyl merge with short-read calls).
- `.claude/skills/genomics-workflow/SKILL.md`: Nextflow profiles, `-resume`, cluster configs for nf-core/methylseq.
- `.claude/skills/genomics/reference.md`: reference FASTA, lambda/pUC19 sequence acquisition, dbSNP downloads.
- `.claude/skills/genomics-benchmark/SKILL.md`: comparing methylation calls against reference datasets.
