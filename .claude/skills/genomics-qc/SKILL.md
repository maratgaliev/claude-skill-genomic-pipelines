---
name: genomics-qc
description: Raw-read QC and trimming (FastQC/fastp/Trim Galore/cutadapt), alignment and coverage metrics (samtools, mosdepth, GATK/Picard metrics), contamination, sample identity and sex checks (VerifyBamID2, somalier, peddy, NGSCheckMate, FastQ Screen, Kraken2), MultiQC aggregation and PASS/WARN/FAIL reporting. Use when the user asks about FASTQ quality, adapters, trimming, duplication, coverage, FREEMIX, sample swaps, sex mismatch, MultiQC, or a QC report. Not for variant-level benchmarking (genomics-benchmark), RNA-seq gene-level QC (genomics-rnaseq), or long-read caller choice (genomics-longread).
argument-hint: "[reads|bam|identity|report] [path]"
paths:
  - "**/multiqc_config.y*ml"
  - "**/*_fastqc.zip"
  - "**/*.fastp.json"
  - "**/*.mosdepth.summary.txt"
  - "**/*.selfSM"
  - "**/*.somalier"
---

# Genomics QC

Expert in sequencing QC: raw reads, alignments, coverage, contamination, identity, sex, batch structure and report generation. Outputs are commands, decision tables, thresholds with rationale, and per-step validation. This skill is the only home for raw-read QC and trimming; sibling skills link here rather than repeating it.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" fastqc fastp trim_galore cutadapt multiqc samtools mosdepth gatk VerifyBamID somalier peddy fastq_screen kraken2 preseq NanoPlot cramino || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| Input | Read | Infer |
|---|---|---|
| FASTQ | `peek.sh R1.fastq.gz 8`; read-name header `@INSTRUMENT:RUN:FLOWCELL:LANE:...` | Platform (NovaSeq/NextSeq = 2-colour → poly-G), lane count, paired layout (`_R1`/`_R2`, `_1`/`_2`), UMI already in name (`:UMI` after read name), phred offset (modern = +33) |
| FASTQ size | `ls -lh`; `zcat R1.fastq.gz \| head -40000 \| awk 'NR%4==2{n++;s+=length($0)} END{print n, s/n}'` on a sample | Approx read length; estimate total reads from size ÷ (≈55 bytes/read at 150 bp gz) |
| BAM/CRAM | `samtools view -H x.bam \| grep -E '^@(SQ\|RG\|PG)' \| head -30` | Assembly (SQ lengths, chr prefix), RG count/platform/library, aligner and whether duplicates already marked (`@PG ID:MarkDuplicates`/`markdup`) |
| BAM pairs | `samtools view x.bam \| head -1000 \| cut -f2 \| sort \| uniq -c` | Paired vs single-end from FLAG bits |
| VCF | `bcftools view -h x.vcf.gz \| tail -3`; `bcftools query -l x.vcf.gz \| wc -l` | Caller, sample count, joint vs single-sample |
| Targets | `peek.sh target.bed`; `cut -f1 target.bed \| sort -u \| head` | WES/panel vs WGS; contig naming matches reference (`ref_sanity.sh`) |
| Expected sex / pedigree | `*.ped`, samplesheet `sex` column | Enables sex and relatedness checks |

Ask only if still ambiguous: assay (WGS/WES/panel/RNA/ChIP/ATAC/methylation/long-read), library kit when adapter autodetection fails, expected species, expected sex/pedigree, and whether downstream caller is DeepVariant/GATK (affects trimming decision).

## Decision table

| Situation | Choice | Why |
|---|---|---|
| Raw Illumina FASTQ, any assay | `fastp` (default) | Adapter autodetect for PE via overlap, poly-G, UMI, JSON for MultiQC, fastest |
| RRBS / bisulfite / methylseq, or need Bismark-compatible trimming | Trim Galore (`--rrbs`, `--clip_R2`) | Purpose-built presets; nf-core/methylseq uses it |
| Custom adapters, primers, linked adapters, demultiplexing | cutadapt directly | Explicit control; `-a`/`-A` per read, `-g` 5', `--pair-adapters` |
| Downstream DeepVariant, BWA-MEM2/minimap2 germline WGS/WES | Do NOT quality-trim; adapter-trim only if adapter content >5% | Aligners soft-clip adapters/low-quality tails; DeepVariant is trained on untrimmed reads; trimming shortens reads and changes insert-size distribution |
| 2-colour chemistry (NovaSeq/NextSeq/MiniSeq) | `fastp --trim_poly_g` (auto-enabled when instrument detected) | Dark cycles read as G; poly-G tails inflate mismatches |
| UMI in read (not in header) | `fastp --umi --umi_loc read1 --umi_len N` (or `per_read`, `index1`) | Moves UMI to read name `:UMI` for `fgbio`/`umi_tools` dedup downstream |
| Species/contamination doubt before alignment | FastQ Screen (subsample vs panel of genomes) then Kraken2 if unresolved | Screen answers "which genome"; Kraken2 answers "what microbial content" |
| Human WGS/WES contamination estimate | VerifyBamID2 (`FREEMIX`) | Ancestry-agnostic; standard in Broad/Illumina pipelines |
| Sample identity across assays (WGS vs RNA vs tumour/normal, re-runs) | somalier (BAM/CRAM/VCF) — NGSCheckMate when only FASTQ or when assays share few sites | somalier is fast and build-aware; NGSCheckMate works on FASTQ via k-mer SNP patterns |
| Pedigree + VCF already joint-called | peddy or `somalier relate --ped` | Checks reported vs inferred sex/relatedness; somalier recommended (peddy README says so) |
| Sex check without variants | mosdepth chrX/chrY ratio + `samtools idxstats` | Works on any BAM, no sites file |
| Coverage WGS | `mosdepth -n --fast-mode` | 5–10× faster than `samtools depth`; summary + distribution files |
| Coverage WES/panel | `mosdepth -n -b targets.bed --thresholds 1,10,20,30` + `gatk CollectHsMetrics` | Per-target % ≥20×, fold-80, on/off-target |
| Library complexity before deep sequencing | `preseq lc_extrap` on a shallow BAM | Predicts distinct reads at higher depth; justifies top-up vs re-prep |
| Long-read FASTQ/BAM | NanoPlot, pycoQC (ONT summary), cramino, `samtools stats` (N50) | Details in `.claude/skills/genomics-longread/SKILL.md` |
| RNA-seq | fastp + `samtools stats` here; RSeQC/strandedness in `.claude/skills/genomics-rnaseq/bulk.md` | Avoid duplication |
| Aggregation | MultiQC 1.35 with a YAML config | One report per project; feeds the PASS/WARN/FAIL table |

## Workflow

Set once per project: `THREADS=8`, `REF=GRCh38.fa`, `OUT=qc`, `SAMPLE=S1`. Every step appends to `${OUT}/versions.txt`.

1. Record versions and validate inputs.

```bash
set -euo pipefail
OUT=qc
mkdir -p "${OUT}"
bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh" fastqc fastp multiqc samtools mosdepth gatk VerifyBamID somalier > "${OUT}/versions.txt"
date -u +%Y-%m-%dT%H:%MZ >> "${OUT}/versions.txt"
for f in *.fastq.gz; do gzip -t "${f}"; done
for b in *.bam; do samtools quickcheck -v "${b}"; done
```

Validate: `cat ${OUT}/versions.txt` shows no `MISSING` for tools you will call.

2. Raw-read QC (FastQC 0.13 + fastp JSON; FastQC `-t` now = threads per file ×4 max per file).

```bash
set -euo pipefail
OUT=qc
THREADS=8
mkdir -p "${OUT}/fastqc"
fastqc -t "${THREADS}" --noextract -o "${OUT}/fastqc" *.fastq.gz
```

Validate: `ls ${OUT}/fastqc/*_fastqc.zip | wc -l` equals FASTQ count; `unzip -p ${OUT}/fastqc/S1_R1_fastqc.zip '*/summary.txt'` lists module PASS/WARN/FAIL.

3. Trim (only if the decision table says so). fastp 1.4.0 paired-end default:

```bash
set -euo pipefail
SAMPLE=S1
THREADS=8
OUT=qc
fastp -i "${SAMPLE}_R1.fastq.gz" -I "${SAMPLE}_R2.fastq.gz" \
  -o "${SAMPLE}_R1.trim.fastq.gz" -O "${SAMPLE}_R2.trim.fastq.gz" \
  --detect_adapter_for_pe --trim_poly_g --correction \
  --qualified_quality_phred 15 --length_required 36 \
  --thread "${THREADS}" --json "${OUT}/${SAMPLE}.fastp.json" --html "${OUT}/${SAMPLE}.fastp.html"
```

UMI variant: add `--umi --umi_loc per_read --umi_len 8 --umi_prefix UMI` (or `read1`/`index1` per kit). Adapter-only mode for DeepVariant/BWA inputs: drop `--qualified_quality_phred`, add `--disable_quality_filtering --disable_length_filtering`.

Trim Galore 2.6.0 (methylation / Bismark): `trim_galore --paired --cores 4 --rrbs -o trimmed/ R1.fastq.gz R2.fastq.gz` (`--fastqc_raw` adds pre-trim FastQC). cutadapt 5.2 explicit adapters: `cutadapt -j ${THREADS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 36 -o R1.trim.fastq.gz -p R2.trim.fastq.gz R1.fastq.gz R2.fastq.gz`.

Validate: `jq '.summary.after_filtering.q30_rate, .adapter_cutting.adapter_trimmed_reads' ${OUT}/${SAMPLE}.fastp.json`; re-run FastQC on trimmed output and confirm Adapter Content PASS.

4. Species / contamination screen on a read subsample.

```bash
set -euo pipefail
OUT=qc
THREADS=8
SAMPLE=S1
KRAKEN_DB=/dbs/k2_standard_08gb
fastq_screen --threads "${THREADS}" --subset 100000 --conf fastq_screen.conf --outdir "${OUT}/fqscreen" "${SAMPLE}_R1.fastq.gz"
seqtk sample -s42 "${SAMPLE}_R1.fastq.gz" 200000 > "${OUT}/${SAMPLE}.sub.fq"
kraken2 --db "${KRAKEN_DB}" --threads "${THREADS}" --report "${OUT}/${SAMPLE}.kraken2.report" --output /dev/null "${OUT}/${SAMPLE}.sub.fq"
```

`fastq_screen.conf` needs `DATABASE <name> <bowtie2 index prefix>` lines; `fastq_screen --get_genomes` downloads the default panel. Kraken2 DB: `<check current: https://benlangmead.github.io/aws-indexes/k2>`.

Validate: `head -5 ${OUT}/${SAMPLE}.kraken2.report` — first line is % unclassified; `grep -c . ${OUT}/fqscreen/*_screen.txt`.

5. Alignment metrics (post-alignment BAM; duplicates marked). Replaces `samtools depth | awk`.

```bash
set -euo pipefail
SAMPLE=S1
BAM="${SAMPLE}.bam"
REF=GRCh38.fa
OUT=qc
THREADS=8
samtools quickcheck "${BAM}"
samtools flagstat -@ "${THREADS}" "${BAM}" > "${OUT}/${SAMPLE}.flagstat"
samtools stats -@ "${THREADS}" -r "${REF}" "${BAM}" > "${OUT}/${SAMPLE}.stats"
samtools idxstats "${BAM}" > "${OUT}/${SAMPLE}.idxstats"
mosdepth -t "${THREADS}" -n --fast-mode "${OUT}/${SAMPLE}" "${BAM}"
gatk CollectAlignmentSummaryMetrics -R "${REF}" -I "${BAM}" -O "${OUT}/${SAMPLE}.alignment_metrics.txt"
gatk CollectInsertSizeMetrics -I "${BAM}" -O "${OUT}/${SAMPLE}.insert_size_metrics.txt" -H "${OUT}/${SAMPLE}.insert_size.pdf"
```

mosdepth forms (`-b`/`--by` are one option; give it once): windows `mosdepth -n --by 500 PFX BAM`; targets `mosdepth -n -b targets.bed --thresholds 1,10,20,30 PFX BAM` (writes `PFX.thresholds.bed.gz`); WGS add `gatk CollectWgsMetrics -R ${REF} -I ${BAM} -O wgs_metrics.txt`; WES add `gatk CollectHsMetrics -R ${REF} -I ${BAM} -O hs_metrics.txt --BAIT_INTERVALS bait.interval_list --TARGET_INTERVALS target.interval_list`; GC `gatk CollectGcBiasMetrics -R ${REF} -I ${BAM} -O gc_bias.txt -CHART gc_bias.pdf -S gc_summary.txt`.

Validate: `grep '^total' ${OUT}/${SAMPLE}.mosdepth.summary.txt` (mean depth); `grep ^SN ${OUT}/${SAMPLE}.stats | grep -E 'reads mapped:|error rate|insert size average'`; `gatk ValidateSamFile -I ${BAM} -M SUMMARY`.

6. Library complexity (shallow pilot BAM or any BAM before deciding to top-up).

```bash
set -euo pipefail
SAMPLE=S1
OUT=qc
preseq lc_extrap -B -P -e 2e9 -s 1e8 -o "${OUT}/${SAMPLE}.lc_extrap.txt" "${SAMPLE}.bam"
```

Validate: `awk 'NR==1 || $1==1e9' ${OUT}/${SAMPLE}.lc_extrap.txt` — expected distinct reads at 1 G total; falls far below the diagonal → saturated library.

7. Contamination (human, GRCh38 with `chr` prefix; VerifyBamID2 2.0.3). The SVD resource bundle ships in the package `resource/` directory as `1000g.phase3.100k.b38.vcf.gz.dat.{UD,mu,bed,V}`; use the `10k` panel for long reads, `exome/` for WES, `--Methylation` for WGBS/EM-seq.

```bash
set -euo pipefail
SAMPLE=S1
REF=GRCh38.fa
OUT=qc
SVD=/opt/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat
VerifyBamID --SVDPrefix "${SVD}" --Reference "${REF}" --BamFile "${SAMPLE}.bam" --NumPC 4 --NumThread 4 --Output "${OUT}/${SAMPLE}.vb2"
```

Validate: `cut -f7 ${OUT}/${SAMPLE}.vb2.selfSM | tail -1` = FREEMIX. The bioconda binary is `VerifyBamID` (the old `verifybamid2`/`--Ref` spelling is invalid). FREEMIX interpretation: <0.01 clean; 0.01–0.03 review (tumour, low depth <5×, or wrong ancestry panel inflate it); >0.03 contaminated for germline; somatic workflows commonly fail >0.05. Without `--Output` files default to `result.*` in CWD.

8. Identity, relatedness and sex (somalier 0.3.7). Sites: `https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz` (also `sites.hg38.nochr.vcf.gz` files/3412454, `sites.GRCh37.vcf.gz` files/3412455, `sites.chm13v2.T2T.vcf.gz` files/9954286, `sites.hg38.rna.vcf.gz` files/4566475 for RNA-seq).

```bash
set -euo pipefail
REF=GRCh38.fa
OUT=qc
SITES=sites.hg38.vcf.gz
mkdir -p "${OUT}/somalier"
for b in *.bam; do somalier extract -d "${OUT}/somalier" --sites "${SITES}" -f "${REF}" "${b}"; done
somalier relate --infer --sites "${SITES}" --ped cohort.ped -o "${OUT}/somalier/cohort" "${OUT}/somalier"/*.somalier
```

Validate: `awk -F'\t' 'NR>1 && $5!=$7' ${OUT}/somalier/cohort.samples.tsv` lists samples whose pedigree sex ≠ inferred sex (columns: 5=sex, 6=phenotype, 7=original_pedigree_sex; check column header names first: `head -1 ... | tr '\t' '\n' | nl`); `awk -F'\t' 'NR>1 && $3>0.8' ${OUT}/somalier/cohort.pairs.tsv` lists identical pairs — expected only for tumour/normal or replicates. RNA-seq: add `--min-ab 0.2`. Cross-assay: extract from every BAM/VCF with the same sites and relate together.

Alternatives: `peddy -p 4 --plot --sites hg38 --prefix ${OUT}/peddy joint.vcf.gz cohort.ped` (VCF+PED; writes `*.sex_check.csv`, `*.ped_check.csv`). NGSCheckMate 1.0.1 from FASTQ: `python ncm_fastq.py -l fastq_list.txt -pt SNP/SNP.pt -O ${OUT}/ncm -N cohort` (FASTQ-only or when assays lack shared called sites).

9. Sex check without a sites file (any BAM): X and Y coverage relative to autosomes.

```bash
set -euo pipefail
SAMPLE=S1
OUT=qc
awk -v s="${SAMPLE}" 'BEGIN{OFS="\t"} $1=="chrX"{x=$4} $1=="chrY"{y=$4} $1~/^chr([0-9]+)$/{a+=$4;n++} END{print s, x/(a/n), y/(a/n)}' "${OUT}/${SAMPLE}.mosdepth.summary.txt"
awk -v s="${SAMPLE}" 'BEGIN{OFS="\t"} $1=="chrY"{y=$3} {t+=$3} END{print s, y/t}' "${OUT}/${SAMPLE}.idxstats"
```

Decision (ratios normalised to diploid autosomal mean): X ≈1.0 and Y <0.05 → XX; X ≈0.5 and Y ≈0.5 → XY; X ≈0.5 and Y <0.05 → X0 or Y-masked reference; X ≈1.0 and Y ≈0.5 → XXY. WES: compute on-target only (`-b targets.bed`, use `*.regions.bed.gz` means). PAR regions and Y-pseudogene reads make Y fraction >0 in XX; use the ratio, not zero.

10. MultiQC 1.35 with a config that works.

```yaml
title: "Project QC"
subtitle: "Run 2026-11"
show_analysis_paths: false
show_analysis_time: false
fastqc_config:
  fastqc_theoretical_gc: "hg38_genome"
table_columns_visible:
  FastQC:
    percent_duplicates: true
    percent_gc: true
    avg_sequence_length: true
  "Samtools: stats":
    error_rate: true
    reads_mapped_percent: true
  mosdepth:
    median_coverage: true
custom_data:
  sex_check:
    file_format: "tsv"
    section_name: "Sex check"
    description: "chrX and chrY depth relative to autosomal mean (mosdepth summary)."
    plot_type: "table"
    pconfig:
      id: "sex_check_table"
      namespace: "Sex check"
    headers:
      x_ratio:
        title: "X ratio"
        min: 0
        max: 1.5
      y_ratio:
        title: "Y ratio"
        min: 0
        max: 1
sp:
  sex_check:
    fn: "*sex_check_mqc.tsv"
```

The custom TSV needs a header line `sample\tx_ratio\ty_ratio` and the filename must match `sp.sex_check.fn`. Run:

```bash
set -euo pipefail
OUT=qc
multiqc "${OUT}" -c multiqc_config.yaml -o "${OUT}/multiqc" -f -m fastqc -m fastp -m samtools -m mosdepth -m picard -m verifybamid -m somalier -m fastq_screen -m kraken -m preseq -m custom_content
```

Module keys are lowercase as printed by `multiqc --list-modules` (`featurecounts`, `star`, `rseqc`, `cutadapt`, `trimgalore` for RNA/other inputs; `featureCounts` is not a valid key). Validate: `ls ${OUT}/multiqc/multiqc_data/multiqc_general_stats.txt` and `grep -c . ${OUT}/multiqc/multiqc_data/multiqc_sources.txt` equals expected file count.

11. Batch structure / outliers (on `multiqc_general_stats.txt` plus a sample→batch map).

```r
library(dplyr)
library(ggplot2)
gs <- read.delim("qc/multiqc/multiqc_data/multiqc_general_stats.txt", check.names = FALSE)
meta <- read.delim("samples.tsv")            # columns: Sample, batch
df <- inner_join(gs, meta, by = "Sample")
num <- df %>% select(where(is.numeric)) %>% select(where(~ sd(.x, na.rm = TRUE) > 0))
num <- num[complete.cases(num), , drop = FALSE]
keep <- complete.cases(df %>% select(all_of(names(num))))
pca <- prcomp(num, scale. = TRUE)
scores <- as.data.frame(pca$x[, 1:2]) %>% mutate(Sample = df$Sample[keep], batch = df$batch[keep])
ggplot(scores, aes(PC1, PC2, colour = batch, label = Sample)) + geom_point(size = 3) + geom_text(vjust = -0.8, size = 2.5) + theme_minimal()
z <- scale(num)
outliers <- data.frame(Sample = df$Sample[keep], n_metrics_gt3sd = rowSums(abs(z) > 3)) %>% filter(n_metrics_gt3sd > 0) %>% arrange(desc(n_metrics_gt3sd))
print(outliers)
```

Validate: `summary(pca)$importance[2, 1:2]` — if PC1 is dominated by one metric (check `pca$rotation[, 1]`), report that metric rather than "batch effect". Batch separation on PC1/PC2 → include batch as a covariate downstream; do not "correct" raw QC values.

## QC thresholds

PASS/WARN/FAIL defaults by assay; the rationale column says when to deviate.

| Metric | Default (PASS / WARN / FAIL) | Rationale / when to deviate |
|---|---|---|
| Q30 bases (fastp `q30_rate`) | ≥0.85 / 0.75–0.85 / <0.75 | Illumina spec ≥80% Q30; long inserts and 2×250 runs legitimately lower; ONT/PacBio use own metrics |
| Adapter-trimmed reads (fastp) | <5% / 5–20% / >20% | >20% implies short inserts (degraded DNA, cfDNA, small RNA, ATAC — expected there; raise FAIL) |
| Mean read length after trim | ≥0.9× nominal / 0.7–0.9× / <0.7× | Aggressive trimming lowers mappability; for DeepVariant inputs expect ≈1.0× |
| Mapping rate (`samtools flagstat` mapped %) | WGS ≥98 / 95–98 / <95; WES ≥98; RNA ≥85 / 75–85 / <75; ChIP/ATAC ≥90 | Low rate + high Kraken2 non-host → contamination; low rate + correct species → wrong reference build or adapters |
| Properly paired | ≥95% / 90–95% / <90% | Lower with SV-rich tumours, long inserts; check insert-size SD first |
| Duplicate rate (marked, not removed) | WGS PCR-free <10 / 10–20 / >20%; WES <25 / 25–40 / >40%; RNA n/a (high is expected); ChIP <30%; ATAC <40% | Duplication rises with depth; compare to preseq expectation rather than a fixed cut; UMI libraries report post-consensus |
| Mean depth (mosdepth `total`) | WGS germline ≥30× / 20–30 / <20; WGS somatic tumour ≥60×; WES ≥80× mean, ≥20× at ≥90% targets; panel ≥500× | Caller requirement, not sequencing spec; DeepVariant OK at 20×; low-pass WGS is a different design |
| % targets ≥20× (`thresholds.bed.gz` or HsMetrics `PCT_TARGET_BASES_20X`) | ≥95 / 90–95 / <90 | Clinical exomes often require ≥98% at 20×; GC-extreme targets drag this down — report fold-80 |
| Fold-80 base penalty (HsMetrics) | ≤2 / 2–3 / >3 | Capture uniformity; >3 → rebalance probes or increase depth |
| Insert size (median, PE) | kit spec ±20% / ±40% / bimodal or <2× read length | <2× read length means adapters read-through: must trim |
| Error rate (`samtools stats`) | <0.5% / 0.5–1% / >1% | Mismatch-dominated error → wrong reference/species or contamination; indel → homopolymer chemistry |
| FREEMIX (VerifyBamID2) | <0.01 / 0.01–0.03 / >0.03 | Germline; somatic normal often enforced <0.02; low depth (<5×) and wrong panel inflate |
| Sex concordance | match / inferred ambiguous / mismatch | Mismatch is a sample swap until proven otherwise; Klinefelter/Turner are rare but real — confirm with somalier het-rate on X |
| Relatedness (somalier, unexpected pairs) | <0.1 / 0.1–0.4 / >0.4 unexpected | >0.9 = duplicate or swap; 0.4–0.6 unexpected = relatives or cross-contamination |
| Kraken2 non-host fraction | <1% / 1–5% / >5% | Saliva/stool/FFPE/microbiome samples legitimately higher; look at genus-level hits |
| preseq distinct at 2× current depth | ≥1.6× current / 1.3–1.6× / <1.3× | <1.3× means top-up sequencing mostly yields duplicates; re-prep instead |
| FastQC per-base quality | all bases ≥Q28 / some <Q20 at 3' / many <Q20 | 3' decay is normal; a mid-read dip suggests bubbles/flowcell issue |
| FastQC overrepresented sequences | <1% / 1–10% / >10% | rRNA in RNA-seq, primers in amplicon are expected — identify before failing |
| rRNA / mito fraction (RNA-seq) | see `.claude/skills/genomics-rnaseq/bulk.md` | Assay-specific |
| Long-read N50 / Q-score | see `.claude/skills/genomics-longread/SKILL.md` | Platform-specific |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `mosdepth: error parsing arguments` | `-b` and `--by` given together (they are one flag) | Use once: `-b targets.bed` or `--by 500` |
| `multiqc` ignores a module | Key spelled like the tool name (`featureCounts`, `Picard`) | Use lowercase keys from `multiqc --list-modules`; check `multiqc_data/multiqc.log` |
| `VerifyBamID: unrecognized option '--Ref'` | Old flag names | `--Reference`, `--BamFile`, `--SVDPrefix`; binary is `VerifyBamID`, not `verifybamid2` |
| FREEMIX 0.02–0.05 on a clean sample | Wrong build panel (b37 vs b38), <5× depth, tumour LOH, long reads on 100k panel | Match panel to build/assay; use 10k panel for long reads; corroborate with `somalier relate` self-pairs and CHARR |
| somalier sex column all `unknown` | `--sites` omitted at `relate`, or BAM lacks chrX reads (targeted panel) | Pass `--sites` to both steps; for panels use idxstats/mosdepth ratio |
| somalier relatedness 0 for everything | Extracted with mismatched sites/build (chr vs no-chr) | Use `sites.hg38.vcf.gz` vs `sites.hg38.nochr.vcf.gz` to match the BAM; `ref_sanity.sh` first |
| DeepVariant/GATK calls drop after "cleanup" | Quality trimming shortened reads, altered insert size, removed soft-clip evidence | Re-run with adapter-only or no trimming; aligners soft-clip |
| Poly-G tails only in R2, low mapping of R2 | 2-colour chemistry dark cycles | `fastp --trim_poly_g` (auto on NovaSeq/NextSeq); do not blame the library |
| Duplicate rate 0% in a WGS BAM | Duplicates never marked (no `MarkDuplicates` in `@PG`) | Run `gatk MarkDuplicates` or `samtools markdup` before QC |
| Duplicate rate 60%+ in RNA-seq | Expected for highly expressed genes | Judge with `dupRadar`/gene-body metrics, not raw %; see rnaseq skill |
| `samtools stats` error rate 2% on human | Mismatched reference (hg19 BAM vs GRCh38 FASTA via `-r`) | `ref_sanity.sh ref.fa sample.bam`; rerun with the aligning reference |
| Sex ratio X≈0.5 but somalier says female | PAR masking / `chrX` aliases (`X` vs `chrX`) in awk | Match contig names from `idxstats`; exclude PAR via `-b` with non-PAR BED |
| FastQ Screen shows 100% "no hits" | Config paths wrong / indexes missing | `fastq_screen --get_genomes`; check `DATABASE` lines point to Bowtie2 prefix |
| MultiQC custom table missing | `sp.<id>.fn` pattern does not match or TSV lacks header | Name file `*_mqc.tsv`, header row first, same id under `custom_data` and `sp` |
| `gatk CollectHsMetrics` fails on BED | Needs `.interval_list` | `gatk BedToIntervalList -I targets.bed -SD ref.dict -O targets.interval_list` |
| `fastqc` `-t 16` slower than expected (0.13) | Threads per file capped at 4; extra threads parallelise files | Pass many files per call |

## nf-core

There is no `nf-core/fastqc`. QC comes from the assay pipeline plus MultiQC:

| Pipeline | QC it produces | Samplesheet header |
|---|---|---|
| `nf-core/sarek` | FastQC, fastp (`--trim_fastq`), mosdepth, samtools stats, MultiQC; `--skip_tools baserecalibrator,...` to shorten | `patient,sample,lane,fastq_1,fastq_2` |
| `nf-core/rnaseq` | FastQC, Trim Galore/fastp (`--trimmer fastp`), RSeQC, Qualimap, dupRadar, preseq, MultiQC; `--skip_alignment` for read-level QC only | `sample,fastq_1,fastq_2,strandedness` |
| `nf-core/methylseq` | Trim Galore (`--rrbs`), Bismark reports, MultiQC | `sample,fastq_1,fastq_2` |
| `nf-core/nanoseq` | NanoPlot, pycoQC, MultiQC for long reads | `group,replicate,barcode,input_file,fasta,gtf` |

Pipeline mechanics (profiles, `-resume`, configs, SLURM): `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

Per-step idioms — run after every producing step, fail the script on non-zero:

```bash
set -euo pipefail
BAM=S1.bam
VCF=S1.vcf.gz
samtools quickcheck -v "${BAM}"
gatk ValidateSamFile -I "${BAM}" -M SUMMARY
bcftools view -H "${VCF}" | head -2
bcftools stats -s - "${VCF}" | grep -E '^SN|^TSTV' | head -20
```

- Truth: on GIAB/HG001–HG007 or Platinum Genomes BAMs, expect FREEMIX <0.005, sex matching NIST metadata, somalier self-relatedness ≈1 across replicates. Variant-level truth-set benchmarking: `.claude/skills/genomics-benchmark/SKILL.md`.
- Sanity: mosdepth `total` mean ≈ (read pairs × 2 × read length × mapped fraction × (1 − dup rate)) / 3.1e9 within ±10%; `samtools stats` `insert size average` ≈ kit spec; Ti/Tv 2.0–2.1 WGS, 2.8–3.3 WES.
- Report every sample as PASS/WARN/FAIL per metric in the table above, with the threshold used, and keep `versions.txt` next to the MultiQC report.

QC report template (`qc/QC_REPORT.md`):

```text
# QC report — <project> — <date>
Tools: see qc/versions.txt. Reference: <REF, build, source URL>.
## Summary: N samples; PASS a; WARN b; FAIL c
## Per-sample table
sample | reads (M) | Q30 | adapter % | mapped % | dup % | mean depth | %≥20x | insert median | FREEMIX | sex (expected/inferred) | status | failing metric(s)
## Failed / flagged samples (one line each: metric, value, threshold, probable cause, action)
## Identity & relatedness: unexpected pairs from somalier pairs.tsv
## Contamination / species: Kraken2 non-host %, FastQ Screen top genome
## Batch structure: PCA PC1/PC2 variance, driving metric, batch separation yes/no
## Actions: re-sequence / top-up / re-prep / investigate swap / proceed
```

## Additional resources

- `.claude/skills/genomics/SKILL.md` — router and eight rules; `scripts/peek.sh`, `ref_sanity.sh`, `check_env.sh`. Load first.
- `.claude/skills/genomics/reference.md` — reference/annotation downloads, GIAB/Platinum truth BAM acquisition.
- `.claude/skills/genomics-rnaseq/bulk.md` — RSeQC, strandedness inference, rRNA/dupRadar QC for RNA-seq.
- `.claude/skills/genomics-longread/SKILL.md` — NanoPlot/pycoQC/cramino, N50 and Q-score thresholds, long-read contamination panels.
- `.claude/skills/genomics-benchmark/SKILL.md` — hap.py/truvari truth-set benchmarking after QC passes.
- `.claude/skills/genomics-workflow/SKILL.md` — nf-core/Nextflow mechanics for running the pipelines named above.
- `.claude/skills/genomics-somatic/SKILL.md` — tumour/normal contamination (`GetPileupSummaries`/`CalculateContamination`) and T/N identity expectations.
