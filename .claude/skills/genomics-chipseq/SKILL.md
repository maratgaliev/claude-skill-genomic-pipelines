---
name: genomics-chipseq
description: ChIP-seq peak calling, QC, differential binding and motif analysis (MACS3, IDR, deepTools, DiffBind, ChIPseeker, HOMER/MEME), with routing to ATAC-seq, CUT&RUN/CUT&Tag and scATAC supporting files. Use when inputs are IP/input FASTQ or BAM pairs, narrowPeak/broadPeak/bedGraph files, or the user mentions antibody, histone mark, Tn5, accessibility, SEACR, IgG, spike-in, footprinting. Not for raw-read QC/trimming (genomics-qc), methylation (genomics-methylation), or RNA-seq (genomics-rnaseq).
argument-hint: "[chip|atac|cutandrun|scatac] [path]"
paths:
  - "**/*.narrowPeak"
  - "**/*.broadPeak"
  - "**/*.gappedPeak"
  - "**/*_peaks.xls"
  - "**/*.bedgraph"
  - "**/*.bdg"
  - "**/*fragments.tsv.gz"
---

# ChIP-seq, ATAC-seq and CUT&RUN analysis

Expert in chromatin profiling: ChIP-seq (TF and histone), ATAC-seq, CUT&RUN/CUT&Tag and single-cell ATAC. Outputs are pinned command blocks, decision tables, ENCODE-derived QC gates and validation checks. This file holds the ChIP-seq core; assay-specific deviations live in `atac.md` and `cutandrun.md` in this directory.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" bowtie2 bwa-mem2 samtools gatk bedtools macs3 idr deeptools bamCoverage plotFingerprint run_spp.R findMotifsGenome.pl meme-chip Rscript nextflow || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| Input | Read | Infer |
|---|---|---|
| FASTQ | `peek.sh R1.fq.gz`; read-name suffixes, read length, presence of `_R2`/`_2` | layout (SE/PE), platform, whether UMIs/barcodes present |
| BAM/CRAM | `samtools view -H` (`@SQ`, `@PG`, `@RG`) | assembly (chr1 vs 1, GRCh38 vs T2T, mm10 vs GRCm39), aligner + flags already applied, dedup already run (`@PG ID:MarkDuplicates`/`samtools markdup`), sorted-by |
| BAM flags | `samtools flagstat`; `samtools view -c -f 1024` | PE vs SE, duplicate marking present, proper-pair rate |
| Peaks | `peek.sh x.narrowPeak`; column count (10 = narrowPeak, 9 = broadPeak) | caller, whether summit column is populated, `-log10(q)` scale |
| Design | sample names, antibody column, control mapping | TF vs histone (narrow vs broad), replicate count, which input/IgG pairs with which IP |

Ask only if still ambiguous: antibody target class (TF / narrow histone / broad histone), which sample is the matched input/IgG, whether a spike-in genome was added, and whether the design is ChIP-seq vs CUT&RUN vs ATAC (fragment-size histogram and no-input design usually answer this).

## Decision table

| Situation | Choice |
|---|---|
| TF or narrow histone (H3K4me3, H3K27ac, H3K9ac, H2A.Z) | `macs3 callpeak` narrow, `-q 0.05`; replicate consistency via IDR |
| Broad histone (H3K27me3, H3K36me3, H3K9me3, H3K4me1, H3K79me2) | `macs3 callpeak --broad --broad-cutoff 0.1`; no IDR, use naive overlap / pseudo-replicates |
| Paired-end | `-f BAMPE` (fragment pileup, `--nomodel` implied); never add `--shift` |
| Single-end | `-f BAM`, let MACS3 build the model; if `predictd` fails use `--nomodel --extsize <phantompeak fragment length>` |
| No replicates | pseudo-replicates (split reads 50/50) + IDR on pseudo-reps; label results as unreplicated |
| Tn5-based accessibility | `atac.md` |
| CUT&RUN / CUT&Tag / IgG control / E. coli spike-in | `cutandrun.md` |
| scATAC / 10x multiome | `atac.md` § scATAC |
| Global occupancy change expected (drug, degron, KO of writer) | spike-in normalised (ChIP-Rx) or `DBA_NORM_LIB` on background bins; never RLE/TMM on consensus peaks |
| Differential binding, ≥2 reps per group | DiffBind v3 (DESeq2 default) |
| Differential binding, 1 rep per group | no statistics; report fold-change ranking only and say so |
| Aligner | Bowtie2 `--very-sensitive` default; BWA-MEM2 acceptable; nf-core/chipseq offers both plus chromap |

## Workflow

Alignment of reads and raw-read QC/trimming are not repeated here; trim per `.claude/skills/genomics-qc/SKILL.md`, then:

### 1. Align and filter

```bash
set -euo pipefail
REF_IDX=/refs/GRCh38/bowtie2/GRCh38
SAMPLE=IP_rep1
THREADS=8
bowtie2 -p "${THREADS}" --very-sensitive -x "${REF_IDX}" \
  -1 "${SAMPLE}_R1.trim.fq.gz" -2 "${SAMPLE}_R2.trim.fq.gz" \
  --no-mixed --no-discordant -X 1000 2> "${SAMPLE}.bowtie2.log" \
  | samtools view -b -q 30 -F 1804 -f 2 - \
  | samtools sort -@ 4 -o "${SAMPLE}.sorted.bam" - && samtools quickcheck "${SAMPLE}.sorted.bam"
samtools index "${SAMPLE}.sorted.bam"
```

Single-end: `-U reads.fq.gz`, filter `-q 30 -F 1796` (drop `-f 2`). `-F 1804` = unmapped, mate unmapped, secondary, QC-fail, duplicate. Validation: `grep overall "${SAMPLE}.bowtie2.log"` (expect >80 %).

### 2. Mark duplicates (choose one), then exclude by flag

Duplicates are **marked** (flag 1024) and excluded downstream with `-F 1804` / `--ignoreDuplicates` (MACS3 ignores them by default); the marked BAM keeps the reads for complexity metrics and for assays where "duplicates" are real (CUT&RUN/Tag, see `cutandrun.md`). Physically removing them (`--REMOVE_DUPLICATES true`, `markdup -r`) is a data-discarding step: do it only on the user's instruction.

```bash
set -euo pipefail
SAMPLE=IP_rep1
gatk MarkDuplicates -I "${SAMPLE}.sorted.bam" -O "${SAMPLE}.dedup.bam" \
  -M "${SAMPLE}.dup_metrics.txt" --REMOVE_DUPLICATES false --CREATE_INDEX true
```

```bash
set -euo pipefail
SAMPLE=IP_rep1
samtools collate -@ 4 -O -u "${SAMPLE}.sorted.bam" \
  | samtools fixmate -@ 4 -m -u - - \
  | samtools sort -@ 4 -u - \
  | samtools markdup -@ 4 -s - "${SAMPLE}.dedup.bam"
samtools index "${SAMPLE}.dedup.bam"
```

`samtools markdup` requires `fixmate -m` (MC/ms tags) on name-collated input; skipping it silently under-marks. Validation: `grep -A1 LIBRARY "${SAMPLE}.dup_metrics.txt"` or `samtools flagstat` duplicate line; record PERCENT_DUPLICATION (ENCODE NRF/PBC1/PBC2 targets in QC table).

### 3. Duplicate-flag and blacklist filtering, re-index

`filt.bam` is the analysis BAM: marked duplicates dropped by flag (the marked `dedup.bam` stays on disk), blacklist regions removed.

```bash
set -euo pipefail
SAMPLE=IP_rep1
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
samtools view -b -F 1024 "${SAMPLE}.dedup.bam" \
  | bedtools intersect -v -abam stdin -b "${BLACKLIST}" > "${SAMPLE}.filt.bam"
samtools index "${SAMPLE}.filt.bam"
samtools quickcheck "${SAMPLE}.filt.bam"
```

ENCODE blacklist v2 (Amemiya 2019): `https://github.com/Boyle-Lab/Blacklist/tree/master/lists` — `hg38-blacklist.v2.bed.gz`, `hg19-blacklist.v2.bed.gz`, `mm10-blacklist.v2.bed.gz`, `dm6-blacklist.v2.bed.gz`, `ce11-blacklist.v2.bed.gz`. Verified 2026-11-22: no GRCm39 list exists in that repo; for GRCm39 either liftOver mm10 v2 with `CrossMap`/UCSC `liftOver` (record chain file) or use Bioconductor `excluderanges` (`<check current: https://bioconductor.org/packages/excluderanges>`), and state which in the methods. hg38 v3 (`ENCFF356LFX`) is available from the ENCODE portal; v2 remains the default in nf-core. The `bedtools intersect -abam` output loses the index: `plotFingerprint`, `bamCoverage` and MACS3 `BAMPE` fail or crawl without `samtools index`. Validation: `samtools idxstats "${SAMPLE}.filt.bam" | head -3`.

### 4. Enrichment QC

```bash
set -euo pipefail
SAMPLE=IP_rep1
CONTROL=INPUT_rep1
plotFingerprint -b "${SAMPLE}.filt.bam" "${CONTROL}.filt.bam" --labels IP Input \
  -p 8 --outQualityMetrics "${SAMPLE}.fingerprint.tsv" -o "${SAMPLE}.fingerprint.png"
run_spp.R -c="${SAMPLE}.filt.bam" -savp="${SAMPLE}.spp.pdf" -out="${SAMPLE}.spp.txt" -p=8
```

phantompeakqualtools (`run_spp.R`) wants SE-style BAM; for PE, pass R1 only (`samtools view -b -f 64`). Columns 9–11 of `spp.txt` are NSC, RSC, quality tag. Validation: `cut -f9-11 "${SAMPLE}.spp.txt"`.

### 5. Peak calling (MACS3 3.1.0)

```bash
set -euo pipefail
SAMPLE=IP_rep1
CONTROL=INPUT_rep1
GSIZE=hs
macs3 callpeak -t "${SAMPLE}.filt.bam" -c "${CONTROL}.filt.bam" -f BAMPE -g "${GSIZE}" \
  -n "${SAMPLE}" --outdir peaks -q 0.05 -B --SPMR --keep-dup all
macs3 callpeak -t "${SAMPLE}.filt.bam" -c "${CONTROL}.filt.bam" -f BAMPE -g "${GSIZE}" \
  -n "${SAMPLE}.broad" --outdir peaks --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all
```

Run only the line matching the target class. `-g`: `hs` 2.7e9, `mm` 1.87e9, `dm` 1.2e8, `ce` 9e7; for other genomes supply the effective size (e.g. from `khmer unique-kmers.py`). `--keep-dup all` because duplicates were already removed; otherwise MACS3 re-deduplicates to 1 per position. For IDR input, call relaxed peaks with `-p 0.01` instead of `-q 0.05`. Validation: `wc -l peaks/${SAMPLE}_peaks.narrowPeak` and FRiP:

```bash
set -euo pipefail
SAMPLE=IP_rep1
TOTAL=$(samtools view -c -F 260 "${SAMPLE}.filt.bam")
INPEAK=$(bedtools sort -i "peaks/${SAMPLE}_peaks.narrowPeak" | bedtools merge -i - | bedtools intersect -u -a "${SAMPLE}.filt.bam" -b - | samtools view -c -F 260 -)
awk -v a="${INPEAK}" -v b="${TOTAL}" 'BEGIN{printf "FRiP\t%.4f\n", a/b}'
```

### 6. Replicate consistency

```bash
set -euo pipefail
mkdir -p idr
idr --samples peaks/IP_rep1_peaks.narrowPeak peaks/IP_rep2_peaks.narrowPeak \
  --input-file-type narrowPeak --rank p.value --output-file idr/rep1_vs_rep2.idr.txt \
  --plot --log-output-file idr/rep1_vs_rep2.idr.log
awk '$5>=540' idr/rep1_vs_rep2.idr.txt > idr/rep1_vs_rep2.idr0.05.narrowPeak
```

Column 5 = `min(int(-125*log2(IDR)),1000)`; `≥540` ⇔ IDR ≤ 0.05 (ENCODE threshold). Rank by `p.value` for MACS3 relaxed peaks (q-values tie). IDR is for narrow peaks only; broad marks use naive overlap (`bedtools intersect -f 0.5 -r`). Pseudo-replicates when unreplicated: `samtools view -b -s 42.5` and `-s 43.5` are **not** disjoint; instead shuffle a name-sorted BED of reads and `split -n l/2`, call peaks on each half, run IDR on the halves. Validation: ENCODE self-consistency and rescue ratios both < 2.

### 7. Signal tracks and profiles (deepTools)

```bash
set -euo pipefail
SAMPLE=IP_rep1
CONTROL=INPUT_rep1
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
bamCoverage -b "${SAMPLE}.filt.bam" -o "${SAMPLE}.cpm.bw" --normalizeUsing CPM --binSize 10 \
  --blackListFileName "${BLACKLIST}" --extendReads -p 8
bamCompare -b1 "${SAMPLE}.filt.bam" -b2 "${CONTROL}.filt.bam" -o "${SAMPLE}.log2ratio.bw" \
  --operation log2 --scaleFactorsMethod readCount --binSize 10 --extendReads -p 8
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 -S "${SAMPLE}.cpm.bw" \
  -R /refs/GRCh38/genes.bed --skipZeros -p 8 -o "${SAMPLE}.tss.mat.gz"
plotHeatmap -m "${SAMPLE}.tss.mat.gz" -o "${SAMPLE}.tss.heatmap.png"
plotProfile -m "${SAMPLE}.tss.mat.gz" -o "${SAMPLE}.tss.profile.png"
```

`--normalizeUsing RPGC --effectiveGenomeSize 2913022398` (GRCh38, 100 bp reads) when 1× coverage is wanted. For SE data give `--extendReads <fragment length>` from step 4. Validation: `bigWigInfo "${SAMPLE}.cpm.bw" | head`.

### 8. Spike-in normalisation (ChIP-Rx, Orlando 2014)

Align to a combined host+spike index or separately to dm6; count spike reads with `samtools view -c -q 30 -F 1804 in.bam $(cut -f1 dm6.fai)`. Scale factor per sample = `1e6 / spike_reads`; apply with `bamCoverage --scaleFactor` and, in DiffBind, `dba.normalize(..., spikein = <dm6 BAMs>)` or `library = <spike counts>`. Requires equal cell numbers per IP and the same antibody lot; otherwise the correction is meaningless. CUT&RUN uses E. coli — see `cutandrun.md`.

### 9. Differential binding (DiffBind ≥ 3.0)

```r
library(DiffBind)
samples <- read.csv("diffbind_samplesheet.csv")
# required columns: SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
dba <- dba(sampleSheet = samples, minOverlap = 2)
dba <- dba.blacklist(dba, blacklist = DBA_BLACKLIST_HG38, greylist = TRUE)
dba <- dba.count(dba, summits = 250, bParallel = TRUE)
dba <- dba.normalize(dba, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL, background = TRUE)
dba <- dba.contrast(dba, categories = DBA_CONDITION, minMembers = 2)
dba <- dba.analyze(dba, method = DBA_DESEQ2)
res <- dba.report(dba, th = 0.05, bCalled = TRUE)
dba.plotMA(dba); dba.plotVolcano(dba); dba.plotPCA(dba, attributes = DBA_CONDITION)
```

Normalisation choice: `DBA_NORM_LIB` + `background = TRUE` (library size over genome-wide bins) is the safe default and the only one that preserves a global loss/gain; `DBA_NORM_RLE`/`DBA_NORM_TMM` on consensus-peak reads assume most peaks are unchanged and will erase a genome-wide shift (DiffBind vignette §7). `minMembers = 2` is the floor for a contrast; with one replicate per condition DiffBind refuses (`minMembers = 1` disables dispersion estimation — do not use it, report fold changes only). `summits = 250` re-centres on ±250 bp for TF data; drop it for broad marks. `bCalled = TRUE` adds per-condition peak-call counts; filter on them before interpretation. Validation: `dba.show(dba, bContrasts = TRUE)` shows DB counts; `dba.plotHeatmap(dba, correlations = TRUE)` clusters by condition.

### 10. Annotation and motifs

```r
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(clusterProfiler)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile("idr/rep1_vs_rep2.idr0.05.narrowPeak")
anno <- annotatePeak(peaks, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoPie(anno); plotDistToTSS(anno)
ego <- enrichGO(unique(as.data.frame(anno)$geneId), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
```

For GRCm39 use `TxDb.Mmusculus.UCSC.mm39.knownGene` (mm10 = legacy alias; do not mix mm10 peaks with mm39 TxDb; `ref_sanity.sh` catches length mismatches).

```bash
set -euo pipefail
PEAKS=idr/rep1_vs_rep2.idr0.05.narrowPeak
REF=/refs/GRCh38/GRCh38.fa
JASPAR=/refs/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
awk 'BEGIN{OFS="\t"}{s=$2+$10; print $1,s-100,s+100,$4}' "${PEAKS}" > summits200.bed
findMotifsGenome.pl summits200.bed "${REF}" homer_out -size given -p 8 -mask
bedtools getfasta -fi "${REF}" -bed summits200.bed -fo summits200.fa
meme-chip -oc memechip_out -meme-nmotifs 5 -db "${JASPAR}" summits200.fa
```

JASPAR 2024 is current (`https://jaspar.elixir.no/downloads/`; MEME-format CORE collections under `download/data/2024/CORE/`). HOMER accepts a FASTA path instead of a preconfigured genome name, which avoids `configureHomer.pl` drift. Validation: `head -20 homer_out/knownResults.txt`; the expected TF family should rank in the top 5 for a TF ChIP.

## QC thresholds

Values are ENCODE (TF/histone standards, 2012 quality-metrics page, data-standards/terms) unless noted. "Usable fragment" = uniquely mapped, non-duplicate, non-blacklist.

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Usable fragments, TF / narrow histone | ≥ 20 M recommended; 10–20 M low; < 5 M extremely low | ENCODE4; per replicate |
| Usable fragments, broad histone | ≥ 45 M recommended; ≥ 35 M acceptable; floor 20 M | ENCODE4; H3K9me3: 45 M total mapped (repeat-rich) |
| Mapping rate | > 80 % | lower suggests contamination/adapter; check `genomics-qc` |
| NRF / PBC1 / PBC2 | > 0.9 / > 0.9 / > 10 ideal; 0.8–0.9 / 0.8–0.9 / 3–10 compliant; < 0.5 / < 0.5 / < 1 severe | ENCODE library-complexity table; low-input ChIP legitimately lower |
| NSC | ≥ 1.05 (≥ 1.1 strong) | ENCODE 2012 metrics; broad marks and low-peak-count factors sit near 1.05 by nature |
| RSC | ≥ 0.8 (≥ 1.0 strong) | `< 0.8` = low signal-to-noise or phantom peak dominates; quality tag −2…2 from `run_spp.R` |
| FRiP, TF (point-source) | ≥ 1 % minimum (Landt 2012); typical good TF 5–30 % | ENCODE4 sets no hard threshold; compare only within identical caller settings. Report the exact value, not a pass/fail word |
| FRiP, narrow histone (H3K4me3, H3K27ac) | ≥ 10 % typical | not an ENCODE standard; derived from practice, state as such |
| FRiP, broad histone | not thresholded | peak-set dependent; use fingerprint elbow + NSC instead |
| IDR self-consistency / rescue ratio | both < 2 ideal; one < 2 acceptable; both > 2 concerning | ENCODE replicate concordance |
| IDR peaks (TF) | report count at IDR ≤ 0.05 | no fixed minimum; < 1000 for a sequence-specific TF warrants a re-look at the antibody |
| Fingerprint | IP curve departs from input; input ≈ diagonal | input that looks enriched implies copy-number or sonication bias |
| Duplicate rate | < 20 % typical | high values with high NRF just mean over-sequencing; high with low NRF = bottleneck |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `macs3 callpeak` finds almost no peaks on PE data with `-f BAM` | reads treated as SE tags; model built on half-fragments | `-f BAMPE` (or `BEDPE`) for paired data |
| Thousands of peaks on chrM / satellite / chrUn | no blacklist filter | step 3; add `--blackListFileName` to deepTools; `dba.blacklist()` in DiffBind |
| `plotFingerprint`: "index not found" after blacklist filtering | `bedtools intersect -abam` writes an unindexed BAM | `samtools index` (step 3) |
| `samtools markdup` marks ~0 % duplicates | missing `fixmate -m` or input not name-collated | step 2 variant B exactly |
| `gatk MarkDuplicates INPUT=…` errors | Picard-style `KEY=value` syntax | `-I/-O/-M` syntax |
| MACS3 peaks halved vs MACS2 run | MACS3 default `--keep-dup 1` on an already deduped BAM is harmless; the difference is usually `-f BAM` vs `BAMPE` or `--nomodel` defaults | compare `_peaks.xls` headers for the parameters actually used |
| IDR aborts with "too few peaks" | input peak sets already thresholded at q 0.05 | call relaxed peaks (`-p 0.01`) for IDR input |
| DiffBind finds no DB sites despite obvious global loss | `DBA_NORM_RLE`/`TMM` on consensus peaks normalised the loss away | `DBA_NORM_LIB` + `background = TRUE`, or spike-in |
| `dba.contrast` error "no contrasts" | fewer than `minMembers` (2) replicates in a group | add replicates; otherwise no statistics |
| All peaks annotate to "Distal Intergenic" | TxDb assembly mismatch (mm10 peaks, mm39 TxDb; or `chr` prefix) | `ref_sanity.sh`; `seqlevelsStyle(peaks) <- "UCSC"` |
| HOMER "genome not found" | preconfigured genome absent | pass the FASTA path as genome |
| `run_spp.R` RSC < 0.5 on PE data | mates treated as independent tags | run on R1 only |
| Spike-in scaling makes replicate differences larger | unequal cell numbers or spike chromatin per IP | drop spike-in; fall back to `DBA_NORM_LIB` background |

## nf-core

Pipelines verified in `https://nf-co.re/pipelines.json`: `nf-core/chipseq` (2.1.0), `nf-core/atacseq` (2.1.2), `nf-core/cutandrun` (3.2.2). Nextflow/profile/`-resume`/cluster mechanics: `.claude/skills/genomics-workflow/SKILL.md`.

```csv
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_H3K27AC,WT_H3K27AC_R1.fastq.gz,WT_H3K27AC_R2.fastq.gz,1,H3K27AC,WT_INPUT,1
WT_INPUT,WT_INPUT_R1.fastq.gz,WT_INPUT_R2.fastq.gz,1,,,
```

Verified 2026-11-22 (`docs/usage.md`, v2.1.0): seven required columns; `group` was replaced by `sample` in v2.0; `control_replicate` is mandatory alongside `control`; `antibody` required whenever `control` is set.

```bash
set -euo pipefail
nextflow run nf-core/chipseq -r 2.1.0 -profile docker --input samplesheet.csv --outdir results \
  --fasta GRCh38.fa --gtf GRCh38.gtf --blacklist hg38-blacklist.v2.bed --read_length 150 \
  --narrow_peak --aligner bowtie2
```

`--narrow_peak` exists (boolean; default is broad-peak mode with `--broad_cutoff 0.1`); `--macs_gsize` or `--read_length` is required to size the genome. v2.1.0 runs MACS3 internally. Add `--skip_spp` only if phantompeakqualtools is not wanted; `--save_macs_pileup` for bedGraphs.

## Validation

- Positive controls: an H3K4me3 or CTCF track must show the expected TSS / motif enrichment; `plotProfile` around TSS for H3K4me3 peaks at the TSS, H3K36me3 climbs across gene bodies, H3K27me3 anti-correlates with expression.
- Motif recovery: for a sequence-specific TF, the cognate JASPAR 2024 motif ranks in the top 5 of HOMER `knownResults.txt` with p < 1e-50 on the IDR set.
- Known loci: spot-check 3–5 literature-known targets in IGV (`bamCoverage` bigWigs + narrowPeak).
- Replicate correlation: `multiBigwigSummary bins` + `plotCorrelation --corMethod spearman` ≥ 0.8 between biological replicates for TF data.
- Peak set against ENCODE: `bedtools jaccard` with a matched ENCODE IDR set for the same factor/cell line; Jaccard ≥ 0.3 is typical for a reproduced experiment.
- Benchmark mechanics and truth-set handling: `.claude/skills/genomics-benchmark/SKILL.md`.

## Additional resources

- `.claude/skills/genomics-chipseq/atac.md` — ATAC-seq: chrM/MT fraction before dedup, Tn5 shifting, the two correct MACS3 modes, `macs3 hmmratac`, Genrich, TSS enrichment and fragment periodicity QC, NFR/mono split, TOBIAS footprinting, scATAC (ArchR, Signac, chromVAR, Cell Ranger ATAC/ARC), nf-core/atacseq. Load for any Tn5 dataset.
- `.claude/skills/genomics-chipseq/cutandrun.md` — CUT&RUN/CUT&Tag: Bowtie2 dovetail settings, E. coli spike-in scaling, SEACR vs MACS3, IgG handling, fragment-size expectations, nf-core/cutandrun. Load when there is no input, there is an IgG, or fragments are < 120 bp dominated.
- `.claude/skills/genomics-qc/SKILL.md` — raw-read QC, adapter trimming, MultiQC.
- `.claude/skills/genomics/reference.md` — reference genome, index and annotation acquisition (GRCh38, GRCm39, T2T, dm6, E. coli).
- `.claude/skills/genomics-workflow/SKILL.md` — nf-core/Nextflow execution, profiles, resume, SLURM.
- `.claude/skills/genomics-benchmark/SKILL.md` — truth-set and reproducibility benchmarking.
- `.claude/skills/genomics-rnaseq/SKILL.md` — integrating peaks with expression (TF targets vs DE genes).
