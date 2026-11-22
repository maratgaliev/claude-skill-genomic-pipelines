# ATAC-seq (bulk) and scATAC — deviations from the ChIP-seq core

Load after `SKILL.md`. Everything not listed here (alignment filters, blacklist, deepTools, DiffBind, annotation, motifs) is identical to the ChIP-seq workflow. ATAC has no input control; Tn5 cuts twice per fragment and inserts 9 bp staggered, so cut sites are read start +4 (plus strand) / read end −5 (minus strand).

## Inspect first (ATAC-specific)

| Read | Infer |
|---|---|
| `samtools idxstats raw.bam` chrM/MT row | MT fraction (must be computed on the **pre-dedup** BAM: MT reads are massively duplicated, so a post-dedup fraction understates the problem) |
| `samtools view -f 2 raw.bam \| cut -f9 \| head -100000` histogram | NFR peak < 100 bp, mono-nucleosome ~180–247 bp, di ~360 bp; absence of periodicity = over-tagmented or degraded |
| `@PG` lines | whether `alignmentSieve --ATACshift` or a shift awk was already applied (shifted BAMs must never be shifted again) |
| Read length / adapter content | Nextera (CTGTCTCTTATACACATCT) adapters indicate ATAC; see `genomics-qc` for trimming |
| Presence of `fragments.tsv.gz`, `filtered_peak_bc_matrix.h5`, `singlecell.csv` | 10x scATAC / multiome output → scATAC section |

Ask only if still ambiguous: Omni-ATAC vs standard protocol (affects MT expectation), nuclei vs whole cells, FACS-sorted input count (< 5k cells changes complexity expectations).

## Decision table

| Situation | Choice |
|---|---|
| Standard bulk PE ATAC, want accessible regions | Mode A: `macs3 callpeak -f BAMPE` on the **unshifted** dedup BAM (fragment pileup) — ENCODE-compatible and simplest |
| Want Tn5 cut-site resolution (footprinting, nucleosome-flank precision) | Mode B: shifted cut sites → BED → `macs3 callpeak -f BED --nomodel --shift -75 --extsize 150` |
| Want NFR + nucleosome states jointly | `macs3 hmmratac` on the unshifted PE BAM |
| Replicates, no control, want joint calling with chrM/dup handling built in | Genrich `-j` (ATAC mode) |
| TF footprints / differential TF activity | TOBIAS (ATACorrect → ScoreBigwig → BINDetect) on **unshifted** BAM (TOBIAS applies its own shift) |
| Differential accessibility | DiffBind as in `SKILL.md` §9, consensus from Mode A peaks, `DBA_NORM_LIB` + background |
| scATAC 10x | `cellranger-atac count` → ArchR or Signac; multiome → `cellranger-arc count` → Signac/ArchR |
| TSS enrichment | `ATACseqQC::TSSEscore` or ENCODE-style `computeMatrix reference-point` |

## Workflow

### 1. Align, keep MT reads for QC

```bash
set -euo pipefail
REF_IDX=/refs/GRCh38/bowtie2/GRCh38
SAMPLE=atac_rep1
THREADS=8
bowtie2 -p "${THREADS}" --very-sensitive -X 2000 -I 10 --no-mixed --no-discordant \
  -x "${REF_IDX}" -1 "${SAMPLE}_R1.trim.fq.gz" -2 "${SAMPLE}_R2.trim.fq.gz" \
  2> "${SAMPLE}.bowtie2.log" \
  | samtools sort -@ 4 -o "${SAMPLE}.raw.bam" - && samtools quickcheck "${SAMPLE}.raw.bam"
samtools index "${SAMPLE}.raw.bam"
```

No MAPQ/flag filter yet: MT fraction and NRF/PBC must be measured on everything. `-X 2000` keeps di-/tri-nucleosome fragments. Validation: `grep overall "${SAMPLE}.bowtie2.log"`.

### 2. MT fraction, then remove chrM, then dedup (order matters)

```bash
set -euo pipefail
SAMPLE=atac_rep1
MT=chrM
samtools idxstats "${SAMPLE}.raw.bam" \
  | awk -v mt="${MT}" '{t+=$3} $1==mt{m=$3} END{printf "MT_fraction\t%.4f\n", m/t}' \
  | tee "${SAMPLE}.mt_fraction.tsv"
samtools view -b -q 30 -F 1804 -f 2 -e 'rname!="chrM" && mrname!="chrM"' \
  "${SAMPLE}.raw.bam" -o "${SAMPLE}.noMT.bam"
samtools index "${SAMPLE}.noMT.bam"
```

Ensembl-style references name it `MT`: set `MT=MT` and edit the `-e` expression; or build a region list from `idxstats` (`samtools idxstats raw.bam | cut -f1 | grep -v -e '^chrM$' -e '^MT$' -e '^\*$' > keep.txt; samtools view -b raw.bam $(cat keep.txt)`), which also drops unplaced contigs if you grep them out. Never `grep -v chrM` on SAM text: it deletes any read whose mate/aux text contains the string and corrupts headers. Then dedup exactly as `SKILL.md` §2 (`gatk MarkDuplicates` or `fixmate -m` + `markdup`) producing `${SAMPLE}.dedup.bam`, then blacklist (§3) producing `${SAMPLE}.filt.bam`. Validation: `cat "${SAMPLE}.mt_fraction.tsv"`; `samtools idxstats "${SAMPLE}.noMT.bam" | grep -c chrM` prints 0.

### 3. Fragment-size periodicity

```bash
set -euo pipefail
SAMPLE=atac_rep1
samtools view -f 66 "${SAMPLE}.filt.bam" | awk '{l=$9<0?-$9:$9; if(l<=1000) h[l]++} END{for(i in h) print i"\t"h[i]}' \
  | sort -k1,1n > "${SAMPLE}.fraglen.tsv"
bamPEFragmentSize -b "${SAMPLE}.filt.bam" --histogram "${SAMPLE}.fraglen.png" --maxFragmentLength 1000 -p 8
```

Expect a dominant < 100 bp mode, a clear ~200 bp mode, a visible ~400 bp shoulder, and 10.4 bp helical ripple in the sub-100 bp region. Validation: `awk '$1<100{n+=$2} $1>=180&&$1<=247{m+=$2} END{print "NFR/mono", n/m}' "${SAMPLE}.fraglen.tsv"` (> 1 expected for Omni-ATAC; ENCODE only requires both modes to be present).

### 4. Tn5 shift (only when cut-site resolution is needed)

```bash
set -euo pipefail
SAMPLE=atac_rep1
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
alignmentSieve -b "${SAMPLE}.filt.bam" -o "${SAMPLE}.shifted.unsorted.bam" --ATACshift \
  --blackListFileName "${BLACKLIST}" --minMappingQuality 30 -p 8
samtools sort -@ 4 -o "${SAMPLE}.shifted.bam" "${SAMPLE}.shifted.unsorted.bam"
samtools index "${SAMPLE}.shifted.bam"
rm "${SAMPLE}.shifted.unsorted.bam"
```

`alignmentSieve --ATACshift` is the only sanctioned shifter (the old awk `$4±4/5` edits POS without fixing CIGAR, PNEXT or TLEN and breaks `BAMPE`). Its output is unsorted; sort and index before any downstream tool. Validation: `samtools view -H "${SAMPLE}.shifted.bam" | grep -c alignmentSieve` prints 1 — if it prints 2 the BAM was shifted twice; discard it.

### 5. Peak calling — pick exactly one mode

**Mode A (fragment-based, ENCODE-style, default).** Unshifted `filt.bam`; `BAMPE` uses the whole fragment so `--shift`/`--extsize` must not be given (MACS3 ignores `--shift` with `BAMPE` and the combination documents a misunderstanding).

```bash
set -euo pipefail
SAMPLE=atac_rep1
GSIZE=hs
macs3 callpeak -t "${SAMPLE}.filt.bam" -f BAMPE -g "${GSIZE}" -n "${SAMPLE}.modeA" --outdir peaks \
  -q 0.01 --keep-dup all -B --SPMR
```

**Mode B (cut-site based).** Convert the shifted BAM to single-end cut-site BED (each read start is a Tn5 insertion), then centre a 150 bp window on each cut (`--shift -75 --extsize 150`). Never run Mode B on a BAM with `-f BAMPE`, and never apply `--shift` to reads that `alignmentSieve --ATACshift` already moved **and** then also feed `BAMPE`: the shift belongs to the tag-based call only.

```bash
set -euo pipefail
SAMPLE=atac_rep1
GSIZE=hs
bedtools bamtobed -i "${SAMPLE}.shifted.bam" \
  | awk 'BEGIN{OFS="\t"} $6=="+"{print $1,$2,$2+1,$4,$5,$6} $6=="-"{print $1,$3-1,$3,$4,$5,$6}' \
  | sort -k1,1 -k2,2n > "${SAMPLE}.cutsites.bed"
macs3 callpeak -t "${SAMPLE}.cutsites.bed" -f BED -g "${GSIZE}" -n "${SAMPLE}.modeB" --outdir peaks \
  --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.01 -B --SPMR
```

If you skipped step 4, Mode B still works from the unshifted BAM by adding the +4/−5 offsets in the awk (`$2+4` for `+`, `$3-5` for `-`); do not do both. Mode A and Mode B peaks overlap > 85 % on good libraries; a large discrepancy means the shift was applied twice or `-f` was wrong. Validation: `wc -l peaks/${SAMPLE}.mode*_peaks.narrowPeak`.

**HMMRATAC (NFR + nucleosome states).** Unshifted PE BAM; run the cutoff analysis first.

```bash
set -euo pipefail
SAMPLE=atac_rep1
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
macs3 hmmratac -i "${SAMPLE}.filt.bam" -f BAMPE -n "${SAMPLE}" --outdir hmmratac -e "${BLACKLIST}" --cutoff-analysis-only
macs3 hmmratac -i "${SAMPLE}.filt.bam" -f BAMPE -n "${SAMPLE}" --outdir hmmratac -e "${BLACKLIST}" -l 10 -u 20 -c 1.2
```

Set `-l/-u/-c` from `hmmratac/${SAMPLE}_cutoff_analysis.tsv` (lower ≈ 10k peaks of 500–1000 bp; upper captures a few hundred extreme regions). The standalone HMMRATAC jar is legacy; `macs3 hmmratac` (3.1.0) also accepts `-f FRAG` for scATAC fragment files. Validation: `wc -l hmmratac/${SAMPLE}_accessible_regions.narrowPeak`.

**Genrich (joint replicates, built-in chrM/dup handling).** Needs name-sorted BAMs.

```bash
set -euo pipefail
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
samtools sort -n -@ 4 -o rep1.nsort.bam atac_rep1.noMT.bam
samtools sort -n -@ 4 -o rep2.nsort.bam atac_rep2.noMT.bam
Genrich -t rep1.nsort.bam,rep2.nsort.bam -o genrich/joint.narrowPeak -j -r -e chrM -E "${BLACKLIST}" -q 0.05 -v
```

`-j` ATAC mode (centres 100 bp on each cut), `-r` remove PCR duplicates, `-e` exclude chromosomes, `-E` exclude regions. Feed the pre-dedup BAM because Genrich dedups itself. Validation: `wc -l genrich/joint.narrowPeak`.

### 6. TSS enrichment and ATAC QC in R

```r
library(ATACseqQC); library(GenomicAlignments); library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
bamfile <- "atac_rep1.filt.bam"
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
# fragment-size periodicity (uses the BAM path)
fragSizeDist(bamfile, "atac_rep1")
# TSSEscore needs a GAlignments object, not a path
ga <- readBamFile(bamfile, tag = character(0), asMates = TRUE, bigFile = TRUE)
# equivalent: ga <- readGAlignments(bamfile, param = ScanBamParam(flag = scanBamFlag(isProperPair = TRUE)))
tsse <- TSSEscore(ga, txs)
tsse$TSSEscore
plot(100 * (-9:10 - .5), tsse$values, type = "b", xlab = "distance to TSS", ylab = "aggregate TSS score")
```

ENCODE's definition uses ±2000 bp aggregate read depth normalised to the 100 bp end flanks; `TSSEscore` implements the same shape. Alternative without R: `computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 --binSize 100` on a CPM bigWig, then centre / mean of the two outer bins. Validation: value recorded in the QC table below; compare against the matching annotation row.

### 7. NFR / mono-nucleosome split by TLEN

```bash
set -euo pipefail
SAMPLE=atac_rep1
samtools view -b -e 'abs(tlen)<100' "${SAMPLE}.filt.bam" -o "${SAMPLE}.nfr.bam"
samtools view -b -e 'abs(tlen)>=180 && abs(tlen)<=247' "${SAMPLE}.filt.bam" -o "${SAMPLE}.mono.bam"
samtools index "${SAMPLE}.nfr.bam"; samtools index "${SAMPLE}.mono.bam"
```

Equivalent via deepTools: `alignmentSieve --maxFragmentLength 100` / `--minFragmentLength 180 --maxFragmentLength 247`. Call TF-centred peaks on `nfr.bam` (Mode A), nucleosome positions with `macs3 hmmratac` or NucleoATAC on `mono.bam`. Validation: `samtools view -c "${SAMPLE}.nfr.bam"` > `samtools view -c "${SAMPLE}.mono.bam"` for a healthy library.

### 8. Footprinting — TOBIAS 0.17.5

```bash
set -euo pipefail
REF=/refs/GRCh38/GRCh38.fa
BLACKLIST=/refs/blacklists/hg38-blacklist.v2.bed
JASPAR=/refs/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
cat peaks/*modeA_peaks.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > tobias/merged_peaks.bed
for S in cond1 cond2; do
  TOBIAS ATACorrect --bam "${S}.filt.bam" --genome "${REF}" --peaks tobias/merged_peaks.bed \
    --blacklist "${BLACKLIST}" --outdir "tobias/${S}" --cores 8
  TOBIAS ScoreBigwig --signal "tobias/${S}/${S}.filt_corrected.bw" --regions tobias/merged_peaks.bed \
    --output "tobias/${S}/${S}_footprints.bw" --cores 8
done
TOBIAS BINDetect --motifs "${JASPAR}" --signals tobias/cond1/cond1_footprints.bw tobias/cond2/cond2_footprints.bw \
  --genome "${REF}" --peaks tobias/merged_peaks.bed --outdir tobias/bindetect --cond_names cond1 cond2 --cores 8
```

Feed the **unshifted** `filt.bam`: ATACorrect models and corrects Tn5 insertion bias itself and expects raw read ends. `BINDetect` accepts JASPAR, MEME or PFM motif files; use a merged peak set across conditions so footprint scores are comparable. Validation: `head tobias/bindetect/bindetect_results.txt` sorted by `cond1_cond2_change`; known lineage TFs should appear at the extremes.

## scATAC-seq

| Step | Command / package | Notes |
|---|---|---|
| 10x scATAC | `cellranger-atac count --id S --reference refdata-cellranger-arc-GRCh38-2024-A --fastqs fq/ --sample S` | Cell Ranger ATAC 2.2 (2025-04-30, release notes verified); outputs `fragments.tsv.gz`, `singlecell.csv`, `filtered_peak_bc_matrix.h5`. Reference bundle name: `<check current: https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads>` |
| 10x multiome | `cellranger-arc count --id S --reference <arc ref> --libraries libraries.csv` | `libraries.csv` columns `fastqs,sample,library_type` with `Chromatin Accessibility` and `Gene Expression` rows |
| ArchR 1.0.3 | `createArrowFiles(inputFiles = frags, sampleNames, minTSS = 4, minFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)` → `ArchRProject` → `addIterativeLSI` → `addClusters` → `addGroupCoverages` → `addReproduciblePeakSet` (calls `macs3`/`macs2`; set `pathToMacs2 <- findMacs2()`) → `addPeakMatrix` → `addMotifAnnotations(motifSet = "JASPAR2024")` → `addDeviationsMatrix` (chromVAR) | Arrow files keep everything on disk; fastest for > 50k cells |
| Signac (Seurat ecosystem) | `CreateChromatinAssay(counts, fragments = "fragments.tsv.gz", sep = c(":","-"))` → `NucleosomeSignal()` → `TSSEnrichment()` → filter `TSS.enrichment > 3`, `nucleosome_signal < 4` → `RunTFIDF` → `FindTopFeatures` → `RunSVD` (drop LSI component 1 if it tracks depth) → `CallPeaks(macs2.path = Sys.which("macs3"))` → `AddMotifs(pfm = getMatrixSet(JASPAR2024, list(species = 9606)))` → `RunChromVAR` | best for multiome co-embedding with RNA (`FindMultiModalNeighbors`) |
| chromVAR 1.34 | `getBackgroundPeaks` → `computeDeviations` → `computeVariability` | GC-matched background peaks; needs `addGCBias` first |
| Peak calling per cluster | `macs3 hmmratac -i fragments.tsv.gz -f FRAG --barcodes cluster1_barcodes.txt --max-count 2` or `macs3 callpeak -f BEDPE` on per-cluster fragment subsets | cluster-wise calling recovers rare-population peaks that a pooled call misses |
| Doublets | `ArchR::addDoubletScores` / `scDblFinder` (`aggregateFeatures = TRUE`) / AMULET | AMULET uses multiplet overlap counts on fragments, independent of clustering |

Signac version not pinned here: `<check current: https://stuartlab.org/signac/>`. scATAC QC gates (per cell): TSS enrichment ≥ 4 (ArchR default), fragments ≥ 1000, FRiP ≥ 0.15–0.2 (10x-reported `peak_region_fragments / passed_filters`), nucleosome signal < 4, blacklist fraction < 0.05.

## QC thresholds

ENCODE ATAC-seq standards (`https://www.encodeproject.org/atac-seq/`, updated July 2020) unless noted.

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Non-duplicate, non-MT aligned reads | ≥ 50 M PE reads (= 25 M fragments) per replicate; 25 M for SE | ENCODE; discovery of rare elements scales with depth |
| Alignment rate | > 95 % ideal; > 80 % acceptable | ENCODE |
| MT fraction (pre-dedup) | < 10 % Omni-ATAC; < 20–30 % standard protocol; > 50 % concerning | not an ENCODE gate; Omni-ATAC (Corces 2017) typically 1–5 %. Compute before dedup |
| NRF / PBC1 / PBC2 | > 0.9 / > 0.9 / > 3 ideal; 0.7–0.9 / 0.7–0.9 / 1–3 acceptable; < 0.7 / < 0.7 / < 1 concerning | ENCODE ATAC library-complexity table (note PBC2 > 3, not > 10 as for ChIP) |
| TSS enrichment, GRCh38 RefSeq | < 5 concerning; 5–7 acceptable; > 7 ideal | ENCODE; hg19 RefSeq: < 6 / 6–10 / > 10; mm10 RefSeq: < 10 / 10–15 / > 15; GRCm39 not tabulated by ENCODE — use the mm10 row and say so |
| FRiP | > 0.3 ideal; > 0.2 acceptable | ENCODE; waived for EN-TEx tissues; depends on peak caller — report caller and settings |
| Replicated peaks | > 150 000 ideal; > 100 000 acceptable | ENCODE (naive overlap on pooled vs each replicate ≥ 50 %) |
| IDR peaks | > 70 000 ideal; > 50 000 acceptable | ENCODE; IDR ≤ 0.05 |
| IDR ratios | self-consistency and rescue both < 2 | ENCODE |
| Fragment-size distribution | NFR mode present and mono-nucleosome peak (147–294 bp) present | ENCODE; missing periodicity = over-tagmentation / dead cells |
| NFR / mono ratio | > 1 typical; report value | not an ENCODE standard; protocol-dependent |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| Peaks all 150 bp wide and offset 75 bp from fragment centres | `--shift -75 --extsize 150` applied with `-f BAMPE`, or applied to an already-shifted BAM | Mode A (BAMPE, no shift) or Mode B (BED cut sites, shift once) — never both |
| `alignmentSieve` output rejected by `samtools index` | unsorted output | `samtools sort` first (step 4) |
| MT fraction reported as 2 % but library looked poor | fraction computed after dedup | recompute on `raw.bam` with `idxstats` |
| `samtools view -e 'rname!="chrM"'` leaves mates on chrM | expression lacked `mrname` clause, or contig is `MT` | include `mrname!="chrM"`; check `idxstats` names |
| `grep -v chrM` produced a truncated/invalid BAM | text grep removed header lines or reads with `chrM` in tags | use `-e` expressions or region lists |
| `TSSEscore` error "no slot of name" / "GAlignments expected" | a file path or GRanges was passed | `readBamFile`/`readGAlignments` first (step 6) |
| TSS score 3 on mouse data that looks fine in IGV | compared against hg38 row; mm10 thresholds are higher (10/15) but computed on RefSeq mm10 TSS | use the matching annotation row and the same TSS set |
| HMMRATAC runs for hours then zero peaks | `-l/-u/-c` defaults unsuited to depth | run `--cutoff-analysis-only` first; set from report |
| TOBIAS footprints look inverted | shifted BAM supplied to ATACorrect | supply unshifted `filt.bam` |
| Genrich "not paired" / wrong duplicate counts | position-sorted BAM | name-sort (`samtools sort -n`) |
| Signac LSI component 1 correlates with depth | expected | `DepthCor`; drop dim 1 in `RunUMAP`/`FindNeighbors` |
| ArchR `addReproduciblePeakSet` fails to find MACS | MACS3 binary not named macs2 | `pathToMacs2 = Sys.which("macs3")` |

## nf-core/atacseq (2.1.2)

```csv
sample,fastq_1,fastq_2,replicate
CONTROL,CONTROL_R1.fastq.gz,CONTROL_R2.fastq.gz,1
TREAT,TREAT_R1.fastq.gz,TREAT_R2.fastq.gz,1
```

Verified 2026-11-22 (`https://nf-co.re/atacseq/usage`): four required columns; with `--with_control` add `control,control_replicate`. `--read_length` (or `--macs_gsize`) is mandatory; `--mito_name chrM` (or `MT`) drives the MT filter; `--narrow_peak` switches MACS from the default broad mode; `--skip_merge_replicates` keeps replicate-level calls for IDR.

```bash
set -euo pipefail
nextflow run nf-core/atacseq -r 2.1.2 -profile docker --input samplesheet.csv --outdir results \
  --fasta GRCh38.fa --gtf GRCh38.gtf --blacklist hg38-blacklist.v2.bed --read_length 150 \
  --mito_name chrM --narrow_peak --aligner bowtie2
```

Generic Nextflow mechanics: `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

- TSS enrichment and fragment periodicity recomputed from the final BAM match the QC table; both plots saved with the sample name.
- Mode A vs Mode B peak overlap > 85 % (`bedtools intersect -u | wc -l` over the smaller set); CTCF motif centred within ±10 bp of footprint centres in `TOBIAS BINDetect` output.
- Housekeeping promoters (GAPDH, ACTB) are accessible in every sample; lineage promoters follow expected pattern (e.g. PTPRC in immune samples).
- Replicate Spearman correlation of counts in consensus peaks ≥ 0.9.
- Benchmark set: ENCODE ATAC IDR peaks for the same cell type (`bedtools jaccard` ≥ 0.4 typical); mechanics in `.claude/skills/genomics-benchmark/SKILL.md`.
