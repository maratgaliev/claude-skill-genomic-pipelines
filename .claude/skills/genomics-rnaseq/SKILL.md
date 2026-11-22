---
name: genomics-rnaseq
description: Bulk, single-cell and spatial RNA-seq - strandedness inference, STAR/Salmon/kallisto quantification, DESeq2/edgeR/limma DE, Seurat/Scanpy, pseudobulk, nf-core/rnaseq and nf-core/scrnaseq. Use when inputs are RNA-seq FASTQ, quant.sf, ReadsPerGene.out.tab, featureCounts tables, Cell Ranger/Space Ranger outs, .h5ad/.rds objects, or the user says "differential expression", "clustering", "Visium", "Xenium", "fusion", "splicing". Not for raw-read QC/trimming (genomics-qc), long-read isoform/direct-RNA analysis (genomics-longread), or scATAC/multiome chromatin (genomics-chipseq/atac.md).
argument-hint: "[bulk|sc|spatial] [path]"
paths:
  - "**/quant.sf"
  - "**/ReadsPerGene.out.tab"
  - "**/filtered_feature_bc_matrix/**"
  - "**/raw_feature_bc_matrix.h5"
  - "**/*.h5ad"
  - "**/spatial/tissue_positions*"
---

# RNA-seq

Expert in transcriptome quantification and inference: bulk RNA-seq (alignment or pseudoalignment → counts → DE → pathways), droplet single-cell (Cell Ranger/STARsolo/kb → QC → integration → annotation → pseudobulk DE) and spatial (Space Ranger/Xenium → Squidpy/Seurat). Outputs are commands, decision tables and validation checks, not narrative.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" STAR salmon kallisto featureCounts samtools infer_experiment.py cellranger spaceranger kb Rscript python3 nextflow || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| What | How | Decides |
|---|---|---|
| Layout, read length | `peek.sh R1.fq.gz`; `zcat R1.fq.gz \| head -4000 \| awk 'NR%4==2{print length($0)}' \| sort -n \| uniq -c` | `--sjdbOverhang` (= readlen−1), rMATS `--readLength`, 10x chemistry |
| 10x chemistry | R1 length **26** → 3' v2; **28** → 3' v3/v3.1/v4 and 5' v2/v3; R1 = cDNA and R2 ≈ 28 → check file labels swapped; `cellranger count` auto-detects | `--protocol 10XV2/10XV3`, kb `-x 10XV2/10XV3`, STARsolo `--soloUMIlen 10/12` |
| UMI presence | read names with `_UMI` / RX tag in BAM (`samtools view -c -d RX x.bam`); 10x R1 = 16 bc + 10/12 UMI | UMI dedup (`umi_tools`, `--with_umi` in nf-core) vs none |
| Strandedness (bulk) | **Before alignment:** `salmon quant -l A ...` on 1M reads → `lib_format_counts.json` `expected_format` (ISR/ISF/IU, SR/SF/U). **With a BAM:** `infer_experiment.py -r genes.bed12 -i x.bam` | `featureCounts -s`, Picard `STRAND_SPECIFICITY`, rMATS `--libType`, kallisto `--rf-stranded`, nf-core `strandedness` |
| Prior processing | BAM `@PG` lines (`samtools view -H x.bam \| grep ^@PG`): STAR version, `--quantMode`, MarkDuplicates | reuse vs realign |
| Reference / annotation | `samtools view -H x.bam \| grep ^@SQ \| head`; `peek.sh genes.gtf` → `chr1` vs `1`, GENCODE vs Ensembl, release number; `ref_sanity.sh` across FASTA/GTF/BAM | contig-name mismatch is the #1 cause of 0 % assigned reads |
| Sample design | samplesheet / `colData`: n per group, batch, paired subjects, covariates | DE method (table below), pseudobulk feasibility |
| Single-cell inputs | `peek.sh` on `metrics_summary.csv`, `web_summary.html` size, matrix dims (`zcat barcodes.tsv.gz \| wc -l`), raw vs filtered matrix present | CellBender needs **raw** matrix; demux needs HTO/feature matrix |

Strandedness mapping (library type → tool flag):

| Salmon `expected_format` | RSeQC | featureCounts | Picard | kallisto | nf-core |
|---|---|---|---|---|---|
| `IU` / `U` | ~50/50 | `-s 0` | `NONE` | (none) | `unstranded` |
| `ISF` / `SF` | `1++,1--,2+-,2-+` ≫ | `-s 1` | `SECOND_READ_TRANSCRIPTION_STRAND` | `--fr-stranded` | `forward` |
| `ISR` / `SR` (dUTP, TruSeq stranded, most kits) | `1+-,1-+,2++,2--` ≫ | `-s 2` | `FIRST_READ_TRANSCRIPTION_STRAND` | `--rf-stranded` | `reverse` |

Ask only if still ambiguous: contrast of interest and reference level; paired/blocked design; whether samples span batches not recorded in the samplesheet; nuclei vs whole cells (affects intronic counting and `percent.mt`); whether FFPE/probe-based (Visium CytAssist needs the probe set).

## Decision table

| Situation | Choice | Why |
|---|---|---|
| Bulk, gene-level DE, standard organism | **Salmon** (decoy-aware, `--gcBias --seqBias`) → `tximport` → DESeq2 | fastest, bias-aware, isoform-aware gene counts; no BAM to store |
| Need BAMs (splicing, fusions, variants, coverage QC, IGV) | **STAR** 2-pass + `featureCounts -s <inferred>` (or STAR `--quantMode GeneCounts`) | alignment required downstream; genome-level QC (RSeQC, Picard) |
| Very large cohort, laptop-class resources | **kallisto ≥0.50** (`kallisto quant`) | lowest RAM; no decoy, no bias model |
| Reproducible end-to-end, >6 samples, HPC | **nf-core/rnaseq** `--aligner star_salmon` (+ `--pseudo_aligner salmon` for comparison) | STAR BAM + Salmon quant + full QC in one run |
| Full-length, non-UMI plates (Smart-seq2/3) | nf-core/rnaseq per cell, then `scuttle`/Seurat | nf-core/scrnaseq is barcode-based only |
| Isoform / DTU | Salmon with `--numGibbsSamples 30` → `fishpond::swish` or `DRIMSeq`+`stageR`; DEXSeq for exon usage | transcript uncertainty must be propagated |
| Fusions | Arriba (STAR chimeric) or STAR-Fusion; nf-core/rnafusion for both + report | see `bulk.md` |
| Long-read cDNA / direct RNA | `.claude/skills/genomics-longread/SKILL.md` | different aligners/quantifiers |
| 10x droplet scRNA/snRNA | Cell Ranger (vendor-grade metrics) or STARsolo / kb-python (open, faster) → `singlecell.md` | |
| scATAC / multiome | `.claude/skills/genomics-chipseq/atac.md` | chromatin tooling |
| Visium / Visium HD / Xenium | Space Ranger / Xenium Ranger → Squidpy or Seurat; nf-core/spatialvi for Visium | `singlecell.md` spatial section |

DE method by design:

| n per group | Design | Method | Notes |
|---|---|---|---|
| 2–3 | simple two-group | DESeq2 (`lfcShrink(type="apeglm")`) | shrinkage essential; report shrunken LFC |
| ≥4 | any, incl. batch/blocked | edgeR quasi-likelihood (`glmQLFit`/`glmQLFTest`) | best FDR control, robust dispersion |
| ≥10 or large cohort | complex, random effects, many contrasts | limma-voom (`voomWithQualityWeights`, `duplicateCorrelation` for repeated measures) | fast, flexible linear modelling |
| any | gene set rather than gene | `fgsea` on ranked stat; `camera` (limma) | see `bulk.md` |
| single-cell, condition effect | **pseudobulk** per sample → edgeR/DESeq2 | cell-level Wilcoxon pseudoreplicates; see `singlecell.md` |

## Workflow

Minimal happy path (bulk, paired-end, human). Full variants, splicing, fusions, pathways: `bulk.md`.

1. Infer strandedness and build the decoy-aware Salmon index.

```bash
set -euo pipefail
GENOME=GRCh38.primary_assembly.genome.fa
TX=gencode.v50.transcripts.fa
IDX=salmon_idx_gencode50
THREADS=16
grep '^>' "${GENOME}" | cut -d ' ' -f1 | sed 's/>//' > decoys.txt
cat "${TX}" "${GENOME}" > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i "${IDX}" -k 31 -p "${THREADS}" --gencode
rm gentrome.fa
```

Validate: `test -s "${IDX}/versionInfo.json" && salmon --version`

2. Quantify every sample with `-l A`; read the inferred library type back.

```bash
set -euo pipefail
IDX=salmon_idx_gencode50
SAMPLE=S1
R1=fastq/${SAMPLE}_R1.fastq.gz
R2=fastq/${SAMPLE}_R2.fastq.gz
THREADS=8
salmon quant -i "${IDX}" -l A -1 "${R1}" -2 "${R2}" -p "${THREADS}" --gcBias --seqBias --numGibbsSamples 30 -o "quant/${SAMPLE}"
jq -r '.expected_format, .compatible_fragment_ratio' "quant/${SAMPLE}/lib_format_counts.json"
```

Validate: `jq '.percent_mapped' quant/*/aux_info/meta_info.json` (expect ≥70 %; all samples the same `expected_format`).

3. Import to gene level and run DESeq2 with an explicit reference level.

```r
suppressPackageStartupMessages({library(tximeta); library(DESeq2); library(apeglm)})
coldata <- read.csv("samplesheet.csv")            # columns: names, files, condition, batch
coldata$files <- file.path("quant", coldata$names, "quant.sf")
stopifnot(all(file.exists(coldata$files)))
se  <- tximeta(coldata)                            # GENCODE/Ensembl index is recognised automatically
gse <- summarizeToGene(se)
gse$condition <- relevel(factor(gse$condition), ref = "control")
dds <- DESeqDataSet(gse, design = ~ batch + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
print(resultsNames(dds))                           # confirm the coef name before shrinking
res <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")
write.csv(as.data.frame(res[order(res$padj), ]), "deseq2_treated_vs_control.csv")
```

Validate: `Rscript -e 'r<-read.csv("deseq2_treated_vs_control.csv",row.names=1); print(summary(r$padj<0.05)); print(range(r$log2FoldChange,na.rm=TRUE))'` — padj NA count small, LFC range not absurd (|LFC|>15 means a zero-count artefact survived filtering).

4. QC the design: PCA on `vst(dds, blind=TRUE)` coloured by condition and batch; a batch-dominated PC1 means batch must stay in the design (never regress it out of counts before DE).

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Reads per sample (bulk DE) | ≥20 M PE | 10–20 M acceptable for gene-level DE in well-expressed genes; splicing/fusion need ≥50 M |
| Salmon mapping rate | ≥70 % | <50 % → wrong species/annotation, heavy rRNA, or degraded FFPE; check `num_decoy_fragments` too |
| STAR uniquely mapped | ≥75 % | <60 % → contamination, adapter read-through, wrong genome; multimappers >20 % → rRNA/pseudogenes |
| Strandedness certainty (`compatible_fragment_ratio` or RSeQC dominant fraction) | ≥0.8 | 0.5–0.8 = mixed library/gDNA; 0.45–0.55 = unstranded |
| rRNA fraction | <10 % | poly-A: <5 %; ribo-depletion: 5–15 % tolerable; >20 % → failed depletion, drop or down-weight |
| Exonic rate (RSeQC/Picard `PCT_MRNA_BASES`) | ≥60 % poly-A | total RNA / nuclei: 30–50 % intronic is normal; count introns (`--quantMode` default includes them in STARsolo `GeneFull`) |
| Median 5'→3' bias (Picard `MEDIAN_5PRIME_TO_3PRIME_BIAS`) | 0.5–1.5 | FFPE/degraded skews 3'; use 3' kits' expectation; do not compare across kits |
| Duplication | <50 % unique-fragment duplication informative only with UMIs | high duplication in RNA-seq is expression, not artefact; never dedup bulk without UMIs |
| Genes detected (≥10 counts) | ≥12 000 (human tissue) | cell lines ~10 000; blood lower; sudden drop → degraded |
| Sample correlation (log counts, Spearman) | within-group ≥0.9 | outlier below group median−0.1 → inspect, do not silently drop |
| Sex check | XIST vs RPS4Y1/DDX3Y agree with metadata | mismatch → sample swap, stop before DE |

Single-cell thresholds are tissue-dependent and MAD-based: `singlecell.md`.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| featureCounts assigns ~50 % or near 0 % | wrong `-s` | infer strandedness first; `-s 2` for dUTP; never guess |
| featureCounts 0 % assigned, no error | GTF contig names (`1`) ≠ BAM (`chr1`) | `ref_sanity.sh`; use the GTF matching the FASTA used for the index |
| Salmon mapping 30–50 %, many `num_decoy_fragments` | genomic/intronic reads (nuclei, total RNA, gDNA) | expected; switch to STAR + intron-inclusive counting if intronic signal matters |
| `lfcShrink` error "coef not found" | `resultsNames(dds)` order differs from assumption; no `relevel` | print `resultsNames(dds)`; `relevel()` before `DESeq()` |
| DESeq2 LFC huge, padj NA | all-zero in one group, outlier replaced | keep count filter; inspect `cooksCutoff`; report shrunken LFC |
| PCA separates by batch, DE finds nothing | batch not in design | `~ batch + condition`; if batch confounded with condition, no fix exists - say so |
| `removeBatchEffect` matrix used for DE | visualisation-only transform fed to a model | DE on raw counts with batch covariate; ComBat-seq only when a covariate is impossible (and never on confounded designs) |
| tximport "transcripts missing from tx2gene" | GTF release ≠ transcript FASTA release, or version suffixes | same GENCODE release for both; `ignoreTxVersion = TRUE` / `ignoreAfterBar = TRUE`; prefer `tximeta` |
| STAR `--sjdbOverhang` warning / lower junction rate | index built for different read length | rebuild with readlen−1 (100 is fine for 101–150 bp) |
| STAR exits with "not enough memory" | human index needs ~32 GB RAM | `--genomeLoad LoadAndKeep` for many samples; or Salmon |
| nf-core/rnaseq `--skip_alignment` run still quantifies | `--skip_alignment` alone only drops STAR; Salmon pseudo-alignment still runs | QC-only run needs `--skip_alignment --skip_pseudo_alignment` |
| nf-core/rnaseq strandedness check "fail" | samplesheet value ≠ RSeQC/Salmon inference | set per-sample `strandedness` to the inferred value and `-resume` |
| 10x run: "chemistry could not be auto-detected" | R1/R2 swapped or trimmed R1 | check R1 length (26/28); never trim 10x R1 |
| Massive "cells" count, low UMI/cell | ambient / empty droplets called as cells | CellBender on raw matrix; inspect barcode-rank knee |
| Cluster = sample | no integration or true batch | Harmony/scVI, then verify by sample mixing per cluster |

## nf-core

Pipeline mechanics (`-profile`, `-resume`, SLURM, configs): `.claude/skills/genomics-workflow/SKILL.md`. Reference downloads: `.claude/skills/genomics/reference.md`.

**nf-core/rnaseq** (3.28.0, 2026-11). Samplesheet header: `sample,fastq_1,fastq_2,strandedness` (`strandedness` ∈ `unstranded|forward|reverse|auto`; `auto` → Salmon on 1 M reads, thresholds `--stranded_threshold 0.8 --unstranded_threshold 0.1`, cross-checked by RSeQC in MultiQC "Strandedness checks").

```bash
set -euo pipefail
FASTA=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
nextflow run nf-core/rnaseq -r 3.28.0 -profile singularity --input samplesheet.csv --outdir results --fasta "${FASTA}" --gtf "${GTF}" --gencode --aligner star_salmon --pseudo_aligner salmon --fail_on_undetermined_strandedness
```

Notes: `--aligner star_salmon` is the default (STAR BAM + Salmon quant of the transcriptome BAM); `--pseudo_aligner salmon` adds decoy-aware Salmon from FASTQ when `--fasta` is given. `--skip_alignment` alone switches to pseudo-alignment only (or reuses `genome_bam` columns) - it is **not** QC-only; QC-only is `--skip_alignment --skip_pseudo_alignment`. Salmon-only: `--skip_alignment --pseudo_aligner salmon`. UMIs: `--with_umi --umitools_extract_method regex --umitools_bc_pattern '...'`. Outputs for R: `results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv` or `results/star_salmon/tx2gene.tsv` + per-sample `quant.sf` → `tximport`. Downstream DE from these counts: **nf-core/differentialabundance** (2.0.0).

**nf-core/scrnaseq** (4.2.0). Header: `sample,fastq_1,fastq_2` (+ optional `expected_cells`). Aligners: `--aligner simpleaf` (default; piscem + alevin-fry), `star` (STARsolo), `kallisto` (kb), `cellranger`, `cellrangermulti` (+ `feature_type` column). `--protocol 10XV2|10XV3|10XV4|auto` (auto only for Cell Ranger). Cell hashing/genetic demux is not in scope (hadge pipeline). Smart-seq2 → nf-core/rnaseq.

**nf-core/rnafusion** (4.2.0). Header: `sample,fastq_1,fastq_2,strandedness` (`forward|reverse` required; or `bam`/`cram`/`junctions` columns). GRCh38 only; references must be downloaded first (`--references_only` or S3 sync) and rebuilt per release; `--tools arriba,starfusion,fusionreport,fusioninspector` (or `all`); `--no_cosmic` without a COSMIC licence.

**nf-core/spatialvi** (no tagged release as of 2026-11-22; `dev` only - pin `-r <commit>`). Raw header: `sample,fastq_dir,image,slide,area` (`cytaimage` instead of `image` for CytAssist; `--spaceranger_probeset` for FFPE). Processed header: `sample,spaceranger_dir`. Visium v1/v2/HD (`--hd_bin_size`, default 8 µm).

**nf-core/rnasplice** (2.0.0): rMATS, DEXSeq, SUPPA2, edgeR DEU from the same `sample,fastq_1,fastq_2,strandedness,condition` sheet.

**nf-core/smrnaseq** (2.4.1 as of 2026-11-22; verify current release in `https://nf-co.re/pipelines.json`). Small RNA / miRNA from single-end FASTQ. Header: `sample,fastq_1`. Library protocol is a config profile since 2.4 (`-profile illumina|qiaseq|nextflex|cats`, or set `--three_prime_adapter`/`--clip_r1` by hand); `--mirtrace_species hsa|mmu|...` is required; `--mirgenedb` switches miRBase → MirGeneDB; Bowtie1 quantification, `edger_qc/` count matrices, `mirtrace/` QC.

## Validation

- **Strandedness agreement:** Salmon `expected_format` (pre-alignment) and RSeQC `infer_experiment.py` (post-alignment) must name the same library type for every sample; record both in the report.
- **Count sanity:** library sizes within 3× of each other; `colSums(counts)` vs Salmon `num_mapped`; housekeeping genes (ACTB, GAPDH) in the top 1 % of every sample; sex genes agree with metadata.
- **Model sanity:** p-value histogram flat with a spike at 0 (U-shape → wrong dispersion/outliers; hill → unmodelled covariate); positive-control genes (known targets, spike-ins/ERCC if present) move as expected; `plotDispEsts` monotone decreasing.
- **Reproducibility:** rerun Salmon on one sample with a different `-p` and confirm `quant.sf` TPM identical to 4 decimals (Salmon 2.x documents thread-independent output; Gibbs samples are `-p`-independent since 2.12.0); DESeq2 result identical given fixed `dds` input.
- **Single-cell:** doublet rate ≈ 0.8 % per 1 000 cells loaded; marker genes recover canonical types before any annotation tool is trusted; pseudobulk PCA separates condition not sample index.
- **Truth sets / benchmarks:** SEQC/MAQC-III (sample A/B titrations, TaqMan truth), ERCC spike-in dose-response (slope ≈1, R² ≥0.9); mechanics in `.claude/skills/genomics-benchmark/SKILL.md`.

## Additional resources

- `bulk.md` - STAR/featureCounts with inferred `-s`, Salmon/kallisto/tximport/tx2gene, DESeq2/edgeR-QL/limma-voom code with batch, ComBat-seq vs covariate, splicing (rMATS-turbo, DEXSeq, SUPPA2), fusions (Arriba, STAR-Fusion), pathways (clusterProfiler, fgsea, msigdbr), RSeQC/Picard QC, visualisation. Load for any bulk analysis beyond the happy path.
- `singlecell.md` - Cell Ranger/STARsolo/kb, CellBender/SoupX, scDblFinder/Scrublet, MAD-based QC, Seurat v5 layers/IntegrateLayers and Scanpy/scvi-tools, SingleR/CellTypist/Azimuth, pseudobulk DE, trajectories, CITE-seq/hashing demux, spatial (Space Ranger, Visium HD, Xenium, Squidpy, Seurat spatial). Load for any droplet, plate or spatial dataset.
- `.claude/skills/genomics-qc/SKILL.md` - FastQC/fastp/MultiQC raw-read QC and trimming; load before alignment when FASTQ quality is unknown.
- `.claude/skills/genomics-longread/SKILL.md` - ONT/PacBio cDNA and direct-RNA (minimap2 `-ax splice`, IsoQuant, bambu, lr-kallisto, nf-core/nanoseq); load when read length >500 bp.
- `.claude/skills/genomics-chipseq/atac.md` - scATAC/multiome, peak calling, chromatin QC.
- `.claude/skills/genomics-workflow/SKILL.md` - Nextflow/nf-core execution, profiles, resume, HPC.
- `.claude/skills/genomics/reference.md` - GENCODE/Ensembl/10x reference acquisition and version pinning.
- `.claude/skills/genomics-benchmark/SKILL.md` - truth-set benchmarking mechanics.
