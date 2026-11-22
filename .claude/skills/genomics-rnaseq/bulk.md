# Bulk RNA-seq: quantification, DE, splicing, fusions, pathways

Supporting file for `.claude/skills/genomics-rnaseq/SKILL.md`. Raw-read QC/trimming: `.claude/skills/genomics-qc/SKILL.md`. Reference downloads: `.claude/skills/genomics/reference.md`. Every strand-dependent flag below is derived from the inferred library type (SKILL.md "Inspect first" mapping table), never from the kit name alone.

Versions verified 2026-11-22: STAR 2.7.11b; Salmon 2.12.0 (Rust rewrite; C++ 1.12.0 is `salmon-cpp`, indices are not interchangeable); kallisto 0.52.0; nf-core/rnaseq 3.28.0; nf-core/rnafusion 4.2.0; nf-core/rnasplice 2.0.0.

## 1. Strandedness inference (do this first)

```bash
set -euo pipefail
IDX=salmon_idx_gencode50
R1=fastq/S1_R1.fastq.gz
R2=fastq/S1_R2.fastq.gz
salmon quant -i "${IDX}" -l A -1 <(zcat "${R1}" | head -4000000) -2 <(zcat "${R2}" | head -4000000) -p 8 -o strand_check/S1
jq '{expected_format, compatible_fragment_ratio, num_compatible_fragments, num_assigned_fragments}' strand_check/S1/lib_format_counts.json
```

Validate: `compatible_fragment_ratio ≥ 0.8`; `expected_format` one of `IU ISF ISR` (PE) / `U SF SR` (SE). Convert to tool flags with the mapping table in SKILL.md and export once for the whole run:

```bash
set -euo pipefail
LIBTYPE=$(jq -r .expected_format strand_check/S1/lib_format_counts.json)
case "${LIBTYPE}" in
  IU|U)   FC_S=0; PICARD_STRAND=NONE;                             KALLISTO_STRAND=""; RMATS_LIB=fr-unstranded ;;
  ISF|SF) FC_S=1; PICARD_STRAND=SECOND_READ_TRANSCRIPTION_STRAND; KALLISTO_STRAND="--fr-stranded"; RMATS_LIB=fr-secondstrand ;;
  ISR|SR) FC_S=2; PICARD_STRAND=FIRST_READ_TRANSCRIPTION_STRAND;  KALLISTO_STRAND="--rf-stranded"; RMATS_LIB=fr-firststrand ;;
  *) echo "unexpected libtype ${LIBTYPE}" >&2; exit 1 ;;
esac
mkdir -p strand_check
printf 'LIBTYPE=%s\nFC_S=%s\nPICARD_STRAND=%s\nKALLISTO_STRAND=%q\nRMATS_LIB=%s\n' "${LIBTYPE}" "${FC_S}" "${PICARD_STRAND}" "${KALLISTO_STRAND}" "${RMATS_LIB}" | tee strand_check/strand.env
```

With an existing BAM instead of FASTQ (RSeQC):

```bash
set -euo pipefail
BAM=S1.Aligned.sortedByCoord.out.bam
BED12=gencode.v50.bed12
infer_experiment.py -r "${BED12}" -i "${BAM}" -s 400000 > strand_check/S1.infer_experiment.txt
cat strand_check/S1.infer_experiment.txt
```

`"1++,1--,2+-,2-+"` dominant → `-s 1`; `"1+-,1-+,2++,2--"` dominant → `-s 2`; ≈0.5/0.5 → `-s 0`. BED12 from GTF: `gtfToGenePred genes.gtf stdout | genePredToBed stdin genes.bed12` (UCSC tools) or `paftools.js gff2bed` (minimap2).

## 2. STAR + featureCounts

```bash
set -euo pipefail
GENOME=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
STAR_IDX=star_idx_gencode50_oh100
READLEN=101
THREADS=16
STAR --runMode genomeGenerate --runThreadN "${THREADS}" --genomeDir "${STAR_IDX}" --genomeFastaFiles "${GENOME}" --sjdbGTFfile "${GTF}" --sjdbOverhang $((READLEN - 1))
```

Validate: `test -s "${STAR_IDX}/SA" && grep -E 'sjdbOverhang|versionGenome' "${STAR_IDX}/genomeParameters.txt"` (human index ≈ 30 GB RAM, ~1 h).

```bash
set -euo pipefail
STAR_IDX=star_idx_gencode50_oh100
SAMPLE=S1
R1=fastq/${SAMPLE}_R1.fastq.gz
R2=fastq/${SAMPLE}_R2.fastq.gz
THREADS=16
STAR --runThreadN "${THREADS}" --genomeDir "${STAR_IDX}" --readFilesIn "${R1}" "${R2}" --readFilesCommand zcat --outFileNamePrefix "star/${SAMPLE}." --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI AS nM NM MD XS --outSAMstrandField intronMotif --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMax 1000000 --alignMatesGapMax 1000000
samtools index "star/${SAMPLE}.Aligned.sortedByCoord.out.bam"
samtools quickcheck "star/${SAMPLE}.Aligned.sortedByCoord.out.bam"
```

Validate: `grep -E 'Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped: too short' star/${SAMPLE}.Log.final.out`. ENCODE-style filters above; drop `--outFilterMismatchNoverReadLmax` for highly polymorphic samples. `--quantMode TranscriptomeSAM` yields `Aligned.toTranscriptome.out.bam` for `salmon quant -t tx.fa -l A -a` or RSEM.

STAR's own counts: `ReadsPerGene.out.tab` columns 2/3/4 = unstranded / `-s 1` / `-s 2`; pick the column by inferred strandedness, same rule as featureCounts. Never pick "the column with more counts" - antisense counting also gives a large number.

```bash
set -euo pipefail
GTF=gencode.v50.primary_assembly.annotation.gtf
source strand_check/strand.env
featureCounts -T 8 -p --countReadPairs -s "${FC_S}" -t exon -g gene_id -a "${GTF}" -o counts/featureCounts.txt star/*.Aligned.sortedByCoord.out.bam
cut -f1,7- counts/featureCounts.txt | tail -n +2 > counts/gene_counts.tsv
```

Validate: `awk 'NR>1 && /Assigned|Unassigned_NoFeatures|Unassigned_Ambiguity/' counts/featureCounts.txt.summary` → Assigned ≥60 % of PE fragments; NoFeatures >40 % means wrong `-s`, wrong GTF build, or intron-rich (nuclei/total RNA) library. Single-end: drop `-p --countReadPairs`. Options worth knowing: `-O` (count multi-overlap), `-M --fraction` (multimappers, usually off for DE), `--extraAttributes gene_name,gene_type`, `-t gene` for intron-inclusive counting of nuclei.

## 3. Salmon (decoy-aware) and kallisto

Index (full genome as decoy, GENCODE headers trimmed with `--gencode`): SKILL.md Workflow step 1. Salmon 2.x removed `--validateMappings` (selective alignment is the only mode; the flag is accepted and ignored) and removed `salmon alevin` (use alevin-fry/simpleaf). Quantify:

```bash
set -euo pipefail
IDX=salmon_idx_gencode50
SAMPLE=S1
R1=fastq/${SAMPLE}_R1.fastq.gz
R2=fastq/${SAMPLE}_R2.fastq.gz
salmon quant -i "${IDX}" -l A -1 "${R1}" -2 "${R2}" -p 8 --gcBias --seqBias --posBias --numGibbsSamples 30 -o "quant/${SAMPLE}"
```

Validate: `jq '{percent_mapped, num_decoy_fragments, library_types}' quant/${SAMPLE}/aux_info/meta_info.json`. `--gcBias --seqBias` are recommended for every Illumina library (cost ~15 %); `--posBias` for degraded/3'-biased samples. Single-end: `-r "${R1}"`. Ensembl FASTA + Ensembl GTF: add `--ignoreTxVersion` or `ignoreTxVersion = TRUE` in tximport. Alignment mode from a STAR transcriptome BAM: `salmon quant -t "${TX}" -l A -a star/S1.Aligned.toTranscriptome.out.bam -o quant/S1 --gcBias --seqBias`.

kallisto ≥0.50 (new index format; rebuild indices from 0.4x):

```bash
set -euo pipefail
TX=gencode.v50.transcripts.fa
KIDX=kallisto_idx_gencode50.idx
SAMPLE=S1
R1=fastq/${SAMPLE}_R1.fastq.gz
R2=fastq/${SAMPLE}_R2.fastq.gz
source strand_check/strand.env
kallisto index -i "${KIDX}" -t 8 "${TX}"
kallisto quant -i "${KIDX}" -o "kallisto/${SAMPLE}" -t 8 -b 30 ${KALLISTO_STRAND} "${R1}" "${R2}"
```

Validate: `jq '{n_processed, n_pseudoaligned, p_pseudoaligned}' kallisto/${SAMPLE}/run_info.json`. Single-end requires `--single -l 200 -s 20` (fragment length is not estimable). kallisto has no decoy or bias model; prefer Salmon when a genome is available.

## 4. tx2gene and import

Preferred: `tximeta` (SKILL.md step 3) - it recognises GENCODE/Ensembl indices by checksum and attaches gene ranges. Manual tx2gene when tximeta cannot identify the index (custom transcriptome):

```r
suppressPackageStartupMessages({library(GenomicFeatures); library(txdbmaker); library(tximport); library(readr)})
gtf <- "gencode.v50.primary_assembly.annotation.gtf"
txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
tx2gene <- AnnotationDbi::select(txdb, keys = keys(txdb, "TXNAME"), keytype = "TXNAME", columns = "GENEID")
write_tsv(tx2gene, "tx2gene.tsv")
samples <- read.csv("samplesheet.csv")              # names, condition, batch
files <- setNames(file.path("quant", samples$names, "quant.sf"), samples$names)
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE, countsFromAbundance = "no")
# kallisto: type = "kallisto", files = .../abundance.h5 ; use ignoreTxVersion = TRUE for Ensembl
```

Validate: `length(setdiff(rownames(read.delim(files[1]))[1:5], tx2gene$TXNAME)) == 0` after bar/version handling; `txi$counts` row count ≈ number of genes in GTF. `countsFromAbundance = "lengthScaledTPM"` only when feeding a method that cannot take offsets (e.g. limma without `voom` offsets, exported matrices); DESeq2 uses `DESeqDataSetFromTximport(txi, colData, design)` which applies length offsets itself. Python equivalent: `pytximport`.

## 5. Differential expression

Design-first rules: reference level set explicitly; batch and known nuisance covariates in the formula; paired subjects as a blocking factor (`~ subject + condition`); interaction only when the question is an interaction. Filter genes before testing (`filterByExpr` or ≥10 counts in ≥ smallest group size).

### DESeq2 (any n; standard for n ≤ 3)

```r
suppressPackageStartupMessages({library(DESeq2); library(apeglm)})
coldata <- read.csv("samplesheet.csv", row.names = "names")
coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch <- factor(coldata$batch)
dds <- DESeqDataSetFromTximport(txi, colData = coldata[colnames(txi$counts), ], design = ~ batch + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= min(table(coldata$condition)), ]
dds <- DESeq(dds)
rn <- resultsNames(dds); print(rn)                  # e.g. "condition_treated_vs_control"
coef <- grep("^condition_.*_vs_control$", rn, value = TRUE)[1]
res <- lfcShrink(dds, coef = coef, type = "apeglm")           # apeglm needs coef, not contrast
res_unshrunk <- results(dds, name = coef, alpha = 0.05)       # for p-values / summary()
summary(res_unshrunk)
# arbitrary contrast (non-reference pair): results(dds, contrast = c("condition","B","A")); shrink with type = "ashr"
out <- as.data.frame(res); out$padj <- res_unshrunk$padj[match(rownames(out), rownames(res_unshrunk))]
write.csv(out[order(out$padj), ], "deseq2_results.csv")
```

Validate: `hist(res_unshrunk$pvalue, breaks = 50)` flat with spike at 0; `plotDispEsts(dds)`; `sum(is.na(res$padj))` small. Continuous covariate (RIN, age): centre it and keep it in the design. Multi-factor time course: `DESeq(dds, test = "LRT", reduced = ~ batch + time)`.

### edgeR quasi-likelihood (n ≥ 4, blocked / batch designs)

```r
suppressPackageStartupMessages(library(edgeR))
cts <- txi$counts                                         # or read featureCounts matrix
coldata <- read.csv("samplesheet.csv", row.names = "names")[colnames(cts), ]
condition <- relevel(factor(coldata$condition), ref = "control")
batch <- factor(coldata$batch)
design <- model.matrix(~ batch + condition)
y <- DGEList(cts, group = condition)
y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)                                      # edgeR ≥4 name for calcNormFactors (TMM)
fit <- glmQLFit(y, design, robust = TRUE)                 # edgeR ≥4: estimateDisp no longer required first
qlf <- glmQLFTest(fit, coef = "conditiontreated")
tt  <- topTags(qlf, n = Inf)$table
write.csv(tt, "edger_ql_results.csv")
# effect-size threshold test instead of post-hoc |logFC| filter:
tr <- glmTreat(fit, coef = "conditiontreated", lfc = log2(1.5)); topTags(tr)
```

Validate: `plotQLDisp(fit)`; `summary(decideTests(qlf))`; `plotMDS(y, col = as.integer(batch))` shows batch structure you modelled. For tximport input with offsets use `y <- edgeR::catchSalmon(dirs)` (divides out overdispersion from Gibbs samples) or `DGEList(counts = txi$counts, ...)` plus `scaleOffset(y, log(txi$length))`.

### limma-voom (large n, repeated measures, many contrasts)

```r
suppressPackageStartupMessages({library(limma); library(edgeR)})
design <- model.matrix(~ 0 + condition + batch)            # cell-means for easy contrasts
colnames(design) <- sub("^condition", "", colnames(design))
y <- DGEList(cts); y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
v <- voomWithQualityWeights(y, design, plot = TRUE)        # down-weights low-quality libraries
# repeated measures: corfit <- duplicateCorrelation(v, design, block = coldata$subject); then block=, correlation=corfit$consensus in lmFit
fit <- lmFit(v, design)
cm <- makeContrasts(treated_vs_control = treated - control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE)
tt <- topTable(fit2, coef = "treated_vs_control", n = Inf)
write.csv(tt, "limma_voom_results.csv")
```

Validate: voom mean-variance trend monotone (hook at low counts = under-filtered); `plotSA(fit2)`.

### Batch handling rules

| Goal | Do | Never |
|---|---|---|
| DE with batch | `~ batch + condition` (first choice; all three frameworks) | pre-correct counts then test without batch (inflated FDR) |
| Counts for a tool that takes no covariates (some classifiers, deconvolution) | `sva::ComBat_seq(counts, batch = batch, group = condition)` | ComBat_seq when batch is confounded with condition |
| PCA/heatmap display | `limma::removeBatchEffect(vst_or_logcpm, batch = batch, design = model.matrix(~condition))` | feed this matrix to DESeq2/edgeR/limma - visualisation only |
| Unknown batch | `sva::svaseq` / `RUVSeq::RUVg` surrogate variables added to the design | more SVs than n/4 |
| Batch == condition (each condition sequenced separately) | report that the effect is unidentifiable; no statistical rescue | any "correction" |

## 6. Splicing

| Question | Tool | Input | Note |
|---|---|---|---|
| Event-level PSI (SE, A5SS, A3SS, MXE, RI), 2 groups | **rMATS-turbo** | STAR BAMs or FASTQ | `--libType` from inferred strandedness; `--readLength` must match |
| Exon-usage, complex designs | **DEXSeq** (`dexseq_prepare_annotation.py` → `dexseq_count.py` → `DEXSeqDataSetFromHTSeq`) | BAMs, flattened GTF | handles batch in `reducedModel` |
| Fast transcript-based PSI from Salmon/kallisto TPM | **SUPPA2** (`generateEvents` → `psiPerEvent` → `diffSplice`) | quant.sf | no alignment needed |
| Isoform-level DTU with uncertainty | `DRIMSeq` + `stageR`, or `fishpond::swish` with Gibbs samples | Salmon `--numGibbsSamples` | two-stage FDR control |
| All of the above reproducibly | **nf-core/rnasplice** 2.0.0 | `sample,fastq_1,fastq_2,strandedness,condition` | |

```bash
set -euo pipefail
GTF=gencode.v50.primary_assembly.annotation.gtf
READLEN=101
source strand_check/strand.env
ls star/ctrl*.Aligned.sortedByCoord.out.bam | paste -sd, > rmats/b1.txt
ls star/treat*.Aligned.sortedByCoord.out.bam | paste -sd, > rmats/b2.txt
rmats.py --b1 rmats/b1.txt --b2 rmats/b2.txt --gtf "${GTF}" -t paired --readLength "${READLEN}" --variable-read-length --libType "${RMATS_LIB}" --nthread 8 --od rmats/out --tmp rmats/tmp --novelSS
```

Validate: `awk -F'\t' 'NR>1 && $20<0.05' rmats/out/SE.MATS.JC.txt | wc -l` (column 20 = FDR in JC files); sashimi-plot (`rmats2sashimiplot`) 3 top events before believing any. Novel-junction STAR runs should use `--twopassMode Basic`.

## 7. Fusions

| Tool | Input | Strength |
|---|---|---|
| **Arriba** | STAR run with chimeric output (`--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreJunctionNonGTAG 0 --chimScoreDropMax 30 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMtype BAM Unsorted`) | fast, blacklist/known-fusion/protein-domain aware, clinical-grade filtering |
| **STAR-Fusion** | `Chimeric.out.junction` or FASTQ + CTAT genome lib | FusionInspector validation, coding-effect annotation |
| **nf-core/rnafusion** 4.2.0 | FASTQ (GRCh38 only) | Arriba + STAR-Fusion + FusionCatcher + fusion-report consensus |

```bash
set -euo pipefail
GENOME=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
ARRIBA_DB=arriba_db
SAMPLE=S1
arriba -x "star_chim/${SAMPLE}.Aligned.out.bam" -o "fusions/${SAMPLE}.fusions.tsv" -O "fusions/${SAMPLE}.fusions.discarded.tsv" -a "${GENOME}" -g "${GTF}" -b "${ARRIBA_DB}/blacklist_hg38_GRCh38_v2.5.0.tsv.gz" -k "${ARRIBA_DB}/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz" -t "${ARRIBA_DB}/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz" -p "${ARRIBA_DB}/protein_domains_hg38_GRCh38_v2.5.0.gff3"
```

Validate: `awk -F'\t' 'NR>1 && $15=="high"' fusions/${SAMPLE}.fusions.tsv | cut -f1,2,15 | head` (column 15 = confidence); every reported fusion has split reads **and** discordant mates, breakpoints at exon boundaries, and is absent from the matched-normal/blacklist. Database filenames follow the Arriba release you installed - list the directory rather than assuming the version string. Report fusions for research use only; clinical reporting needs orthogonal confirmation.

## 8. Pathway and gene-set analysis

Rules: use the **universe** of tested genes (not the genome) for ORA; rank by a signed statistic (`stat` from DESeq2, `t` from limma, `F`-signed logFC from edgeR) for GSEA, never by p-value alone; map IDs once (`ENSG` → ENTREZ/SYMBOL via `org.Hs.eg.db`), drop duplicates by max |stat|; report normalised enrichment and leading edge, not just p.

```r
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr); library(fgsea); library(enrichplot)})
res <- read.csv("deseq2_results.csv", row.names = 1)
res$ensembl <- sub("\\..*$", "", rownames(res))
map <- bitr(res$ensembl, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
res <- merge(res, map, by.x = "ensembl", by.y = "ENSEMBL")
universe <- unique(res$ENTREZID)
sig <- unique(res$ENTREZID[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1])
ego  <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego  <- simplify(ego)                                   # collapse redundant GO terms
ekeg <- enrichKEGG(sig, organism = "hsa", universe = universe)   # needs network access to KEGG
# GSEA with fgsea on MSigDB Hallmark (msigdbr >= 10 uses collection =; older versions category =)
hall <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hall$entrez_gene, hall$gs_name)
res_stat <- res[!is.na(res$stat), ]; res_stat <- res_stat[order(-abs(res_stat$stat)), ]; res_stat <- res_stat[!duplicated(res_stat$ENTREZID), ]
ranks <- sort(setNames(res_stat$stat, res_stat$ENTREZID), decreasing = TRUE)
set.seed(1); gsea <- fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500, eps = 0)
gsea <- gsea[order(gsea$padj), ]
write.csv(as.data.frame(gsea[, -"leadingEdge"]), "fgsea_hallmark.csv")
dotplot(ego, showCategory = 20); gseaplot2(GSEA(ranks, TERM2GENE = hall[, c("gs_name", "entrez_gene")]), geneSetID = 1)
```

Validate: `ego@universe` length ≈ number of tested genes (not ~20 000 if you tested 14 000); `gsea$padj` not all identical (indicates ties from a p-value ranking); the top Hallmark sets make biological sense for a positive control contrast. If `stat` is missing (limma: use `t`; edgeR: `sign(logFC) * -log10(PValue)`). `msigdbr` also serves mouse (`species = "Mus musculus"`) and other collections (`C2` with `subcollection = "CP:REACTOME"`, `C5` GO).

## 9. QC after alignment

RSeQC and Picard both need the strand setting from step 1.

```bash
set -euo pipefail
BAM=star/S1.Aligned.sortedByCoord.out.bam
BED12=gencode.v50.bed12
REFFLAT=gencode.v50.refFlat.txt
RRNA_INTERVALS=gencode.v50.rRNA.interval_list
source strand_check/strand.env
infer_experiment.py -r "${BED12}" -i "${BAM}" > qc/S1.infer_experiment.txt
read_distribution.py -r "${BED12}" -i "${BAM}" > qc/S1.read_distribution.txt
geneBody_coverage.py -r "${BED12}" -i "${BAM}" -o qc/S1
junction_annotation.py -r "${BED12}" -i "${BAM}" -o qc/S1 2> qc/S1.junction_annotation.txt
gatk CollectRnaSeqMetrics -I "${BAM}" -O qc/S1.rnaseq_metrics.txt --REF_FLAT "${REFFLAT}" --STRAND_SPECIFICITY "${PICARD_STRAND}" --RIBOSOMAL_INTERVALS "${RRNA_INTERVALS}"
gatk CollectInsertSizeMetrics -I "${BAM}" -O qc/S1.insert_size.txt -H qc/S1.insert_size.pdf
multiqc -o qc/multiqc qc star quant
```

Validate: `grep -A2 '^## METRICS' qc/S1.rnaseq_metrics.txt | tail -2 | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[i]=$i} NR==2{for(i=1;i<=NF;i++)print h[i]"\t"$i}' | grep -E 'PCT_(RIBOSOMAL|CODING|UTR|INTRONIC|INTERGENIC|MRNA|CORRECT_STRAND)_|5PRIME_TO_3PRIME'` → `PCT_RIBOSOMAL_BASES`, `PCT_MRNA_BASES`, `PCT_INTRONIC_BASES`, `MEDIAN_5PRIME_TO_3PRIME_BIAS` in the SKILL.md table range. A wrong `STRAND_SPECIFICITY` shows as `PCT_CORRECT_STRAND_READS` ≈ 0.5 or ≈ 0 - fix the flag, do not report the number.

refFlat from GTF: `gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf stdout | awk 'BEGIN{OFS="\t"}{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > refFlat.txt`. rRNA interval list: `samtools view -H "${BAM}" > rRNA.interval_list; grep -E 'gene_type "(rRNA|rRNA_pseudogene|Mt_rRNA)"' "${GTF}" | awk 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4,$5,$7,$10}' >> rRNA.interval_list`.

rRNA fraction without Picard: quantify against an rRNA FASTA only (`salmon quant -i rrna_idx -l A ...` → `percent_mapped`), or `sortmerna` (nf-core `--remove_ribo_rna`). Threshold in SKILL.md QC table.

## 10. Visualisation

```r
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(pheatmap); library(EnhancedVolcano)})
vsd <- vst(dds, blind = TRUE)                                           # blind for QC
p <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
ggplot(p, aes(PC1, PC2, colour = condition, shape = batch)) + geom_point(size = 3) + coord_fixed()
mat_vis <- limma::removeBatchEffect(assay(vsd), batch = vsd$batch, design = model.matrix(~ condition, colData(vsd)))  # display only
top <- head(order(res_unshrunk$padj), 50)
pheatmap(mat_vis[top, ], scale = "row", annotation_col = as.data.frame(colData(dds)[, c("condition", "batch")]), show_colnames = FALSE)
EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1, title = coef)   # shrunken LFC on x, padj on y
plotMA(res, ylim = c(-5, 5)); plotCounts(dds, gene = rownames(res)[which.min(res$padj)], intgroup = "condition")
```

Validate: sample distance heatmap (`dist(t(assay(vsd)))`) groups replicates; a PC1 that tracks library size or RIN means a covariate is missing from the design, not a reason to remove samples.
