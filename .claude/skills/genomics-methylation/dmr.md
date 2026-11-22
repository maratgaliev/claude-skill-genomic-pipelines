# Methylation: DMR callers, designs, arrays

Supporting file for `.claude/skills/genomics-methylation/SKILL.md`. Inputs assumed: per-sample Bismark coverage (`chr start end %meth countM countU`) or MethylDackel `--methylKit` output. Record `sessionInfo()` with every run.

## Caller choice

| Caller | Input | Model | Best for | Avoid when |
| --- | --- | --- | --- | --- |
| DSS | cov | Beta-binomial, shrinkage, optional smoothing | 1–3 reps/group WGBS; multi-factor via `DMLfit.multiFactor` | Need region-level FDR guarantee |
| dmrseq (bsseq) | BSseq | Region-level permutation, controls FDR on regions | ≥3 reps/group WGBS, honest FDR | <2 reps/group (permutations impossible) |
| methylKit | methylKit-format or cov | Logistic regression / Fisher, overdispersion correction, covariates | RRBS/targeted, many samples, batch covariates | Smoothing wanted for sparse WGBS |
| metilene | tab matrix | Segmentation + MWU/2D-KS, no R | Fast genome-wide WGBS, two groups | Covariates or paired designs |
| limma on M-values | β matrix | Linear model | Arrays (EPIC), per-probe | Sequencing counts with variable depth |

## methylKit (RRBS/targeted, covariates)

```r
suppressPackageStartupMessages(library(methylKit))
files <- list("meth/C1.cov.gz","meth/C2.cov.gz","meth/T1.cov.gz","meth/T2.cov.gz")
ids   <- list("C1","C2","T1","T2")
obj <- methRead(files, sample.id=ids, assembly="GRCh38", treatment=c(0,0,1,1),
                pipeline="bismarkCoverage", context="CpG", mincov=10)
obj <- filterByCoverage(obj, lo.count=10, hi.perc=99.9)
obj <- normalizeCoverage(obj, method="median")
tiles <- tileMethylCounts(obj, win.size=1000, step.size=1000, cov.bases=5)
meth  <- unite(tiles, min.per.group=2L)
covs  <- data.frame(batch=c("a","b","a","b"), age=c(40,55,43,60))
diff  <- calculateDiffMeth(meth, covariates=covs, overdispersion="MN", test="Chisq", mc.cores=8)
dmrs  <- getMethylDiff(diff, difference=10, qvalue=0.05)
write.table(getData(dmrs), "dmr/methylKit_tiles.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```
Checks: `getCorrelation(meth, plot=FALSE)` ≥0.85 within group; `PCASamples(meth)` PC1 separates condition not batch; `diff` object `nrow` before/after `getMethylDiff`.

## dmrseq (region-level FDR, ≥3 reps)

```r
suppressPackageStartupMessages({library(bsseq); library(dmrseq)})
bs <- read.bismark(files=c("meth/C1.cov.gz","meth/C2.cov.gz","meth/C3.cov.gz","meth/T1.cov.gz","meth/T2.cov.gz","meth/T3.cov.gz"),
                   colData=data.frame(cond=rep(c("ctrl","case"), each=3), row.names=c("C1","C2","C3","T1","T2","T3")),
                   rmZeroCov=TRUE, strandCollapse=TRUE)
keep <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0)==0)
bs   <- bs[keep,]
regions <- dmrseq(bs, testCovariate="cond", cutoff=0.05, minNumRegion=5, BPPARAM=BiocParallel::MulticoreParam(8))
sig <- regions[regions$qval < 0.05]
rtracklayer::export(sig, "dmr/dmrseq_sig.bed")
```
Checks: `length(regions)` vs `length(sig)`; `table(sign(sig$stat))`; `plotDMRs(bs, regions=sig[1:3], testCovariate="cond")`.

## metilene (two groups, no R)

```bash
set -euo pipefail
OUT=dmr
metilene_input.pl --in1 meth/C1.bedGraph,meth/C2.bedGraph --in2 meth/T1.bedGraph,meth/T2.bedGraph --h1 ctrl --h2 case --out "${OUT}/metilene_input.tsv"
metilene -a ctrl -b case -t 8 -M 300 -m 10 -d 0.1 "${OUT}/metilene_input.tsv" > "${OUT}/metilene_raw.tsv"
metilene_output.pl -q "${OUT}/metilene_raw.tsv" -o "${OUT}/metilene" -p 0.05 -c 10 -d 0.1 -a ctrl -b case
```
Validate: `wc -l ${OUT}/metilene_qval.0.05.out`; columns are `chr start end q-value mean_diff nCpG p(MWU) p(2D-KS) mean_g1 mean_g2`.

## Multi-factor DSS

```r
suppressPackageStartupMessages(library(DSS))
design <- data.frame(cond=factor(c("ctrl","ctrl","case","case")), sex=factor(c("F","M","F","M")))
fit <- DMLfit.multiFactor(bs, design=design, formula=~cond+sex, smoothing=TRUE)
res <- DMLtest.multiFactor(fit, coef="condcase")
dmrs <- callDMR(res, p.threshold=0.01, minlen=50, minCG=3, dis.merge=100)
```
Note: multi-factor mode has no `delta` argument; filter DMRs post hoc by mean β difference computed from `getMeth(bs, regions=dmrs, type="raw")`.

## Annotating DMRs

```bash
set -euo pipefail
GTF=ref/gencode.v50.annotation.gtf.gz
CGI=ref/cpgIslandExt.GRCh38.bed
DMR=dmr/DSS_dmrs.bed
bedtools intersect -wa -wb -a "${DMR}" -b "${CGI}" | cut -f1-3 | sort -u > dmr/dmr_in_cgi.bed
zcat "${GTF}" | awk '$3=="gene"' | bedtools closest -d -a <(sort -k1,1 -k2,2n "${DMR}") -b <(sort -k1,1 -k4,4n - | awk 'BEGIN{OFS="\t"}{match($0,/gene_name "[^"]+"/); print $1,$4-1,$5,substr($0,RSTART+11,RLENGTH-12)}') > dmr/dmr_nearest_gene.tsv
```
Use `annotatr` (R) for promoter/shore/shelf/enhancer categories with GRCh38/GRCm39 builds; GTF acquisition in `.claude/skills/genomics/reference.md`.

## EPIC v2 / EPIC v1 / 450k arrays (SeSAMe)

```r
suppressPackageStartupMessages({library(sesame); library(SummarizedExperiment)})
sesameDataCache()
idat_dir <- "idats"
betas <- openSesame(idat_dir, prep="QCDPB", BPPARAM=BiocParallel::MulticoreParam(8))
qcs   <- openSesame(idat_dir, prep="", func=sesameQC_calcStats, funs=c("detection","intensity"))
qc_df <- do.call(rbind, lapply(qcs, as.data.frame))
sex   <- vapply(searchIDATprefixes(idat_dir), function(p) inferSex(getBetas(prepSesame(readIDATpair(p), "QCDPB"))), "")
write.csv(betas, "array/betas_sesame.csv")
write.csv(qc_df, "array/qc_sesame.csv")
sessionInfo()
```
Notes:
- `prep="QCDPB"` = Q (mask quality), C (channel inference), D (dye bias NL), P (pOOBAH detection masking), B (noob background). This is the SeSAMe default recommendation for EPIC/EPIC v2.
- EPIC v2 has replicate probes (suffix `_BC11`, `_TC21`, …); collapse with `betasCollapseToPfx()` before comparing to EPIC v1/450k coordinates.
- Platform is inferred from the IDAT; if `inferPlatformFromProbeIDs` disagrees with the metadata, stop and ask.
- Differential analysis: convert to M-values (`BetaValueToMValue`) and use `limma` with sex/age/cell-type (e.g. `EpiDISH`/`FlowSorted.Blood.EPIC`) covariates; DMR level via `DMRcate` with EPIC v2 annotation (`EPICv2` array type).
- minfi still works for 450k/EPIC v1 (`read.metharray.exp`, `preprocessNoob`, `detectionP`); for EPIC v2 prefer SeSAMe or `minfi` with the `IlluminaHumanMethylationEPICv2manifest` package.
- Sample QC: drop samples with mean detection p>0.01 on >10 % probes, sex mismatch, or median intensity log2 <10.5 (`sesameQC_plotBar`).

## Reporting checklist

- Caller, version, model (smoothing on/off, covariates, formula), thresholds (Δβ, FDR, minCpG, minlen).
- Number of CpGs tested after depth/SNP filters per sample.
- Direction split of DMRs; top loci with coverage plots.
- Assembly string (GRCh38 / GRCm39; name legacy mm10 explicitly if a cohort used it).
