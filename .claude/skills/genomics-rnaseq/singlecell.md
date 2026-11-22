# Single-cell and spatial RNA-seq

Supporting file for `.claude/skills/genomics-rnaseq/SKILL.md`. Covers droplet scRNA/snRNA (10x and open alternatives), QC, integration, annotation, pseudobulk DE, trajectories, multimodal demux and spatial. scATAC / multiome chromatin: `.claude/skills/genomics-chipseq/atac.md`. Long-read single-cell: `.claude/skills/genomics-longread/SKILL.md`. Reference bundles: `.claude/skills/genomics/reference.md`.

Versions verified 2026-11-22: Cell Ranger 10.1.0; Space Ranger 4.1.0; 10x human GEX reference `refdata-gex-GRCh38-2024-A` (GRCm39 mouse `refdata-gex-GRCm39-2024-A`) - still the current public bundle per nf-core/scrnaseq docs, confirm at `<check current: https://www.10xgenomics.com/support/software/cell-ranger/downloads>`; CellBender 0.4.0 (identical model to 0.3.x, Python 3.10–3.12); Seurat 6.0.0 on CRAN (v5 layer/`IntegrateLayers` API unchanged; 5.6.0 is the last 5.x); Scanpy 1.13.0 (`scanpy.external` deprecated → native `sc.pp.harmony_integrate`, `sc.pp.scrublet`, `sc.pp.bbknn`); kallisto 0.52.0 / kb-python; nf-core/scrnaseq 4.2.0; nf-core/spatialvi has no tagged release (dev only).

## 1. Inspect first (single-cell specifics)

| Check | Command | Decision |
|---|---|---|
| Chemistry | R1 read length: 26 → 3' v2, 28 → 3' v3/v3.1/v4 or 5'; `cellranger count` auto-detects (`--chemistry auto`) | whitelist / `-x` string / `--soloUMIlen` |
| R1/R2 orientation | `peek.sh *_R1_001.fastq.gz`: R1 should be barcode+UMI (short), R2 cDNA | swapped files break every tool silently |
| Nuclei vs cells | metadata; `pct_counts_mt` ≈ 0 and high intronic fraction → nuclei | intron-inclusive counting (default in Cell Ranger ≥7, `GeneFull` in STARsolo, `nac` in kb), `percent.mt` threshold ≈ 5 % |
| Existing outputs | `ls outs/`: `raw_feature_bc_matrix.h5` (needed for CellBender), `filtered_feature_bc_matrix`, `molecule_info.h5` (for `cellranger aggr`), `metrics_summary.csv` | skip re-counting if `cellranger --version` in `_versions` matches the plan |
| Multiplexing | `feature_reference.csv`, `Antibody Capture` / `Multiplexing Capture` rows in `features.tsv.gz` | HTO/CMO demux step required |
| Expected cells | `metrics_summary.csv` "Estimated Number of Cells" vs loading plan | >2× expected → ambient/empty calls; <0.5× → wetting failure or wrong chemistry |

## 2. Counting

### Cell Ranger (vendor pipeline; required for Loupe/10x support cases)

```bash
set -euo pipefail
REF=/refs/refdata-gex-GRCh38-2024-A
SAMPLE=S1
FASTQS=/data/fastq/S1
cellranger count --id="${SAMPLE}" --transcriptome="${REF}" --fastqs="${FASTQS}" --sample="${SAMPLE}" --create-bam=false --localcores=16 --localmem=64
test -s "${SAMPLE}/outs/filtered_feature_bc_matrix.h5"
```

Validate: `cut -d, -f1-6 "${SAMPLE}/outs/metrics_summary.csv"`; "Fraction Reads in Cells" ≥0.7, "Reads Mapped Confidently to Transcriptome" ≥0.3 (nuclei: lower, but "to Genome" ≥0.8), "Sequencing Saturation" ≥0.5, "Valid Barcodes" ≥0.9. `--create-bam` is mandatory since v8 (set `true` only if you need velocyto/souporcell). Introns are counted by default since v7 (`--include-introns false` to revert). Multi-library / multiplexed (CMO, hashing via Feature Barcode, Flex probes, VDJ): `cellranger multi --csv config.csv`. Multi-sample merge: `cellranger aggr --csv aggr.csv --normalize=none` (normalise downstream, not by read subsampling).

### STARsolo (open, same output layout, ~2–3× faster)

```bash
set -euo pipefail
STAR_IDX=star_idx_gencode50_oh100
WL=3M-february-2018.txt
SAMPLE=S1
R1=$(ls fastq/${SAMPLE}_S*_R1_001.fastq.gz | paste -sd,)
R2=$(ls fastq/${SAMPLE}_S*_R2_001.fastq.gz | paste -sd,)
STAR --runThreadN 16 --genomeDir "${STAR_IDX}" --readFilesIn "${R2}" "${R1}" --readFilesCommand zcat --outFileNamePrefix "solo/${SAMPLE}." --soloType CB_UMI_Simple --soloCBwhitelist "${WL}" --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloFeatures Gene GeneFull Velocyto --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --outSAMtype None
```

Validate: `cat solo/${SAMPLE}.Solo.out/GeneFull/Summary.csv`. Read order is **cDNA then barcode** (`R2 R1`). Whitelists: v2 `737K-august-2016.txt` (`--soloUMIlen 10`), v3/v3.1 `3M-february-2018.txt`, v4 `3M-3pgex-may-2023.txt`; 5' kits have their own files - always copy the whitelist from the installed Cell Ranger (`lib/python/cellranger/barcodes/`) rather than from memory. `Velocyto` output feeds scVelo. Use `GeneFull` for nuclei, `Gene` for cells.

### kb-python (kallisto|bustools; lowest resources; `nac` workflow gives nascent/mature)

```bash
set -euo pipefail
GENOME=GRCh38.primary_assembly.genome.fa
GTF=gencode.v50.primary_assembly.annotation.gtf
SAMPLE=S1
kb ref -i kb/index.idx -g kb/t2g.txt -f1 kb/cdna.fa -f2 kb/nascent.fa -c1 kb/cdna_t2c.txt -c2 kb/nascent_t2c.txt --workflow nac "${GENOME}" "${GTF}"
kb count -i kb/index.idx -g kb/t2g.txt -c1 kb/cdna_t2c.txt -c2 kb/nascent_t2c.txt --workflow nac -x 10XV3 -t 8 --h5ad --sum total -o "kb/${SAMPLE}" fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz
```

Validate: `jq '{n_processed, n_pseudoaligned, p_pseudoaligned}' kb/${SAMPLE}/run_info.json`; `kb --list` shows technology strings. Standard (cells) workflow: drop `--workflow nac`, `-f2/-c1/-c2`. Alternative with identical philosophy: `simpleaf` (piscem + alevin-fry; nf-core/scrnaseq default). `salmon alevin` no longer exists in Salmon 2.x.

## 3. Ambient RNA

| Tool | Input | When |
|---|---|---|
| **CellBender** `remove-background` | **raw** (unfiltered) matrix `.h5` | default for droplet data; GPU strongly recommended; learns ambient profile + calls cells |
| **SoupX** | raw + filtered + clusters | CPU-only, fast; needs clustering first; good when CellBender unavailable |
| `DropletUtils::emptyDrops` | raw matrix | cell calling only (no ambient removal) |

```bash
set -euo pipefail
SAMPLE=S1
cellbender remove-background --cuda --input "${SAMPLE}/outs/raw_feature_bc_matrix.h5" --output "cellbender/${SAMPLE}.h5" --fpr 0.01 --epochs 150
ptrepack --complevel 5 "cellbender/${SAMPLE}_filtered.h5:/matrix" "cellbender/${SAMPLE}_filtered_seurat.h5:/matrix"
```

Validate: open `cellbender/${SAMPLE}_report.html`; training loss converged (no late upswing); fraction of counts removed per cell typically 1–10 % (>30 % → overloaded/failed run, lower `--fpr` or inspect). `--expected-cells/--total-droplets-included` are optional since 0.3 but set them when the barcode-rank plot has no clear knee. `ptrepack` (PyTables) makes the output readable by `Seurat::Read10X_h5`; Scanpy reads it directly (`sc.io.read_10x_h5` / `cellbender.remove_background.downstream.anndata_from_h5`). Remember CellBender outputs are *denoised* counts: still integers, usable for all downstream steps including pseudobulk.

SoupX (R):

```r
suppressPackageStartupMessages({library(SoupX); library(Seurat)})
sc <- load10X("S1/outs")                 # needs raw + filtered + analysis/clustering from Cell Ranger
sc <- autoEstCont(sc)                    # rho = contamination fraction; typical 0.02-0.10
adj <- adjustCounts(sc, roundToInt = TRUE)
```

## 4. QC: thresholds are tissue-dependent, use MAD outliers

Rules: compute per sample; flag with median ± 3 MAD on log scale (`sum`, `detected`) and upper 3 MAD (or fixed cap) for mito; **look** at the distributions before dropping; never apply PBMC cut-offs (`percent.mt < 5`) to heart, kidney, liver or tumours.

| Tissue / prep | `percent.mt` default upper | `nFeature` lower | Notes |
|---|---|---|---|
| PBMC / immune | 10 % | 200 | neutrophils have low nFeature (~100–300) - do not filter them away if present |
| Solid tissue (lung, gut, skin, tumour) | 15–20 % | 300 | stressed epithelia carry 15–25 % mito legitimately |
| Heart, kidney, liver, muscle | 25–40 % (cardiomyocytes >50 %) | 300 | use MAD within cell type after a first clustering |
| Nuclei (any tissue) | 5 % (ideally <1 %) | 200 | mito reads in nuclei = cytoplasmic debris; also filter `percent.ribo` |
| Cell lines | 15 % | 500 | |
| Upper `nCount` / doublet guard | >5 MAD above median, or leave to doublet caller | | hard caps bias against large cell types |

```r
suppressPackageStartupMessages({library(scuttle); library(scDblFinder); library(SingleCellExperiment)})
sce <- DropletUtils::read10xCounts(dirs, sample.names = ids, col.names = TRUE)   # filtered or CellBender matrices
rownames(sce) <- rowData(sce)$ID
is_mito <- grepl("^MT-", rowData(sce)$Symbol)                       # mouse: "^mt-"
sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is_mito))
qc <- perCellQCFilters(colData(sce), sub.fields = "subsets_Mito_percent", nmads = 3, batch = sce$Sample)
qc$high_mito_fixed <- sce$subsets_Mito_percent > 20                 # tissue-dependent cap, see table
table(sce$Sample, qc$discard | qc$high_mito_fixed)
sce <- sce[, !(qc$discard | qc$high_mito_fixed)]
sce <- scDblFinder(sce, samples = "Sample", clusters = TRUE)         # run per sample, BEFORE integration
table(sce$Sample, sce$scDblFinder.class)
sce <- sce[, sce$scDblFinder.class == "singlet"]
```

Validate: doublet rate per sample ≈ 0.8 % per 1 000 cells recovered (10x), so 10 k cells → ~8 %; a 25 % call rate means overloading or a cluster of real intermediate cells - check `scDblFinder.score` distribution. Scanpy equivalent: `sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True)` + `sc.pp.scrublet(adata, batch_key='sample')` (expected doublet rate 0.06 default; set `expected_doublet_rate` from loading).

## 5. Normalisation, integration, clustering

### Seurat (v5 API; Seurat 6.0.0 keeps it)

```r
suppressPackageStartupMessages({library(Seurat); library(SeuratWrappers)})
options(Seurat.object.assay.version = "v5")
mats <- lapply(setNames(dirs, ids), Read10X_h5)                      # or CellBender *_filtered_seurat.h5
obj <- merge(CreateSeuratObject(mats[[1]], project = ids[1], min.cells = 3), y = mapply(CreateSeuratObject, mats[-1], project = ids[-1]), add.cell.ids = ids)
obj$sample <- obj$orig.ident
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
# v5 keeps one counts layer per sample (counts.S1, counts.S2 ...); pipelines run per layer
obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 50)
# option: SCTransform v2 instead of the three steps above (vst.flavor = "v2" is the default in v5)
# obj <- SCTransform(obj, vars.to.regress = "percent.mt", vst.flavor = "v2"); RunPCA(obj)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
# alternatives: CCAIntegration / RPCAIntegration (new.reduction = "integrated.cca"), scVIIntegration (needs reticulate + scvi-tools; normalization.method = "SCT" when SCT was used)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) |> FindClusters(resolution = 0.5, algorithm = 4) |> RunUMAP(reduction = "harmony", dims = 1:30)
obj <- JoinLayers(obj)                                               # collapse layers before FindMarkers / export
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
```

Validate: `Layers(obj)` shows a single `counts`/`data` after `JoinLayers`; per-cluster sample composition (`prop.table(table(obj$seurat_clusters, obj$sample), 1)`) not dominated by one sample unless biology says so; `algorithm = 4` (Leiden) requires `leidenalg` via reticulate - fall back to default Louvain (1) if absent. `ElbowPlot(obj, ndims = 50)` to justify `dims`.

### Scanpy / scvi-tools

```python
import scanpy as sc, anndata as ad
adatas = {s: sc.io.read_10x_h5(f"cellbender/{s}_filtered.h5") for s in samples}   # scanpy <1.13: sc.read_10x_h5
for a in adatas.values(): a.var_names_make_unique()
adata = ad.concat(adatas, label="sample", index_unique="-")
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts", batch_key="sample")   # seurat_v3 needs raw counts
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")              # scanpy 1.13: mask=
sc.pp.harmony_integrate(adata, batch="sample")                         # scanpy >=1.13 (1.12: sc.external.pp.harmony_integrate(adata, key="sample")) -> obsm["X_pca_harmony"]
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False, key_added="leiden_0.5")
sc.tl.rank_genes_groups(adata, "leiden_0.5", method="wilcoxon", layer=None)
adata.write_h5ad("integrated.h5ad")
```

scVI when batches are strong or many (>10 samples, multiple chemistries):

```python
import scvi
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample", categorical_covariate_keys=["chemistry"], continuous_covariate_keys=["pct_counts_mt"])
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb"); model.train()
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI"); sc.tl.umap(adata); sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
# scANVI (label transfer) builds on the trained SCVI model: scvi.model.SCANVI.from_scvi_model(model, labels_key=..., unlabeled_category="Unknown")
```

Validate integration (any method): `scib-metrics` (`Benchmarker`) or at minimum iLISI/kBET-style sample mixing per cluster **and** preservation of known markers/cell types; over-integration erases condition-specific states - compare clusters with and without integration before trusting a disappearance.

## 6. Annotation

| Tool | Mode | Reference | Use |
|---|---|---|---|
| **SingleR** + `celldex` | per-cell correlation | `celldex::fetchReference("hpca", "2024-02-26")`, `"blueprint_encode"`, `"monaco_immune"`, `"mouse_rnaseq"` | fast first pass on any tissue; `labels = ref$label.main` then `label.fine` |
| **CellTypist** | logistic-regression models | `celltypist.models.download_models(model=["Immune_All_Low.pkl", ...])` | immune and many tissue models; `majority_voting=True` with your clusters |
| **Azimuth** | reference mapping (Seurat) | `pbmcref`, `lungref`, `kidneyref`, `bonemarrowref`, ... (`SeuratData`) | PBMC/lung with full hierarchy and QC mapping scores |
| Marker tables / `scType`, `decoupler` + PanglaoDB | cluster-level scoring | curated markers | sanity check; always report the markers used |
| LLM-based (GPTCelltype etc.) | cluster-level | none | hypothesis only, never final |

```r
suppressPackageStartupMessages({library(SingleR); library(celldex)})
ref <- celldex::fetchReference("hpca", "2024-02-26")
pred <- SingleR(test = GetAssayData(obj, layer = "data"), ref = ref, labels = ref$label.main, clusters = obj$seurat_clusters)
obj$singler_cluster <- pred$labels[match(obj$seurat_clusters, rownames(pred))]
plotScoreHeatmap(pred)
```

```python
import celltypist
celltypist.models.download_models(model=["Immune_All_Low.pkl"])
# CellTypist expects log1p(normalize_total target_sum=1e4) in .X
pred = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True, over_clustering="leiden_0.5")
adata = pred.to_adata()
```

Validate: ≥80 % of cells in a cluster share the label; `delta.next` (SingleR) / `conf_score` (CellTypist) distributions not bimodal; canonical markers (`PTPRC`, `CD3E`, `MS4A1`, `LYZ`, `EPCAM`, `PECAM1`, `COL1A1`, `MKI67`) agree with labels. Disagreement between two tools = unlabeled, not majority vote.

## 7. Differential expression between conditions: pseudobulk by default

Cells from one sample are not replicates. For any multi-sample design (condition, treatment, genotype) aggregate counts per sample × cell type and test with edgeR/DESeq2; cell-level Wilcoxon/`FindMarkers` across conditions inflates false positives by orders of magnitude and is acceptable only for marker discovery within a sample set, or when n = 1 per group (then report it as exploratory).

```r
suppressPackageStartupMessages({library(scuttle); library(edgeR); library(SingleCellExperiment)})
sce <- as.SingleCellExperiment(JoinLayers(obj))                     # counts assay must be raw/denoised integers
pb  <- aggregateAcrossCells(sce, ids = colData(sce)[, c("sample", "celltype")], statistics = "sum")
pb  <- pb[, pb$ncells >= 20]                                         # drop sample x celltype combos with too few cells
res_list <- list()
for (ct in unique(pb$celltype)) {
  sub <- pb[, pb$celltype == ct]
  if (length(unique(sub$condition)) < 2 || any(table(sub$condition) < 2)) next
  cond  <- relevel(factor(sub$condition), ref = "control")
  des   <- model.matrix(~ cond)                                      # add + batch / + sex when present
  y <- DGEList(counts(sub), group = cond); y <- y[filterByExpr(y, design = des), , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
  fit <- glmQLFit(y, des, robust = TRUE); q <- glmQLFTest(fit, coef = 2)
  res_list[[ct]] <- cbind(celltype = ct, topTags(q, n = Inf)$table)
}
write.csv(do.call(rbind, res_list), "pseudobulk_edgeR.csv", row.names = TRUE)
# one-call alternatives: scran::pseudoBulkDGE(pb, label = pb$celltype, design = ~ condition, coef = "conditiontreated", condition = pb$condition)
# DESeq2: DESeqDataSetFromMatrix(counts(sub), colData(sub), ~ condition) per cell type; Python: decoupler.pp.pseudobulk(adata, sample_col="sample", groups_col="celltype", layer="counts") -> PyDESeq2
```

Validate: pseudobulk PCA per cell type separates condition rather than sample index; genes DE in *every* cell type with the same sign → ambient/contamination or batch, inspect; number of DE genes scales with number of cells aggregated (power), not with cell-type "importance". Compositional change (proportion shifts) is a separate test: `speckle::propeller`, `scCODA`, or `sccomp`.

## 8. Trajectory, velocity, communication (pointers)

| Question | Tool | Input requirement |
|---|---|---|
| Pseudotime along a continuous process | **slingshot** (lineages on UMAP/PCA from clusters), **Monocle3** (`learn_graph` → `order_cells`, root chosen biologically) | a genuinely continuous population; subset to the lineage first |
| Direction from splicing kinetics | **scVelo** (`dynamical` mode) → **CellRank 2** for fate probabilities | spliced/unspliced layers from STARsolo `Velocyto`, kb `nac`, or `velocyto run10x` (needs BAM) |
| Ligand–receptor | **LIANA+** (consensus), **CellChat v2** | annotated cell types, ideally per condition |
| Gene regulatory | **pySCENIC/SCENIC+** | counts + (for +) ATAC |

Trajectory output is a hypothesis: validate root choice with known early markers, and report that velocity on nuclei or 3' v2 data is unreliable (low unspliced fraction).

## 9. Multimodal: CITE-seq, hashing, genetic demux

| Task | Tool | Note |
|---|---|---|
| HTO / CMO demultiplexing | `Seurat::HTODemux` (CLR-normalised HTO assay), `DropletUtils::hashedDrops`, `cellranger multi` (CMO) | compare two callers; discard "Negative"; "Doublet" calls also train your doublet threshold |
| ADT (CITE-seq) normalisation | `NormalizeData(normalization.method = "CLR", margin = 2)` or `dsb` (uses empty droplets) | WNN (`FindMultiModalNeighbors`) for joint clustering |
| Genotype demux (pooled donors) | `vireo` (with/without genotypes), `souporcell`, `demuxlet/freemuxlet` | needs BAM (`--create-bam true`) |
| Hashing + genotype | `hadge` pipeline (nf-core-adjacent) | nf-core/scrnaseq does not demultiplex |

```r
obj[["HTO"]] <- CreateAssay5Object(counts = hto_counts[, colnames(obj)])
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
table(obj$HTO_classification.global)                                 # Singlet / Doublet / Negative
obj <- subset(obj, HTO_classification.global == "Singlet")
```

Validate: singlet fraction ≥70 %; `RidgePlot(obj, assay = "HTO", features = rownames(obj[["HTO"]]))` bimodal per hashtag; singlets per hashtag roughly match pooling ratios.

## 10. Spatial transcriptomics

| Platform | Primary pipeline | Unit | Downstream |
|---|---|---|---|
| Visium (v1/v2, CytAssist, FFPE probe-based) | `spaceranger count` (4.1.0) | 55 µm spots (multi-cell) | Seurat `Load10X_Spatial`, Squidpy, deconvolution (cell2location, RCTD/`spacexr`) |
| **Visium HD** | `spaceranger count` with HD slide; outputs 2/8/16 µm bins (`binned_outputs/square_008um`) + cell segmentation (4.0+) | bins or segmented cells | Seurat ≥5.4 `Load10X_Spatial(bin.size = 8)`, Squidpy/`spatialdata`; `--hd_bin_size` in nf-core/spatialvi |
| **Xenium** (in situ, targeted panels; 5K panel) | on-instrument + Xenium Ranger (`<check current: https://www.10xgenomics.com/support/software/xenium-ranger>`) | single cells (segmentation) | Seurat `LoadXenium`, `spatialdata-io` → Squidpy; no deconvolution needed |
| MERSCOPE / CosMx | vendor outputs | single cells | `spatialdata-io`, Seurat `ReadVizgen`/`LoadNanostring` |
| Slide-seq / Stereo-seq / Curio | vendor outputs | beads/bins | `spatialdata-io`, Seurat `ReadSlideSeq`, `LoadCurioSeeker` |

```bash
set -euo pipefail
REF=/refs/refdata-gex-GRCh38-2024-A
PROBESET=/opt/spaceranger/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv
SAMPLE=V1
spaceranger count --id="${SAMPLE}" --transcriptome="${REF}" --probe-set="${PROBESET}" --fastqs=/data/fastq/V1 --sample="${SAMPLE}" --cytaimage=/data/img/V1_cytassist.tif --image=/data/img/V1_HE.tif --slide=V53A06-123 --area=A1 --create-bam=false --localcores=16 --localmem=64
test -s "${SAMPLE}/outs/spatial/tissue_positions.csv"
```

Validate: `web_summary.html` tissue detection overlays tissue correctly (manual alignment via Loupe when it does not); "Fraction Reads in Spots Under Tissue" ≥0.5; probe-based runs report "Reads Mapped Confidently to Probe Set". Probe set file name/path follows the installed Space Ranger - `ls $(dirname $(which spaceranger))/../probe_sets/`. Fresh-frozen poly-A: drop `--probe-set` and `--cytaimage`. Visium HD: add `--slide`/`--area` from the HD slide ID; bins appear under `outs/binned_outputs/`.

```python
import squidpy as sq, scanpy as sc
adata = sq.read.visium("V1/outs")                                    # Visium HD: spatialdata_io.visium_hd("V1/outs", bin_size=8)
adata.var_names_make_unique(); adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_counts=500); sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata); sc.pp.log1p(adata); sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.pca(adata); sc.pp.neighbors(adata); sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2, directed=False)
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)       # Visium hex grid; "generic" for Xenium/HD cells
sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, n_jobs=8)  # spatially variable genes
sq.gr.nhood_enrichment(adata, cluster_key="leiden"); sq.pl.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.spatial_scatter(adata, color=["leiden", "total_counts"])
```

Seurat: `obj <- Load10X_Spatial("V1/outs"); obj <- SCTransform(obj, assay = "Spatial") |> RunPCA() |> FindNeighbors(dims = 1:30) |> FindClusters() |> RunUMAP(dims = 1:30); SpatialDimPlot(obj); FindSpatiallyVariableFeatures(obj, selection.method = "moransi")`. Spot deconvolution with a matched scRNA reference: `spacexr::run.RCTD` (R) or `cell2location` (Python, GPU) - required before any cell-type-specific claim on 55 µm spots; HD 8 µm bins and Xenium cells can be clustered/annotated directly (CellTypist/Azimuth/SingleR as in section 6). nf-core/spatialvi runs Space Ranger + Scanpy/Squidpy QC/clustering/SVG for Visium (samplesheet in SKILL.md).

Validate: spatial clusters follow histology (overlay on H&E); Moran's I top genes include known structural markers of the tissue; deconvolution proportions per spot sum to 1 and cell types land where histology expects them.

## 11. Pitfalls (single-cell specific)

| Symptom | Cause | Fix |
|---|---|---|
| `cellranger count` error "--create-bam" required | Cell Ranger ≥8 | add `--create-bam=false|true` |
| Chemistry auto-detect fails | R1 trimmed, R1/R2 swapped, or 5' v3/Flex | check R1 length; set `--chemistry` explicitly |
| Seurat `FindMarkers` error about multiple layers | v5 object not joined | `JoinLayers(obj)` after integration |
| `IntegrateLayers` error "fewer than two layers" | object already joined or single sample | `split(obj[["RNA"]], f = obj$sample)` first; a single sample needs no integration |
| `highly_variable_genes(flavor='seurat_v3')` warns about non-integers | ran on log data | pass `layer="counts"` |
| Leiden very slow / different results between runs | `leidenalg` default flavor | `flavor="igraph", n_iterations=2, directed=False`; fix `random_state` |
| Scanpy ≥1.13 deprecation warnings for `sc.external`, `use_rep`, `layer`, `copy=` | API migration to `use=`/`out=` accessors | follow the warning's named replacement; `pip install 'scanpy<1.13'` only to reproduce old results |
| CellBender h5 unreadable in Seurat | 0.3+ h5 layout | `ptrepack ... :/matrix` as above |
| Thousands of DE genes between conditions from `FindMarkers` | pseudoreplication | pseudobulk (section 7) |
| Mito-high cluster with few genes | dying cells or true high-mito type | compare with tissue table; subcluster; check `nFeature` jointly |
| Cluster of cells expressing every marker | doublets or ambient | doublet scores + `scDblFinder` before clustering; CellBender |
| Sample-specific clusters after Harmony | residual batch or real biology | scVI with `categorical_covariate_keys`; verify with marker genes before concluding |
| Visium spots over empty areas | tissue detection failure | manual alignment JSON (`--loupe-alignment`) and rerun |
