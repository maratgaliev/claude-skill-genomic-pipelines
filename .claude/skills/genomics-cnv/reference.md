# genomics-cnv reference: full command chains

Companion to `SKILL.md` (decision logic lives there). All blocks assume `set -euo pipefail`, GRCh38 with `chr` prefixes, duplicate-marked coordinate-sorted BAM/CRAM with index, and `ref_sanity.sh` already run. Replace `${UPPER_CASE}` placeholders at the top of each block.

Sections: 1 GATK gCNV (cohort, postprocess, case) · 2 GATK somatic CNV · 3 CNVkit step-by-step · 4 SV callers (smoove population, DELLY, Tiddit, dysgu, duphold, GRIDSS note, Jasmine) · 5 ExomeDepth R · 6 Resource files

## 1. GATK germline CNV (gCNV)

Interval preparation (`PreprocessIntervals`, `AnnotateIntervals`) is in `SKILL.md` step 1. Verified against the Broad `cnv_germline_cohort_workflow.wdl` defaults (padding 250, bin 1000, `OVERLAPPING_ONLY` everywhere, FilterIntervals GC 0.1–0.9, mappability ≥0.9, segdup ≤0.5, low-count 10 reads in 50% of samples, extreme-count 1–99th percentile in 90% of samples).

### 1a. Cohort mode, sharded (WGS or large WES cohorts)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
COHORT=cohort01
PRIORS=ploidy_priors.tsv
SHARD_SIZE=5000
THREADS=8
export MKL_NUM_THREADS="${THREADS}" OMP_NUM_THREADS="${THREADS}"
INPUTS=$(printf -- '-I %s ' counts/*.counts.hdf5)
gatk DetermineGermlineContigPloidy -L intervals.filtered.interval_list ${INPUTS} \
  --contig-ploidy-priors "${PRIORS}" --interval-merging-rule OVERLAPPING_ONLY \
  -O ploidy --output-prefix "${COHORT}"
gatk IntervalListTools -I intervals.filtered.interval_list --SUBDIVISION_MODE INTERVAL_COUNT \
  --SCATTER_CONTENT "${SHARD_SIZE}" -O shards
i=0
for SHARD in shards/temp_*/scattered.interval_list; do
  gatk GermlineCNVCaller --run-mode COHORT -L "${SHARD}" ${INPUTS} \
    --contig-ploidy-calls "ploidy/${COHORT}-calls" --annotated-intervals intervals.annotated.tsv \
    --interval-merging-rule OVERLAPPING_ONLY -O "gcnv/shard_${i}" --output-prefix "${COHORT}"
  i=$((i+1))
done
```

Validate: every `gcnv/shard_*/${COHORT}-calls/SAMPLE_*/` has `baseline_copy_number_t.tsv` and `log_q_c_tc.tsv`; `gcnv/shard_*/${COHORT}-tracking/` ELBO converged (last lines of `main_elbo_history.tsv` flat). Memory: ~1 GB per 1,000 intervals × samples/50; wall time hours for 100 WGS samples per shard. Run shards in parallel on a scheduler rather than this serial loop.

### 1b. Postprocess (one call per sample, all shards at once)

```bash
set -euo pipefail
COHORT=cohort01
SAMPLE_INDEX=0
MODELS=$(printf -- '--model-shard-path %s ' gcnv/shard_*/"${COHORT}"-model)
CALLS=$(printf -- '--calls-shard-path %s ' gcnv/shard_*/"${COHORT}"-calls)
SAMPLE=$(cat "gcnv/shard_0/${COHORT}-calls/SAMPLE_${SAMPLE_INDEX}/sample_name.txt")
gatk PostprocessGermlineCNVCalls ${MODELS} ${CALLS} \
  --contig-ploidy-calls "ploidy/${COHORT}-calls" --sample-index "${SAMPLE_INDEX}" \
  --autosomal-ref-copy-number 2 --allosomal-contig chrX --allosomal-contig chrY \
  --output-genotyped-intervals "calls/${SAMPLE}.intervals.vcf.gz" \
  --output-genotyped-segments "calls/${SAMPLE}.segments.vcf.gz" \
  --output-denoised-copy-ratios "calls/${SAMPLE}.denoisedCR.tsv"
bcftools view -i 'GT!="0" && QS>=50' "calls/${SAMPLE}.segments.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%CN\t%QS\t%NP]\n' > "calls/${SAMPLE}.cnv.tsv"
```

Sample index = position in the `-I` list given to `GermlineCNVCaller` (0-based); `sample_name.txt` guards against misassignment. Segment VCF FORMAT: `CN` copy number, `NP` points, `QA`/`QS`/`QSS`/`QSE` Phred qualities. Defaults to justify: `QS ≥ 50` for WES exon-level, `QS ≥ 100` WGS; require `NP ≥ 2` for exon deletions unless the gene has a single targeted exon. Validate: `bcftools view -H calls/${SAMPLE}.segments.vcf.gz | wc -l` within cohort range; no segments on contigs absent from `intervals.filtered.interval_list`.

### 1c. Case mode (new sample against an existing cohort model)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
COHORT=cohort01
S=case01
BAM=bams/case01.bam
gatk CollectReadCounts -I "${BAM}" -R "${REF}" -L intervals.preprocessed.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY -O "counts_case/${S}.counts.hdf5"
gatk DetermineGermlineContigPloidy --model "ploidy/${COHORT}-model" -I "counts_case/${S}.counts.hdf5" \
  -O "ploidy_case/${S}" --output-prefix "${S}"
for SHARD in gcnv/shard_*; do
  gatk GermlineCNVCaller --run-mode CASE --model "${SHARD}/${COHORT}-model" -I "counts_case/${S}.counts.hdf5" \
    --contig-ploidy-calls "ploidy_case/${S}/${S}-calls" -O "gcnv_case/${S}/$(basename "${SHARD}")" --output-prefix "${S}"
done
MODELS=$(printf -- '--model-shard-path %s ' gcnv/shard_*/"${COHORT}"-model)
CALLS=$(printf -- '--calls-shard-path %s ' gcnv_case/"${S}"/shard_*/"${S}"-calls)
gatk PostprocessGermlineCNVCalls ${MODELS} ${CALLS} --contig-ploidy-calls "ploidy_case/${S}/${S}-calls" \
  --sample-index 0 --autosomal-ref-copy-number 2 --allosomal-contig chrX --allosomal-contig chrY \
  --output-genotyped-intervals "calls/${S}.intervals.vcf.gz" --output-genotyped-segments "calls/${S}.segments.vcf.gz" \
  --output-denoised-copy-ratios "calls/${S}.denoisedCR.tsv"
```

Case mode uses the *preprocessed* intervals for counting but the model's own filtered interval list internally; the counts file must contain every interval the model was trained on (same `PreprocessIntervals` output, same kit). Validate: `ploidy_case/${S}/${S}-calls/SAMPLE_0/contig_ploidy.tsv` sex matches metadata. nf-core/raredisease consumes a cohort model built this way (`gcnvcaller_model` TSV; `nf-core/createpanelrefs` builds it).

## 2. GATK somatic CNV (PoN → denoise → ModelSegments → plots)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
DICT=ref/GRCh38.dict
SNPS=resources/common_biallelic_snps.GRCh38.interval_list
PON=cnv_pon.hdf5
TUMOR=tumor.bam
NORMAL=normal.bam
S=tumor
for BAM in normals/*.bam; do
  N=$(basename "${BAM}" .bam)
  gatk CollectReadCounts -I "${BAM}" -R "${REF}" -L intervals.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O "pon_counts/${N}.counts.hdf5"
done
PON_INPUTS=$(printf -- '-I %s ' pon_counts/*.counts.hdf5)
gatk --java-options "-Xmx16g" CreateReadCountPanelOfNormals ${PON_INPUTS} \
  --annotated-intervals intervals.annotated.tsv --minimum-interval-median-percentile 5.0 -O "${PON}"
gatk CollectReadCounts -I "${TUMOR}" -R "${REF}" -L intervals.preprocessed.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY -O "${S}.counts.hdf5"
gatk DenoiseReadCounts -I "${S}.counts.hdf5" --count-panel-of-normals "${PON}" \
  --standardized-copy-ratios "${S}.standardizedCR.tsv" --denoised-copy-ratios "${S}.denoisedCR.tsv"
gatk CollectAllelicCounts -I "${TUMOR}" -R "${REF}" -L "${SNPS}" -O "${S}.allelicCounts.tsv"
gatk CollectAllelicCounts -I "${NORMAL}" -R "${REF}" -L "${SNPS}" -O normal.allelicCounts.tsv
gatk ModelSegments --denoised-copy-ratios "${S}.denoisedCR.tsv" --allelic-counts "${S}.allelicCounts.tsv" \
  --normal-allelic-counts normal.allelicCounts.tsv --output-prefix "${S}" -O segments
gatk CallCopyRatioSegments -I "segments/${S}.cr.seg" -O "segments/${S}.called.seg"
gatk PlotDenoisedCopyRatios --standardized-copy-ratios "${S}.standardizedCR.tsv" \
  --denoised-copy-ratios "${S}.denoisedCR.tsv" --sequence-dictionary "${DICT}" --output-prefix "${S}" -O plots
gatk PlotModeledSegments --denoised-copy-ratios "${S}.denoisedCR.tsv" --allelic-counts "segments/${S}.hets.tsv" \
  --segments "segments/${S}.modelFinal.seg" --sequence-dictionary "${DICT}" --output-prefix "${S}" -O plots
```

Notes and defaults to justify:

| Item | Default | Why / when to change |
|---|---|---|
| PoN size and composition | ≥ 20 normals, same kit/protocol, sex-aware (build one PoN per sex or exclude allosomes) | Fewer normals leave GC/batch noise; mixed-sex PoN creates chrX artefacts |
| `--annotated-intervals` in PoN | include when `AnnotateIntervals` was run | Enables explicit GC correction in denoising |
| `CreateReadCountPanelOfNormals --number-of-eigensamples` | tool default (20) | Lower (5–10) for small PoNs to avoid over-fitting |
| `ModelSegments --number-of-changepoints-penalty-factor` | 1.0 | Raise (2–4) for FFPE/noisy WES to reduce hypersegmentation |
| `ModelSegments --minimum-total-allele-count-case` | 0 (uses normal to pick hets) | Tumour-only: set `--minimum-total-allele-count-case 30` and drop `--normal-allelic-counts` |
| Common SNP list | gnomAD/1000G biallelic SNPs AF 0.05–0.95 (`genomics/reference.md`) | Sparse lists starve BAF segmentation on WES |
| `PlotDenoisedCopyRatios` `.deltaMAD` / `.scaledDeltaMAD` | denoised MAD < standardized MAD | If denoising does not lower MAD, PoN is mismatched |

Outputs: `segments/${S}.modelFinal.seg` (CR + minor-allele fraction per segment), `segments/${S}.called.seg` (`+`/`-`/`0` calls on copy ratio only). Purity/ploidy and absolute copy number from these segments: `.claude/skills/genomics-somatic/SKILL.md`. Validate: `plots/${S}.denoised.png` shows chromosome-arm-level structure, not periodic GC waves; `segments/${S}.hets.tsv` row count ≥ 1,000 (WES) or ≥ 100,000 (WGS).

## 3. CNVkit step-by-step

Equivalent to `cnvkit.py batch` (hybrid capture shown). Use when you need to persist a pooled reference, change bin sizes or inspect intermediate files.

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TARGETS=targets.bed
REFFLAT=resources/refFlat.GRCh38.txt
ACCESS=resources/access-5kb-mappable.GRCh38.bed
EXCLUDE=resources/GRCh38_lowmappability.bed
THREADS=8
cnvkit.py access "${REF}" -x "${EXCLUDE}" -o "${ACCESS}"
cnvkit.py autobin normal/*.bam tumor/*.bam -t "${TARGETS}" -g "${ACCESS}" --annotate "${REFFLAT}" --short-names
for BAM in normal/*.bam tumor/*.bam; do
  S=$(basename "${BAM}" .bam)
  cnvkit.py coverage "${BAM}" targets.target.bed -p "${THREADS}" -o "cov/${S}.targetcoverage.cnn"
  cnvkit.py coverage "${BAM}" targets.antitarget.bed -p "${THREADS}" -o "cov/${S}.antitargetcoverage.cnn"
done
cnvkit.py reference cov/normal*.targetcoverage.cnn cov/normal*.antitargetcoverage.cnn --fasta "${REF}" -o pooled_reference.cnn
for BAM in tumor/*.bam; do
  S=$(basename "${BAM}" .bam)
  cnvkit.py fix "cov/${S}.targetcoverage.cnn" "cov/${S}.antitargetcoverage.cnn" pooled_reference.cnn -o "cnvkit/${S}.cnr"
  cnvkit.py segment "cnvkit/${S}.cnr" -m cbs --drop-low-coverage -p "${THREADS}" -o "cnvkit/${S}.cns"
  cnvkit.py segmetrics "cnvkit/${S}.cnr" -s "cnvkit/${S}.cns" --ci --sem -o "cnvkit/${S}.ci.cns"
  cnvkit.py call "cnvkit/${S}.ci.cns" --center median --filter ci -m threshold -o "cnvkit/${S}.call.cns"
  cnvkit.py export seg "cnvkit/${S}.call.cns" -o "cnvkit/${S}.seg"
  cnvkit.py export vcf "cnvkit/${S}.call.cns" -i "${S}" -o "cnvkit/${S}.cnv.vcf"
done
cnvkit.py metrics cnvkit/*.cnr -s cnvkit/*.cns > cnvkit/metrics.tsv
cnvkit.py heatmap cnvkit/*.call.cns -d -o cnvkit/heatmap.pdf
```

Variants:

| Situation | Change |
|---|---|
| Germline samples | `cnvkit.py segment` without `--drop-low-coverage`; `call -m threshold -t=-1.1,-0.4,0.3,0.7` or `-m clonal` |
| Male reference / sex handling | add `-y` to `reference`, `fix`, `call`, `export` consistently; `cnvkit.py sex *.cnr` to check inferred sex |
| WGS | `autobin -m wgs -b 50000 -g ${ACCESS} --annotate ${REFFLAT}`; no antitargets; `reference --no-edge`; `fix --no-edge`; `segment -t 1e-6` or `-m hmm` |
| Amplicon | `target ${TARGETS} --split -o targets.split.bed`; `coverage ... -p 0`; empty antitarget file (`touch MT`); `reference --no-edge`; `fix ... MT ref.cnn --no-edge`; **unmarked duplicates** |
| Tumour-only, no normals | `cnvkit.py reference -o flat_reference.cnn -f ${REF} -t targets.target.bed -a targets.antitarget.bed` (flat); report reduced sensitivity |
| BAF and purity | `bcftools mpileup -f ${REF} -a AD,DP ${TUMOR} \| bcftools call -mv -Oz -o ${S}.vcf.gz`; `segment --vcf ${S}.vcf.gz -i ${S}`; `call -m clonal --purity 0.6 --ploidy 2 -v ${S}.vcf.gz -i ${S}` |
| Known FFPE noise | `fix --no-edge`; `segment -m hmm-tumor` or raise `-t`; accept MAD up to ~0.4 |

Validate: `cnvkit/metrics.tsv` MAD/`bivar` per sample (SKILL.md QC table); `cnvkit.py breaks cnvkit/${S}.cnr cnvkit/${S}.call.cns` lists genes split by a breakpoint; `cnvkit.py genemetrics cnvkit/${S}.cnr -s cnvkit/${S}.call.cns -t 0.2 -m 3` lists gene-level gains/losses.

## 4. SV callers

### 4a. smoove population workflow (N > ~40; parallel per sample)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
EXCLUDE=resources/exclude.cnvnator_100bp.GRCh38.20170403.bed
COHORT=cohort01
for BAM in bams/*.bam; do
  S=$(basename "${BAM}" .bam)
  smoove call --outdir results-smoove --exclude "${EXCLUDE}" --name "${S}" --fasta "${REF}" -p 1 --genotype "${BAM}"
done
smoove merge --name "${COHORT}" -f "${REF}" --outdir ./ results-smoove/*.genotyped.vcf.gz
for BAM in bams/*.bam; do
  S=$(basename "${BAM}" .bam)
  smoove genotype -d -x -p 1 --name "${S}-joint" --outdir results-genotyped --fasta "${REF}" --vcf "${COHORT}.sites.vcf.gz" "${BAM}"
done
smoove paste --name "${COHORT}" results-genotyped/*.vcf.gz
smoove annotate --gff resources/gencode.GRCh38.gff3.gz "${COHORT}.smoove.square.vcf.gz" | bgzip -c > "${COHORT}.smoove.square.anno.vcf.gz"
tabix -p vcf "${COHORT}.smoove.square.anno.vcf.gz"
```

`genotype -d` runs duphold (adds `DHFFC`, `DHBFC`, `DHSP`); `-x` removes `PRPOS`/`PREND`. Small cohorts (< 40) can use one `smoove call -x --name ${COHORT} --exclude ... --genotype --duphold bams/*.bam`. Set `TMPDIR` to a large disk. First-pass filters: `INFO/MSHQ > 3`; DEL `DHFFC < 0.7`; DUP `DHBFC > 1.3`. Validate: `bcftools query -l ${COHORT}.smoove.square.vcf.gz | wc -l` equals sample count; `bcftools view -H ... | wc -l` identical across `results-genotyped/*.vcf.gz`.

### 4b. DELLY (v2.x subcommands: `sr`, `lr`, `cnv`, `merge`, `filter`)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
EXCL=resources/human.hg38.excl.tsv
S=sample
TUMOR=tumor.bam
NORMAL=normal.bam
delly sr -g "${REF}" -x "${EXCL}" -o "${S}.delly.bcf" "bams/${S}.bam"
delly sr -g "${REF}" -x "${EXCL}" -o tn.bcf "${TUMOR}" "${NORMAL}"
printf 'tumor\ttumor\nnormal\tcontrol\n' > samples.tsv
delly filter -f somatic -o tn.somatic.bcf -s samples.tsv tn.bcf
delly merge -o sites.bcf bams_delly/*.bcf
delly sr -g "${REF}" -v sites.bcf -o "${S}.geno.bcf" "bams/${S}.bam"
bcftools merge -m id -Ob -o cohort.geno.bcf geno/*.geno.bcf
delly filter -f germline --hwe 0.000001 -o cohort.germline.bcf cohort.geno.bcf
delly cnv -g "${REF}" -m resources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz -c "${S}.cov.gz" -o "${S}.cnv.bcf" -u "${S}.seg.bed" "bams/${S}.bam"
```

`human.hg38.excl.tsv` ships in the DELLY repo (`excludeTemplates/`); sample IDs in `samples.tsv` must equal the BAM `@RG SM` values. `delly cnv` is WGS-only; mappability map from <https://gear-genomics.embl.de/data/delly/> (`<check current>` for the exact filename). `delly call` is the pre-2.0 name of `delly sr`. Validate: `bcftools view -H ${S}.delly.bcf | cut -f7 | sort | uniq -c`; somatic BCF has `FORMAT/FT` PASS in tumour only.

### 4c. Tiddit

```bash
set -euo pipefail
REF=ref/GRCh38.fa
S=sample
tiddit --sv --bam "bams/${S}.bam" --ref "${REF}" --threads 8 -o "tiddit/${S}"
grep -E '^#|PASS' "tiddit/${S}.vcf" > "tiddit/${S}.pass.vcf"
tiddit --cov --bam "bams/${S}.bam" --ref "${REF}" -z 10000 -o "tiddit/${S}.cov"
```

Needs `bwa` + `fermi2` + `ropebwt2` for local assembly (or `--skip_assembly`). Tumour: `--p_ratio 0.10`. Merge across callers/samples with SVDB (`svdb --merge --vcf a.vcf b.vcf --bnd_distance 500 --overlap 0.6`), which is what nf-core/raredisease does. Validate: `tiddit/${S}.ploidies.tab` shows ~2 on autosomes and expected X/Y.

### 4d. dysgu

```bash
set -euo pipefail
REF=ref/GRCh38.fa
S=sample
dysgu run -p 4 --clean "${REF}" "dysgu_tmp_${S}" "bams/${S}.bam" > "dysgu/${S}.vcf"
dysgu filter "dysgu/${S}.vcf" > "dysgu/${S}.filt.vcf"
dysgu filter --normal-vcf "dysgu/normal.vcf" "dysgu/tumor.vcf" > "dysgu/tumor.somatic.vcf"
```

Diploid model by default; `--diploid False` for tumours with CNAs; `dysgu merge` for cohorts. Quality in `INFO/PROB` (ML probability) and `FILTER`. Temp dir needs > 5 GB per 30× genome. Validate: `grep -vc '^#' dysgu/${S}.filt.vcf` in the 5,000–12,000 range for 30× WGS.

### 4e. duphold standalone (annotate any DEL/DUP VCF with depth support)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
S=sample
duphold -t 4 -v "${S}.manta.vcf" -b "bams/${S}.bam" -f "${REF}" -o "${S}.manta.duphold.vcf"
bcftools view -i '(INFO/SVTYPE="DEL" && FMT/DHFFC[0]<0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0]>1.3) || (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")' \
  -Oz -o "${S}.manta.duphold.filt.vcf.gz" "${S}.manta.duphold.vcf"
```

Optional `-s snps.vcf.gz` adds het-site allele-balance (`DHBFC` uses GC-matched bins, `DHFFC` flanks 1 kb each side). Published performance: DEL `DHFFC < 0.7` keeps ~98–99% of GIAB true DELs and removes 60–85% of false ones.

### 4f. GRIDSS (assembly-based; somatic breakpoints)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
GRIDSS_JAR=/opt/gridss/gridss.jar
gridss --reference "${REF}" --output gridss/tn.vcf.gz --assembly gridss/tn.assembly.bam --threads 8 \
  --jvmheap 30g --workingdir gridss/work --jar "${GRIDSS_JAR}" normal.bam tumor.bam
gridss_somatic_filter --pondir resources/gridss_pon --input gridss/tn.vcf.gz --output gridss/tn.somatic.vcf.gz \
  --fulloutput gridss/tn.full.vcf.gz --scriptdir /opt/gridss -n 1 -t 2
```

Order of BAMs sets sample indices (`-n 1 -t 2` = normal first). Use the HMF PON for GRCh38 `<check current: https://github.com/hartwigmedical/hmftools>`; in practice run GRIDSS through PURPLE/LINX via nf-core/oncoanalyser (`.claude/skills/genomics-somatic/SKILL.md`).

### 4g. Multi-caller / multi-sample merging

```bash
set -euo pipefail
S=sample
printf '%s\n' "${S}.manta.vcf" "${S}.smoove.vcf" "${S}.dysgu.filt.vcf" > "${S}.vcf_list.txt"
SURVIVOR merge "${S}.vcf_list.txt" 1000 2 1 1 0 50 "${S}.sv.merged.vcf"
jasmine file_list="${S}.vcf_list.txt" out_file="${S}.sv.jasmine.vcf" max_dist=500 min_support=2 --output_genotypes --normalize_type --dup_to_ins threads=8
SURVIVOR stats "${S}.sv.merged.vcf" -1 -1 -1 "${S}.sv.merged.stats"
```

SURVIVOR arguments: max breakpoint distance (bp), minimum supporting callers, require same type (1), require same strand (1), estimate distance (0), minimum size (bp). SURVIVOR needs uncompressed VCFs and keeps the first file's genotype; Jasmine is breakpoint-distance-aware and handles INS better (`--dup_to_ins` treats DUP as INS when merging with insertion callers). Validate: `grep -c SUPP=2 ${S}.sv.merged.vcf` vs total; the two mergers should agree within ~10% on DEL/DUP counts.

## 5. ExomeDepth (raw R; what DECoN wraps)

```r
suppressPackageStartupMessages(library(ExomeDepth))
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]        # 4 columns: chromosome start end name, same contig names as BAMs
bam_list <- args[2]        # text file, one BAM path per line
fasta    <- args[3]
out_dir  <- args[4]
bed <- read.table(bed_file, header = FALSE, col.names = c("chromosome", "start", "end", "name"))
bams <- readLines(bam_list)
counts <- getBamCounts(bed.frame = bed, bam.files = bams, include.chr = FALSE, referenceFasta = fasta)
counts <- as.data.frame(counts)
count_mat <- as.matrix(counts[, grep("\\.bam$", colnames(counts))])
colnames(count_mat) <- sub("\\.bam$", "", basename(bams))
bin_len <- (counts$end - counts$start) / 1000
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (i in seq_len(ncol(count_mat))) {
  test_name <- colnames(count_mat)[i]
  choice <- select.reference.set(test.counts = count_mat[, i], reference.counts = count_mat[, -i, drop = FALSE],
                                 bin.length = bin_len, n.bins.reduced = 10000)
  ref_sum <- rowSums(count_mat[, choice$reference.choice, drop = FALSE])
  ed <- new("ExomeDepth", test = count_mat[, i], reference = ref_sum, formula = "cbind(test, reference) ~ 1")
  ed <- CallCNVs(x = ed, transition.probability = 1e-4, chromosome = counts$chromosome,
                 start = counts$start, end = counts$end, name = counts$name)
  calls <- ed@CNV.calls
  calls$sample <- test_name
  calls$n_ref_samples <- length(choice$reference.choice)
  calls$ref_correlation <- max(choice$summary.stats$correlations)
  write.csv(calls[order(-calls$BF), ], file.path(out_dir, paste0(test_name, ".exomedepth.csv")), row.names = FALSE)
}
```

Run: `Rscript exomedepth.R targets.4col.bed bams.txt ref/GRCh38.fa exomedepth/`. Interpretation defaults: Bayes factor `BF ≥ 10` and `reads.ratio` < 0.75 (del) or > 1.25 (dup); reference set of 5–20 samples with correlation ≥ 0.98; `transition.probability` 1e-4 for exomes, 1e-2 (DECoN default `--transProb 0.01`) for small panels where single-exon sensitivity matters. Sex: build reference sets within sex or exclude chrX/chrY; ExomeDepth has no ploidy model. Validate: per-sample call count within 2× cohort median; known positive controls recovered with `BF ≥ 10`.

## 6. Resource files referenced above

| File | Source | Notes |
|---|---|---|
| Mappability track for `AnnotateIntervals` | Umap k100 single-read mappability (GRCh38) `<check current: https://bismap.hoffmanlab.org/>` | Convert to BED, bgzip, tabix; GATK needs `--mappability-track` as an indexed BED |
| Segmental duplications | UCSC `genomicSuperDups` table for hg38, exported as BED | bgzip + tabix |
| CNV blacklist / exclude | GATK: centromeres + gaps as `-XL`; smoove: `exclude.cnvnator_100bp.GRCh38.20170403.bed` (hall-lab/speedseq); DELLY: `excludeTemplates/human.hg38.excl.tsv` | Different formats; do not interchange |
| Common SNP sites (`CollectAllelicCounts`) | gnomAD v4.1 or 1000G biallelic SNPs, AF 0.05–0.95, as `.interval_list` | Download mechanics: `.claude/skills/genomics/reference.md` |
| `refFlat.txt` (CNVkit `--annotate`) | UCSC `refFlat` for hg38 | Gene labels for bins |
| SVAFotate population BED | `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz` (Zenodo record 11642574) | Aggregates gnomAD-SV, CCDG, 1000G, TOPMed |
| gnomAD SVs v4.1.0 / CNVs v4.1.0 | <https://gnomad.broadinstitute.org/downloads> | SV = genomes; CNV = exome depth calls; both GRCh38 |
| ClinGen dosage sensitivity | `<check current: https://ftp.clinicalgenome.org/>` (gene and region curation lists, GRCh38) | Record download date with every ACMG score |
| AnnotSV annotations | `make PREFIX=. install-human-annotation` in the AnnotSV install | Bundle date goes in `results/versions.txt` |
| GIAB HG002 SV truth | `*_HG002-T2TQ100-V1.1_stvar.vcf.gz` + `.benchmark.bed` under the NIST GIAB FTP | Benchmarking: `.claude/skills/genomics-benchmark/SKILL.md` |
