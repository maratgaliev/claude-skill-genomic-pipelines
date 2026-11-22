---
name: genomics-cnv
description: Short-read copy-number (CNV) and structural-variant (SV) calling, filtering, annotation and ACMG/ClinGen CNV scoring for germline WES/WGS/panels and somatic copy number. Use when the user mentions GATK gCNV/ModelSegments, CNVkit, ExomeDepth/DECoN, CNVpytor, Manta, DELLY, smoove/duphold, Tiddit, dysgu, SURVIVOR/Jasmine, AnnotSV, SVAFotate, .cns/.cnr/.seg files or "deletion/duplication/translocation". Not for long-read SV callers (genomics-longread), somatic SNV/indel or purity/ploidy modelling (genomics-somatic), or truth-set benchmarking mechanics (genomics-benchmark).
argument-hint: "[germline-wes|germline-wgs|somatic|sv|annotate] [path]"
paths:
  - "**/*.cns"
  - "**/*.cnr"
  - "**/*.cnn"
  - "**/*.seg"
  - "**/*.counts.hdf5"
  - "**/diploidSV.vcf.gz"
  - "**/somaticSV.vcf.gz"
  - "**/*-smoove.genotyped.vcf.gz"
  - "**/*.annotated.tsv"
---

# Genomics: CNV and SV (short read)

Expert in read-depth CNV and paired/split-read SV detection from Illumina WES, WGS and panel data: germline cohort and case calling, somatic copy-ratio/BAF segmentation, SV filtering with caller-specific evidence fields, population-frequency and dosage-sensitivity annotation, and ACMG/ClinGen 2020 CNV scoring. Outputs are exact commands, decision tables and validation steps; thresholds are defaults to justify, not truths.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" gatk cnvkit.py configManta.py delly smoove duphold tiddit dysgu cnvpytor SURVIVOR jasmine AnnotSV svafotate bcftools samtools mosdepth samplot Rscript || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| Question | Where | Command |
|---|---|---|
| Assembly, contig naming, alt/decoy contigs | BAM/CRAM `@SQ` | `samtools view -H ${BAM} \| grep -c '^@SQ'`; `peek.sh ${BAM}` |
| WES vs WGS vs panel | `@PG`/`@CO` lines, target BED in project dir, read density | `samtools idxstats ${BAM} \| head -3`; `mosdepth -n -x --by 1000 q ${BAM}` |
| Duplicates marked (required by Manta/DELLY/GATK; must NOT be marked for amplicon CNVkit) | `@PG ID:MarkDuplicates`/`samblaster` | `samtools view -H ${BAM} \| grep -i -c -E 'MarkDuplicates\|samblaster\|sambamba'` |
| Insert size, paired, read length | first reads | `samtools stats ${BAM} 2>/dev/null \| grep -E '^SN\s+(insert size average\|average length\|reads paired)'` |
| Sex (needed for PoN matching, allosomal ploidy priors) | chrX/chrY depth | `samtools idxstats ${BAM} \| awk '$1=="chrX"\|\|$1=="chrY"'` |
| Cohort size / PoN availability | sample sheet, prior `*.counts.hdf5`, `*.cnn` | `ls *.counts.hdf5 *.cnn 2>/dev/null \| wc -l` |
| Tumour/normal pairing, purity estimate | metadata, pathology | sample sheet `status` column |
| Prior calls | `*.cns`, `*.seg`, `diploidSV.vcf.gz`, `*.bcf` | `peek.sh <file> 20`; `bcftools view -h x.vcf.gz \| grep -E '^##(source\|INFO=<ID=SVTYPE)'` |

Ask only if still ambiguous: capture kit BED when none is on disk; tumour purity when no normal and no pathology estimate; whether the question is diagnostic (ACMG scoring needed) or discovery; consent/data-use limits on sharing calls with web tools.

## Decision table

| Data / question | First choice | Alternative / notes |
|---|---|---|
| WES germline, clinical exon-level CNV, ≥10 samples same kit | ExomeDepth or DECoN (ExomeDepth wrapper with QC + plots) | GATK gCNV cohort (≥30 samples, better for large cohorts); CNVkit hybrid as second caller |
| WES germline, single sample, existing gCNV model | GATK gCNV CASE mode against cohort model | ExomeDepth with a stored reference count matrix |
| WGS germline CNV (read depth) | GATK gCNV (`--bin-length 1000`) or CNVpytor (fast, single sample) | DELLY `cnv`; CNVkit `--method wgs` for quick looks |
| WGS germline SV (breakpoints) | Manta + smoove (`--genotype --duphold`); merge with SURVIVOR/Jasmine, require 2 callers | dysgu (ML-scored, single binary); Tiddit (nf-core/raredisease default); DELLY `sr` |
| WGS germline SV, cohort / population | smoove population workflow (call → merge → genotype → paste) | DELLY merge/genotype/filter chain |
| Somatic WES/WGS copy number, matched normal | GATK somatic CNV (ModelSegments with allelic counts) | CNVkit hybrid/wgs with pooled normal reference |
| Somatic, tumour-only | CNVkit with PoN reference; GATK somatic CNV with PoN | Flat reference only when no normals exist; expect noise |
| Tumour purity/ploidy and allele-specific CN | ASCAT / PURPLE / FACETS → `.claude/skills/genomics-somatic/SKILL.md` | one decision row here; full treatment there |
| Somatic SV | Manta tumour/normal (`somaticSV.vcf.gz`) | DELLY `sr` + `delly filter -f somatic`; GRIDSS for assembly-level breakpoints |
| Panel / amplicon (no off-target reads) | CNVkit `--method amplicon` (unmarked duplicates) | ExomeDepth/DECoN for exon-level panels |
| Low-pass WGS (<1x) | CNVpytor or QDNAseq at 50-100 kb bins | not GATK gCNV defaults |
| Long-read SV (Sniffles2, cuteSV, pbsv, Severus) | `.claude/skills/genomics-longread/SKILL.md` | DELLY `lr` also documented there |
| STR expansions | ExpansionHunter → `.claude/skills/genomics-annotation/rare-disease.md` | MEI: MELT / xTea (below) |

## Workflow

Full command dumps: `reference.md` in this directory. Below are canonical happy paths.

### 1. Interval preparation (GATK gCNV and somatic CNV share this)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TARGETS=targets.bed
MAPPABILITY_BED=resources/k100.umap.GRCh38.bed.gz
SEGDUP_BED=resources/genomicSuperDups.GRCh38.bed.gz
MODE=WES
if [ "${MODE}" = WES ]; then BIN=0; PAD=250; else BIN=1000; PAD=0; fi
gatk PreprocessIntervals -R "${REF}" -L "${TARGETS}" --bin-length "${BIN}" --padding "${PAD}" \
  --interval-merging-rule OVERLAPPING_ONLY -O intervals.preprocessed.interval_list
gatk AnnotateIntervals -R "${REF}" -L intervals.preprocessed.interval_list \
  --mappability-track "${MAPPABILITY_BED}" --segmental-duplication-track "${SEGDUP_BED}" \
  --interval-merging-rule OVERLAPPING_ONLY -O intervals.annotated.tsv
```

WGS: drop `-L` and pass `-XL blacklist.bed` (centromeres, gaps). If you have no mappability/segdup BEDs, drop both `--*-track` flags **and** the `--minimum-mappability`/`--maximum-segmental-duplication-content` flags in `FilterIntervals` (GC filtering still works). Validate: `grep -vc '^@' intervals.preprocessed.interval_list` (WES ≈ number of targets; WGS ≈ genome/1000).

### 2. Germline CNV — GATK gCNV cohort (≥30 samples; WES or WGS)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
PRIORS=ploidy_priors.tsv
COHORT=cohort01
for BAM in bams/*.bam; do
  S=$(basename "${BAM}" .bam)
  gatk CollectReadCounts -I "${BAM}" -R "${REF}" -L intervals.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O "counts/${S}.counts.hdf5"
done
INPUTS=$(printf -- '-I %s ' counts/*.counts.hdf5)
gatk FilterIntervals -L intervals.preprocessed.interval_list --annotated-intervals intervals.annotated.tsv ${INPUTS} \
  --minimum-gc-content 0.1 --maximum-gc-content 0.9 --minimum-mappability 0.9 \
  --maximum-segmental-duplication-content 0.5 --low-count-filter-count-threshold 10 \
  --low-count-filter-percentage-of-samples 50 --extreme-count-filter-minimum-percentile 1 \
  --extreme-count-filter-maximum-percentile 99 --extreme-count-filter-percentage-of-samples 90 \
  --interval-merging-rule OVERLAPPING_ONLY -O intervals.filtered.interval_list
gatk DetermineGermlineContigPloidy -L intervals.filtered.interval_list ${INPUTS} --contig-ploidy-priors "${PRIORS}" \
  --interval-merging-rule OVERLAPPING_ONLY -O ploidy --output-prefix "${COHORT}"
gatk GermlineCNVCaller --run-mode COHORT -L intervals.filtered.interval_list ${INPUTS} \
  --contig-ploidy-calls "ploidy/${COHORT}-calls" --annotated-intervals intervals.annotated.tsv \
  --interval-merging-rule OVERLAPPING_ONLY -O gcnv --output-prefix "${COHORT}"
```

Ploidy priors file (tab-separated; every contig in the interval list must appear):

```text
CONTIG_NAME	PLOIDY_PRIOR_0	PLOIDY_PRIOR_1	PLOIDY_PRIOR_2	PLOIDY_PRIOR_3
chr1	0.01	0.01	0.97	0.01
chr22	0.01	0.01	0.97	0.01
chrX	0.01	0.49	0.49	0.01
chrY	0.50	0.50	0.00	0.00
```

Validate: `ls gcnv/${COHORT}-calls/SAMPLE_0/` and `grep -c PASS ploidy/${COHORT}-calls/SAMPLE_0/contig_ploidy.tsv`; inspect `ploidy/*-calls/SAMPLE_*/contig_ploidy.tsv` for sex concordance with metadata. Postprocess per sample and CASE mode: `reference.md` §1. WGS cohorts: scatter `-L` into shards of ~5,000–10,000 intervals and run one `GermlineCNVCaller` per shard (same `--output-prefix`, separate `-O`), then give all shard dirs to `PostprocessGermlineCNVCalls`.

### 3. Germline WES CNV — ExomeDepth / DECoN (clinical exon-level standard)

DECoN wraps ExomeDepth with sample-failure QC and plots; run from its install directory (uses `renv`). BED must be 4 columns (chr, start, end, name) with contig names matching the BAM.

```bash
set -euo pipefail
DECON_DIR=/opt/DECoN/Linux
BAMS=bams.txt
BED=targets.4col.bed
REF=ref/GRCh38.fa
cd "${DECON_DIR}"
Rscript ReadInBams.R --bams "${BAMS}" --bed "${BED}" --fasta "${REF}" --out decon/counts
Rscript IdentifyFailures.R --RData decon/counts.RData --mincorr 0.98 --mincov 100 --out decon/qc
Rscript makeCNVcalls.R --RData decon/counts.RData --transProb 0.01 --out decon/calls --plot None
```

Validate: `peek.sh decon/qc_Failures.txt`; `wc -l decon/calls_all.txt`; per-sample call counts within 2× cohort median. Raw ExomeDepth R: `reference.md` §5.

### 4. Germline WGS CNV — CNVpytor (single sample, fast)

```bash
set -euo pipefail
BAM=sample.bam
SAMPLE=sample
REF=ref/GRCh38.fa
cnvpytor -root "${SAMPLE}.pytor" -rd "${BAM}" -T "${REF}" -chrom chr{1..22} chrX chrY
cnvpytor -root "${SAMPLE}.pytor" -his 1000 10000 100000
cnvpytor -root "${SAMPLE}.pytor" -partition 1000 10000 100000
cnvpytor -root "${SAMPLE}.pytor" -call 10000 > "${SAMPLE}.cnvpytor.10k.tsv"
```

Validate: `awk '$5<0.0001 && $9<0.5' "${SAMPLE}.cnvpytor.10k.tsv" | wc -l` (column 5 = e-val1, column 9 = q0 fraction of MAPQ0 reads; keep e-val1 < 1e-4 and q0 < 0.5). Bin 1,000 for 30x+, 10,000 for 5-10x, 100,000 for low-pass.

### 5. Somatic copy number — GATK ModelSegments (matched normal)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TUMOR=tumor.bam
NORMAL=normal.bam
PON=cnv_pon.hdf5
SNPS=resources/common_biallelic_snps.GRCh38.interval_list
S=tumor
gatk CollectReadCounts -I "${TUMOR}" -R "${REF}" -L intervals.preprocessed.interval_list --interval-merging-rule OVERLAPPING_ONLY -O "${S}.counts.hdf5"
gatk DenoiseReadCounts -I "${S}.counts.hdf5" --count-panel-of-normals "${PON}" \
  --standardized-copy-ratios "${S}.standardizedCR.tsv" --denoised-copy-ratios "${S}.denoisedCR.tsv"
gatk CollectAllelicCounts -I "${TUMOR}" -R "${REF}" -L "${SNPS}" -O "${S}.allelicCounts.tsv"
gatk CollectAllelicCounts -I "${NORMAL}" -R "${REF}" -L "${SNPS}" -O normal.allelicCounts.tsv
gatk ModelSegments --denoised-copy-ratios "${S}.denoisedCR.tsv" --allelic-counts "${S}.allelicCounts.tsv" \
  --normal-allelic-counts normal.allelicCounts.tsv --output-prefix "${S}" -O segments
gatk CallCopyRatioSegments -I "segments/${S}.cr.seg" -O "segments/${S}.called.seg"
```

Validate: `grep -vc '^@' segments/${S}.modelFinal.seg` (WES tumours typically 50–500 segments); `PlotModeledSegments` output visually flat in known-diploid regions. PoN build and plots: `reference.md` §2. Purity/ploidy fitting (ASCAT/PURPLE) and absolute CN: `.claude/skills/genomics-somatic/SKILL.md`.

### 6. CNVkit (WES hybrid / WGS / amplicon; germline or somatic)

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TARGETS=targets.bed
REFFLAT=resources/refFlat.GRCh38.txt
ACCESS=resources/access-5kb-mappable.GRCh38.bed
METHOD=hybrid
cnvkit.py access "${REF}" -o "${ACCESS}"
cnvkit.py batch tumor/*.bam --normal normal/*.bam --method "${METHOD}" --targets "${TARGETS}" \
  --annotate "${REFFLAT}" --fasta "${REF}" --access "${ACCESS}" --drop-low-coverage \
  --output-reference pooled_reference.cnn --output-dir cnvkit -p 8 --scatter --diagram
S=tumor1
cnvkit.py segmetrics "cnvkit/${S}.cnr" -s "cnvkit/${S}.cns" --ci -o "cnvkit/${S}.ci.cns"
cnvkit.py call "cnvkit/${S}.ci.cns" --center median --filter ci -m threshold -o "cnvkit/${S}.call.cns"
cnvkit.py export seg "cnvkit/${S}.call.cns" -o "cnvkit/${S}.seg"
cnvkit.py export vcf "cnvkit/${S}.call.cns" -i "${S}" -o "cnvkit/${S}.cnv.vcf"
```

Semantics to get right: `--normal` with BAMs builds a pooled reference; `--normal` with **no** files builds a flat reference (tumour-only fallback, noisy); `-r existing.cnn` reuses a reference and skips normals. `--method wgs` needs `--annotate` and no `--targets`; `--method amplicon` needs `--targets` and **unmarked** duplicates. `--drop-low-coverage` is for tumours (drops zero-depth bins that mimic homozygous deletions); omit for germline. Known purity: `cnvkit.py call -m clonal --purity 0.6 -v tumor.vcf`. Validate: `cnvkit.py metrics cnvkit/*.cnr -s cnvkit/*.cns` (see QC table).

### 7. SV calling — Manta + smoove, merge, depth-annotate

```bash
set -euo pipefail
REF=ref/GRCh38.fa
BAM=sample.bam
SAMPLE=sample
EXCLUDE=resources/exclude.cnvnator_100bp.GRCh38.20170403.bed
REGIONS=resources/primary_chroms.bed.gz
MANTA_DIR=/opt/manta
configManta.py --bam "${BAM}" --referenceFasta "${REF}" --callRegions "${REGIONS}" --runDir "manta_${SAMPLE}"
"manta_${SAMPLE}/runWorkflow.py" -j 8
python "${MANTA_DIR}/libexec/convertInversion.py" "$(command -v samtools)" "${REF}" \
  "manta_${SAMPLE}/results/variants/diploidSV.vcf.gz" > "${SAMPLE}.manta.vcf"
smoove call -x --name "${SAMPLE}" --exclude "${EXCLUDE}" --fasta "${REF}" -p 2 --genotype --duphold --outdir smoove "${BAM}"
bcftools view -Ov -o "${SAMPLE}.smoove.vcf" "smoove/${SAMPLE}-smoove.genotyped.vcf.gz"
printf '%s\n' "${SAMPLE}.manta.vcf" "${SAMPLE}.smoove.vcf" > vcf_list.txt
SURVIVOR merge vcf_list.txt 1000 2 1 1 0 50 "${SAMPLE}.sv.2callers.vcf"
```

WES: add `--exome` to `configManta.py` (and `--callRegions` = padded target BED); do not run smoove/DELLY/Tiddit on exomes. Tumour/normal: `configManta.py --tumorBam T.bam --normalBam N.bam` → `somaticSV.vcf.gz`. `convertInversion.py` turns Manta's 4-breakend inversions back into `<INV>` records (needed by most mergers/annotators); SURVIVOR needs uncompressed VCFs. Jasmine alternative: `jasmine file_list=vcf_list.txt out_file=merged.vcf max_dist=500 min_support=2 --output_genotypes threads=8`. Validate: `cut -f7 "${SAMPLE}.manta.vcf" | grep -v '^#' | sort | uniq -c` and `bcftools query -f '%SVTYPE\n' smoove/*.genotyped.vcf.gz | sort | uniq -c`. Population workflow, DELLY, Tiddit, dysgu, duphold-standalone: `reference.md` §3–4.

### 8. SV filtering — caller-specific fields

Manta `PR`/`SR` are **FORMAT** fields (`ref,alt` pairs per sample); DELLY keeps `PE`/`SR` in INFO and per-sample `DV`/`RV`/`DR`/`RR`/`GQ` in FORMAT; smoove/duphold adds FORMAT `DHFFC`/`DHBFC`; dysgu uses `PROB`.

```bash
set -euo pipefail
SAMPLE=sample
bcftools view -i 'FILTER="PASS" && (FMT/SR[0:1]>=2 || FMT/PR[0:1]>=3) && (INFO/SVTYPE="BND" || abs(INFO/SVLEN)>=50)' \
  -Oz -o "${SAMPLE}.manta.filt.vcf.gz" "${SAMPLE}.manta.vcf"
bcftools view -i 'FILTER="PASS" && INFO/PRECISE=1 && (INFO/PE>=3 || INFO/SR>=2) && FMT/GQ[0]>=15 && FMT/DV[0]+FMT/RV[0]>=3' \
  -Ob -o "${SAMPLE}.delly.filt.bcf" "${SAMPLE}.delly.bcf"
bcftools view -i '(INFO/SVTYPE="DEL" && FMT/DHFFC[0]<0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0]>1.3) || (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")' \
  -Oz -o "${SAMPLE}.smoove.filt.vcf.gz" "smoove/${SAMPLE}-smoove.genotyped.vcf.gz"
```

Validate: `bcftools stats` before/after; expect DEL/DUP counts to fall 40–70% after duphold, INS/BND unchanged.

### 9. Annotation and scoring

```bash
set -euo pipefail
SAMPLE=sample
SVAF_BED=resources/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz
svafotate annotate -v "${SAMPLE}.sv.2callers.vcf" -o "${SAMPLE}.sv.af.vcf" -b "${SVAF_BED}" -f 0.8 -a best
AnnotSV -SVinputFile "${SAMPLE}.sv.af.vcf" -genomeBuild GRCh38 -annotationMode both -outputDir annotsv -outputFile "${SAMPLE}.annotsv.tsv"
AnnotSV -SVinputFile "${SAMPLE}.cnv.bed" -svtBEDcol 4 -samplesidBEDcol 5 -genomeBuild GRCh38 -outputDir annotsv -outputFile "${SAMPLE}.cnv.annotsv.tsv"
```

`-svtBEDcol` only applies to BED input (column holding DEL/DUP); VCF input reads `SVTYPE` and must not carry it. AnnotSV `ACMG_class` (1–5) implements the Riggs 2020 point scheme using ClinGen dosage, gnomAD-SV, DGV, OMIM; borderline calls can be recomputed in the ClinGen CNV calculator (<https://cnvcalc.clinicalgenome.org/>), but that is a web form receiving patient CNV coordinates: a human enters them after confirming the DUA allows it (rule 8); the agent does not submit data there. Dosage map: ClinGen gene/region curation TSVs `<check current: https://ftp.clinicalgenome.org/>`. Population frequency: gnomAD SVs v4.1.0 (genomes, GRCh38) and gnomAD CNVs v4.1.0 (exome CNVs) `<check current: https://gnomad.broadinstitute.org/downloads>`; SVAFotate bundle aggregates gnomAD-SV, CCDG, 1000G, TOPMed (`-f 0.8` = 80% reciprocal overlap; `Max_AF`, `Max_PopMax_AF` added to INFO). Validate: `cut -f1-4 annotsv/${SAMPLE}.annotsv.tsv | head` and the `AnnotSV_ranking_score` column is populated.

ACMG/ClinGen 2020 CNV scoring (Riggs et al., Genet Med 22:245) — five sections summed to a point score: (1) genomic content (protein-coding genes/known elements present?); (2) overlap with established haploinsufficient/triplosensitive genes or regions (HI/TS score 3 → +1.0; benign region → −1.0); (3) gene count (loss: 0–24 genes 0, 25–34 +0.45, ≥35 +0.9; gain: 0–34 0, 35–49 +0.45, ≥50 +0.9); (4) literature/database case evidence (de novo in affected cases up to +0.90; common in gnomAD-SV/DGV controls down to −1.0); (5) inheritance (de novo +0.15–0.45; inherited from unaffected parent −0.30 to −0.45). Totals: ≥0.99 pathogenic; 0.90–0.98 likely pathogenic; −0.89 to 0.89 VUS; −0.90 to −0.98 likely benign; ≤−0.99 benign. Record ClinGen dosage-map date, gnomAD-SV version and AnnotSV version with every score.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

### 10. STR / MEI and mosaic CNV (brief)

- STR: ExpansionHunter with the Clinical-Genomics extended variant catalog → `.claude/skills/genomics-annotation/rare-disease.md` (nf-core/raredisease runs it).
- MEI: `java -jar MELT.jar Single -bamfile ${BAM} -h ${REF} -n ${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed -t mei_list.txt -w melt_${SAMPLE}` (Alu/L1/SVA; `mei_list.txt` lists the `me_refs/Hg38/*.zip`); xTea (`xtea -i bams.txt -r ${REF} -g gencode.gff3 ...`) supports Illumina and long reads; Manta/dysgu report MEIs only as partially assembled INS.
- Mosaic CNV: integer-CN HMMs (gCNV, ExomeDepth) miss events below ~30–40% cell fraction. Signs: segment log2 between 0 and −0.4 (loss) or 0 and +0.3 (gain) with BAF split <0.5/0.5. Use `cnvkit.py call -m none -v ${SAMPLE}.vcf` to keep continuous log2 + BAF, or CNVpytor `-baf`; treat cell fraction $f = 2(1-2^{\log_2 R})$ for a loss as an estimate; confirm with a second tissue or orthogonal assay.

### 11. Visualisation

```bash
set -euo pipefail
SAMPLE=sample
BAM=sample.bam
REF=ref/GRCh38.fa
cnvkit.py scatter "cnvkit/${SAMPLE}.cnr" -s "cnvkit/${SAMPLE}.call.cns" -c chr17:7400000-7800000 -o "${SAMPLE}.TP53.scatter.png"
samplot plot -n "${SAMPLE}" -b "${BAM}" -r "${REF}" -c chr1 -s 1000000 -e 1050000 -t DEL -o "${SAMPLE}.chr1_del.png"
printf 'new\ngenome hg38\nload %s\nload cnvkit/%s.seg\nsnapshotDirectory snapshots\ngoto chr17:7400000-7800000\nsnapshot %s_TP53.png\nexit\n' "${BAM}" "${SAMPLE}" "${SAMPLE}" > igv_batch.txt
xvfb-run -a igv.sh -b igv_batch.txt
```

R: use `CopyNumberPlots` (loads CNVkit/DECoN/GATK segments; `plotCopyNumberCalls()` on a `karyoploteR` ideogram), not the archived `copynumber` package. Validate: files exist and are non-empty (`find . -name "${SAMPLE}*.png" -size +10k`).

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| CNVkit `metrics` MAD of bin log2 (`.cnr`) | WES < 0.30; WGS 1 kb bins < 0.25 | MAD tracks noise independent of true CNVs; > 0.5 means PoN mismatch, GC/FFPE artefacts or low coverage; FFPE tolerates up to ~0.4 with wider segments |
| CNVkit `metrics` segment count | WES germline < 100; tumour 50–500 | Hypersegmentation (>1,000) signals noise or wrong `--method`; raise `segment -t` threshold or drop sample |
| gCNV events per sample (PostprocessGermlineCNVCalls segments VCF, non-ref) | WES < ~100; flag > 2× cohort median | Outliers reflect sample noise or sex-mismatch to priors; GATK WDL enforces `maximum_number_events_per_sample` |
| ExomeDepth/DECoN reference correlation | ≥ 0.98 (DECoN `--mincorr`) | Below: wrong kit/batch in reference set; expected BF for a real single-exon CNV ≥ 10 |
| DECoN minimum median exon coverage | ≥ 100× (`--mincov`) | Below: exons flagged "failed"; report as untested, not negative |
| PoN composition | ≥ 20 normals (GATK), ≥ 10 (CNVkit), same kit, same library protocol, same sequencer, sex-matched or sex-aware | Batch effects dominate WES depth; a mismatched PoN raises MAD and creates recurrent false CNVs at GC extremes; GATK `--minimum-interval-median-percentile 5` drops low-count intervals |
| Sex concordance | gCNV `contig_ploidy.tsv` chrX/chrY ploidy = metadata sex | Mismatch → relabel or exclude before CNV calls on allosomes |
| Allelic counts for ModelSegments | ≥ 1,000 het sites WES, ≥ 100,000 WGS | Fewer → BAF segmentation unstable; use a denser common-SNP list |
| SV evidence (per sample) | Manta `PR_alt ≥ 3` or `SR_alt ≥ 2`; DELLY `GQ ≥ 15`, `DV+RV ≥ 3`; smoove `DHFFC < 0.7` (DEL), `DHBFC > 1.3` (DUP) | duphold thresholds retain ~98% of GIAB DELs while removing most FPs; relax `SR` for IMPRECISE calls in repeats |
| SV size | ≥ 50 bp | Below 50 bp belongs to indel callers; Manta `candidateSmallIndels` feeds Strelka |
| Caller concordance | require 2 of {Manta, smoove, dysgu/DELLY} for discovery; single-caller OK for targeted confirmation | Union is high-recall/low-precision; intersection misses INS |
| Genome-wide SV count (30× WGS, PASS) | Manta 5,000–10,000; smoove 3,000–8,000 | Far outside → contamination, wrong exclude BED or insert-size failure |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `gatk FilterIntervals` fails: "mappability annotation not found" | `AnnotateIntervals` run without `--mappability-track` but filter asked for it | Add track to `AnnotateIntervals` or drop `--minimum-mappability`/`--maximum-segmental-duplication-content` |
| gCNV calls a whole-chromosome CNV on chrX/chrY | Sex not modelled; priors wrong or missing contig | Fix `ploidy_priors.tsv`; check `contig_ploidy.tsv`; set `--allosomal-contig chrX --allosomal-contig chrY` in postprocess |
| `GermlineCNVCaller` memory/time explode on WGS | Unsharded genome | Shard filtered intervals (5–10k per shard); run shards in parallel |
| Recurrent "deletions" in all samples at the same exons | PoN/reference from another kit or batch; GC extremes | Rebuild PoN per kit/batch; `FilterIntervals` GC 0.1–0.9; check capture BED version |
| CNVkit: every amplicon sample shows depth 0 after dedup | Duplicates marked on amplicon libraries | Re-run CNVkit on BAMs without `MarkDuplicates`; `--method amplicon` |
| CNVkit tumour-only: noisy flat-reference profile | `--normal` given with no files | Build a pooled reference from any normals of the same kit, or accept `--drop-low-coverage` + larger `segment -t` and report limited sensitivity |
| CNVkit `call --filter ci` errors | `.cns` lacks CI columns | Run `cnvkit.py segmetrics --ci` first |
| Manta inversions appear as 4 BND records | Manta ≥1.5 default | `convertInversion.py` before merging/annotation |
| Manta runs for hours on GRCh38 | Thousands of alt/decoy contigs | `--callRegions primary_chroms.bed.gz` (bgzipped + tabix) |
| `bcftools view -i 'INFO/PR>=3'` returns nothing on Manta | `PR`/`SR` are FORMAT | Use `FMT/PR[0:1]`, `FMT/SR[0:1]` |
| DELLY: `delly call` unknown command | DELLY ≥2 renamed subcommands | `delly sr` (short read), `delly lr`, `delly cnv`; exclude file `human.hg38.excl.tsv` from the repo `excludeTemplates/` |
| `lumpyexpress` missing / extractSplitReads errors | LUMPY is legacy | Use `smoove call --genotype --duphold` (wraps lumpy + svtyper + duphold) |
| smoove `paste` complains about different variant counts | Samples genotyped against different site lists | Re-run `smoove genotype` on every sample with the same `merged.sites.vcf.gz` |
| AnnotSV: "SV type not found" on BED | Missing `-svtBEDcol` | Add `-svtBEDcol <col>`; never add it for VCF input |
| SURVIVOR merge drops INS | Type matching on INS with different SVLEN | Lower size requirement (`... 1 1 0 30`) or merge INS separately with Jasmine `--dup_to_ins` |
| ACMG class 5 on a common CNV | Stale gnomAD-SV/DGV in AnnotSV annotations | Update AnnotSV annotation bundle; cross-check `Max_AF` from SVAFotate (> 1% → benign-leaning) |
| `nf-core/copynumber` not found | Pipeline does not exist | `nf-core/sarek --tools cnvkit,manta,tiddit` or `nf-core/raredisease` |
| R `library(copynumber)` fails | Package archived from Bioconductor | `CopyNumberPlots` + `karyoploteR` |

## nf-core

| Pipeline | What it gives | Key flags |
|---|---|---|
| `nf-core/sarek` | Germline/somatic SV + CNV: Manta, Tiddit, CNVkit, ASCAT, Control-FREEC | `--tools cnvkit,manta,tiddit` (add `ascat,controlfreec` for tumour/normal); `--wes --intervals targets.bed` |
| `nf-core/raredisease` | Rare-disease WGS/WES: GATK gCNV (needs a pre-built model, `--skip_tools germlinecnvcaller` otherwise), Manta, Tiddit, CNVnator, ExpansionHunter, SVDB merge + annotation, ranking | `--analysis_type wgs|wes`; `--target_bed` for WES |

Samplesheet headers:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
```

```csv
patient,sample,bam,bai
```

```csv
sample,lane,fastq_1,fastq_2,sex,phenotype,paternal_id,maternal_id,case_id
```

(sarek FASTQ; sarek `--step variant_calling` from BAM; raredisease.) Profiles, `-resume`, configs and SLURM: `.claude/skills/genomics-workflow/SKILL.md`. Neither `nf-core/copynumber` nor a standalone nf-core CNV pipeline exists.

## Validation

- Truth sets: GIAB HG002 SV benchmark (assembly-based, HG002-T2T-Q100 v1.1 `*_stvar.vcf.gz` + `.benchmark.bed`, GRCh37/GRCh38/CHM13) with `truvari bench -s 50 --pick single`; CNV truth from the same HG002 set (DEL/DUP ≥50 bp) or array CNV calls. Mechanics, stratifications and expected precision/recall: `.claude/skills/genomics-benchmark/SKILL.md`.
- Internal positives: known pathogenic CNVs in the cohort (e.g. clinically confirmed exon deletions) must be recovered; compute per-exon sensitivity on spiked-in or previously MLPA-confirmed samples before clinical use.
- Sanity stats: per-sample event counts vs cohort distribution; DEL/DUP ratio ≈ 1.5–2.5 in germline WGS SV sets; no CNVs on contigs excluded by the interval list; sex-chromosome calls consistent with inferred sex.
- Orthogonal confirmation: BAF shift inside DEL/DUP (`bcftools query -f '[%AD]' | ...` over het sites), split-read support in IGV/samplot, qPCR/MLPA/array for clinical reporting.
- Reproducibility: record `gatk --version`, `cnvkit.py version`, `smoove --version`, AnnotSV version + annotation bundle date, gnomAD-SV/ClinGen dosage dates, PoN sample list and hash in `results/versions.txt`.

## Additional resources

- `reference.md` (this dir): full GATK gCNV cohort + case + postprocess commands; GATK somatic CNV PoN → ModelSegments → plots; CNVkit step-by-step and BAF/purity; smoove population workflow; DELLY, Tiddit, dysgu, duphold standalone; ExomeDepth R script. Load when running any of these end to end.
- `.claude/skills/genomics-somatic/SKILL.md`: ASCAT/PURPLE/FACETS purity-ploidy, allele-specific CN, tumour-only SNV context. Load for somatic absolute copy number.
- `.claude/skills/genomics-longread/SKILL.md`: Sniffles2, cuteSV, pbsv, Severus, DELLY `lr`. Load for ONT/PacBio SV.
- `.claude/skills/genomics-benchmark/SKILL.md`: Truvari/hap.py, GIAB downloads, stratified metrics. Load when proving sensitivity.
- `.claude/skills/genomics-annotation/SKILL.md`: VEP for SV consequences, HPO/Exomiser ranking, ACMG sequence-variant rules. Load when combining CNVs with SNV/indel interpretation.
- `.claude/skills/genomics-annotation/rare-disease.md`: ExpansionHunter, trio workflows. Load for STR and family-based CNV filtering.
- `.claude/skills/genomics/reference.md`: GRCh38 analysis-set, mappability/segdup tracks, gnomAD/ClinGen download mechanics. Load before building new resources.
- `.claude/skills/genomics-workflow/SKILL.md`: nf-core execution mechanics. Load when running sarek/raredisease on a cluster.
