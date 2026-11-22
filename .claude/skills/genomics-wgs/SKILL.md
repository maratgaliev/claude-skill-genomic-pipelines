---
name: genomics-wgs
description: Germline short-read WGS/WES SNV+indel calling from FASTQ/BAM/CRAM to filtered VCF, single-sample or joint cohort (bwa-mem2, GATK HaplotypeCaller/GenomicsDB/VQSR, DeepVariant+GLnexus, nf-core/sarek). Use when the user mentions HaplotypeCaller, DeepVariant, gVCF, GenomicsDBImport, joint genotyping, VQSR, hard filters, BQSR, sarek germline, exome capture targets, or asks to "call variants" on human/mammalian short reads. Not for somatic tumour/normal (genomics-somatic), long reads (genomics-longread), CNV/SV (genomics-cnv), or variant annotation/ACMG (genomics-annotation).
argument-hint: "[wgs|wes] [cohort] [path]"
paths:
  - "**/*.g.vcf.gz"
  - "**/*.gvcf.gz"
  - "**/genomicsdb*/**"
  - "**/*.recal"
  - "**/*.tranches"
---

# Germline WGS/WES short-variant calling

Expert in germline SNV/indel discovery from Illumina/MGI short reads: alignment, duplicate marking, optional BQSR, HaplotypeCaller or DeepVariant, joint genotyping (GenomicsDB or GLnexus), and filtering (VQSR / NVScoreVariants / hard filters). Outputs are pinned commands, decision tables, resource estimates, and validation steps; never prose tutorials.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" gatk bwa-mem2 samtools bcftools docker singularity glnexus_cli somalier mosdepth VerifyBamID nextflow || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version. Pins used below (verified 2026-11-22): GATK 4.7.0.0, DeepVariant 1.10.0, GLnexus v1.4.1 (`ghcr.io/dnanexus-rnd/glnexus:v1.4.1`), somalier v0.3.7, nf-core/sarek 3.12.0.

## Inspect first

| Input | Command | Infer |
|---|---|---|
| FASTQ | `peek.sh R1.fastq.gz 8`; `ls -l` | Instrument/flowcell/lane from read name (`@A01234:…:HXXXX:1:`), read length, paired layout (`_R1/_R2`, `_1/_2`), lanes to merge, approximate yield (gz bytes ÷ ~2.5 ≈ bases) |
| BAM/CRAM | `samtools view -H x.bam \| grep -E '^@(SQ\|RG\|PG)'` | Assembly + contig style (`chr1` vs `1`; `_alt` contigs present ⇒ alt-aware ref), `@RG` SM/LB/PL completeness, prior MarkDuplicates/BQSR in `@PG`, CRAM `UR:` reference path |
| gVCF/VCF | `bcftools view -h x.g.vcf.gz \| grep -E '^##(contig\|GATKCommandLine\|DeepVariant\|source)'`; `bcftools query -l` | Caller + version, sample names, reference, `<NON_REF>` ⇒ gVCF |
| Targets BED | `peek.sh targets.bed`; `awk '{s+=$3-$2} END{print s}' targets.bed` | WES vs panel (≈30–60 Mb exome; <5 Mb panel ⇒ hard filters only), contig naming |
| Reference | `cut -f1,2 ref.fa.fai \| head -30`; `ls ref.fa.*` | Alt contigs, decoys (`hs38d1`), bwa-mem2 index present (`.0123`, `.bwt.2bit.64`), `.dict` present |

Then run `ref_sanity.sh ${REF} <bam/vcf/bed...>` once. Ask only if still ambiguous: WGS vs WES (no BED, coverage looks targeted); expected sex/pedigree for ploidy and trio logic; whether BQSR/MarkDuplicates were already applied upstream when `@PG` is missing; cohort size target (decides filtering strategy).

## Decision table

| Situation | Choice | Why |
|---|---|---|
| Single sample or ≤10 samples, Illumina/MGI, GRCh38 | DeepVariant 1.10.0 (`--model_type WGS` or `WES`) on the **dedup BAM** | Highest single-sample accuracy; models were trained without BQSR so skip BQSR for DV |
| Cohort to be joint-genotyped with DeepVariant | DeepVariant gVCF per sample → GLnexus `--config DeepVariantWGS|DeepVariantWES` | GLnexus is the supported DV joint caller; GATK GenotypeGVCFs does not accept DV gVCFs |
| Cohort needing GATK annotations/VQSR, or non-human with GATK resources | HaplotypeCaller `-ERC GVCF` → GenomicsDBImport → GenotypeGVCFs | Standard GATK joint pipeline; incremental cohort growth via `--genomicsdb-update-workspace-path` |
| Illumina DRAGEN parity requested, single sample | HaplotypeCaller `--dragen-mode true` + DRAGMAP + `CalibrateDragstrModel`; QUAL-only hard filter | GATK current single-sample best practice; no BQSR needed |
| GPU node available | Parabricks `pbrun germline`/`pbrun deepvariant`; or Sentieon DNAseq/DNAscope (licence) | Same outputs, 10–50× faster; keep version pins in the record |
| Pipeline-grade reproducibility, many samples | nf-core/sarek 3.12.0 `--tools haplotypecaller,deepvariant` | Containers, MultiQC, joint genotyping built in (see `## nf-core`) |
| Filtering: ≥30 WES or ≥1 WGS (GATK minimum; ≥30 WGS preferred) | VQSR (`reference.md`) | Gaussian mixture needs enough variants per annotation |
| Filtering: single WGS/WES sample, GATK ≥4.6 | NVScoreVariants + FilterVariantTranches | Replaces deprecated CNNScoreVariants |
| Filtering: small panels / few samples / no truth resources | Hard filters (`reference.md`) | Only option that needs no training data |
| WES | Always `-L targets.bed -ip 100` (HC) / `--regions targets.bed` (DV); `--wes` + `--intervals` in sarek | Off-target calls inflate FP and break Ti/Tv expectations |
| Male samples, X/Y | DV `--haploid_contigs chrX,chrY --par_regions_bed GRCh38_PAR.bed`; HC `--ploidy 1` on non-PAR X/Y intervals | Diploid genotyping on haploid regions yields spurious hets |
| Somatic tumour/normal | Stop here → `.claude/skills/genomics-somatic/SKILL.md` | Different models, PoN, contamination handling |
| PacBio HiFi / ONT | Stop here → `.claude/skills/genomics-longread/SKILL.md` | DV PACBIO/ONT models, Clair3, phasing |

BQSR: run it for HaplotypeCaller on non-DRAGEN inputs (GATK best practice); skip for DeepVariant and for `--dragen-mode`. Record the decision.

## Workflow

Reference/known-sites acquisition (Broad hg38 bundle, dbSNP, Mills, HapMap/Omni/1000G, Axiom) lives in `.claude/skills/genomics/reference.md`; load it before step 1. NCBI dbSNP VCFs carry `NC_0000xx` contig names and need `bcftools annotate --rename-chrs` before any GATK step. Raw-read QC/trimming is `.claude/skills/genomics-qc/SKILL.md`; modern PCR-free libraries usually need no trimming.

### 1. Align, mark duplicates, index

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
LIB=lib1
FLOWCELL_LANE=HXXXX.1
R1=${SAMPLE}_R1.fastq.gz
R2=${SAMPLE}_R2.fastq.gz
THREADS=16
bwa-mem2 mem -t ${THREADS} -K 100000000 -Y \
  -R "@RG\tID:${FLOWCELL_LANE}\tSM:${SAMPLE}\tLB:${LIB}\tPL:ILLUMINA\tPU:${FLOWCELL_LANE}" \
  ${REF} ${R1} ${R2} \
  | samtools fixmate -@ 4 -m - - \
  | samtools sort -@ 8 -m 2G -T tmp.${SAMPLE} -o ${SAMPLE}.sorted.bam - \
  && samtools quickcheck ${SAMPLE}.sorted.bam
samtools markdup -@ 8 -s -f ${SAMPLE}.markdup.txt ${SAMPLE}.sorted.bam ${SAMPLE}.dedup.bam
samtools index -@ 4 ${SAMPLE}.dedup.bam
```

Validate: `samtools flagstat -@ 4 ${SAMPLE}.dedup.bam | grep -E 'mapped|duplicates'` and `samtools view -H ${SAMPLE}.dedup.bam | grep -c '^@RG'` (= number of lanes).

Notes: `-K 100000000` makes output deterministic regardless of thread count; `-Y` keeps supplementary alignments soft-clipped (required by DRAGEN-mode and many SV callers). One `@RG` per flowcell-lane, same `SM`; merge lanes with `samtools merge` before markdup. Alternative markdup with metrics in Picard format: `gatk MarkDuplicatesSpark -I ${SAMPLE}.sorted.bam -O ${SAMPLE}.dedup.bam -M ${SAMPLE}.dup_metrics.txt --conf 'spark.executor.cores=8'` (sorts internally; needs ~3× BAM size in `--tmp-dir`). CRAM for storage: `samtools view -@ 8 -C -T ${REF} -o ${SAMPLE}.dedup.cram ${SAMPLE}.dedup.bam && samtools index ${SAMPLE}.dedup.cram` (keep `${REF}` immutable; CRAM is unreadable without it).

### 2. Sex inference and contamination (before calling)

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
SITES=/refs/somalier/sites.hg38.vcf.gz
somalier extract -d somalier/ --sites ${SITES} -f ${REF} ${SAMPLE}.dedup.bam
somalier relate --infer -o somalier/cohort somalier/*.somalier
VerifyBamID --SVDPrefix /refs/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat --Reference ${REF} --BamFile ${SAMPLE}.dedup.bam --NumThread 4 --Output ${SAMPLE}.vbid2
```

Validate: `cut -f2,5,7 somalier/cohort.samples.tsv | head` (sample_id, inferred sex 1=male 2=female, pedigree sex; compare to metadata) and `grep -v '^#' ${SAMPLE}.vbid2.selfSM | cut -f7` (FREEMIX). somalier `relate` also flags unexpected duplicates/relatedness in cohorts. `sites.hg38.vcf.gz` ships with the somalier release page.

### 3a. HaplotypeCaller path (with BQSR)

```bash
set -euo pipefail
MODE=WES
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
DBSNP=/refs/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf
MILLS=/refs/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
KNOWN_INDELS=/refs/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz
TARGETS=targets.bed
PAD=100
INTERVAL_ARGS=()
[ "${MODE}" = WES ] && INTERVAL_ARGS=(-L "${TARGETS}" -ip "${PAD}")
gatk --java-options "-Xmx8g" BaseRecalibrator -R ${REF} -I ${SAMPLE}.dedup.bam \
  --known-sites ${DBSNP} --known-sites ${MILLS} --known-sites ${KNOWN_INDELS} \
  "${INTERVAL_ARGS[@]}" -O ${SAMPLE}.recal.table
gatk --java-options "-Xmx8g" ApplyBQSR -R ${REF} -I ${SAMPLE}.dedup.bam \
  --bqsr-recal-file ${SAMPLE}.recal.table -O ${SAMPLE}.bqsr.bam
gatk --java-options "-Xmx8g" HaplotypeCaller -R ${REF} -I ${SAMPLE}.bqsr.bam \
  -D ${DBSNP} "${INTERVAL_ARGS[@]}" -ERC GVCF \
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
  --native-pair-hmm-threads 4 -O ${SAMPLE}.g.vcf.gz
```

Validate: `bcftools view -H ${SAMPLE}.g.vcf.gz | head -1` (first record present) and `bcftools query -f '%ID\n' -i 'ID!="."' ${SAMPLE}.g.vcf.gz | head -1` (dbSNP IDs attached).

WGS: `MODE=WGS` or, faster, scatter by interval shards (`gatk SplitIntervals --scatter-count 24`, one HC per shard with `-L shard.interval_list`, then `gatk GatherVcfs -I … -O ${SAMPLE}.g.vcf.gz` + `tabix`); generation and gather commands in `reference.md`. ~4 GB heap per shard is enough; 8 GB for whole-genome single run. Male X/Y: run non-PAR X/Y intervals separately with `--ploidy 1` and gather (intervals in `reference.md`).

### 3b. DeepVariant path (dedup BAM, no BQSR)

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
MODEL=WGS
THREADS=16
PAR_BED=/refs/GRCh38/GRCh38_PAR.bed
docker run --rm -v "$(dirname ${REF})":/ref -v "${PWD}":/data google/deepvariant:1.10.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=${MODEL} \
  --ref=/ref/$(basename ${REF}) \
  --reads=/data/${SAMPLE}.dedup.bam \
  --output_vcf=/data/${SAMPLE}.dv.vcf.gz \
  --output_gvcf=/data/${SAMPLE}.dv.g.vcf.gz \
  --num_shards=${THREADS} \
  --haploid_contigs=chrX,chrY --par_regions_bed=/ref/$(basename ${PAR_BED}) \
  --intermediate_results_dir=/data/dv_tmp_${SAMPLE}
```

Validate: `bcftools stats -f PASS ${SAMPLE}.dv.vcf.gz | grep -E '^SN|^TSTV' | head -12`.

DeepVariant models (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA) were trained on BAMs without BQSR; feeding BQSR output shifts base-quality channels away from training distribution — always use the dedup BAM. WES: `--model_type=WES --regions=/data/targets.bed` (restricts to capture targets; DV does not pad, pad the BED yourself with `bedtools slop -b 100` if desired). Omit `--haploid_contigs` for female/unknown-sex samples (it forces hom calls on all of X). CRAM input works if the reference is mounted at the path recorded in the CRAM header or passed as `--ref`. Singularity: `singularity exec -B "$(dirname ${REF})":/ref,"${PWD}":/data docker://google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant …`.

### 4. Joint genotyping

GATK gVCFs:

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
DBSNP=/refs/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf
INTERVALS=calling.interval_list
COHORT=cohort
gatk --java-options "-Xmx16g -Xms16g" GenomicsDBImport \
  --sample-name-map ${COHORT}.sample_map.tsv \
  --genomicsdb-workspace-path ${COHORT}_gdb \
  -L ${INTERVALS} --merge-input-intervals \
  --batch-size 50 --reader-threads 4 --genomicsdb-shared-posixfs-optimizations \
  --tmp-dir tmp_gdb
gatk --java-options "-Xmx16g" GenotypeGVCFs -R ${REF} -V gendb://${COHORT}_gdb \
  -D ${DBSNP} -G StandardAnnotation -G AS_StandardAnnotation \
  --only-output-calls-starting-in-intervals -L ${INTERVALS} -O ${COHORT}.raw.vcf.gz
```

Validate: `bcftools query -l ${COHORT}.raw.vcf.gz | wc -l` equals `wc -l < ${COHORT}.sample_map.tsv`.

`sample_map.tsv` is `sample<TAB>/abs/path.g.vcf.gz`, no header. GenomicsDB heap: ≈ batch-size × 50 MB plus 4 GB; set `-Xmx` below container limit by ≥2 GB (native memory). Workspace path must not pre-exist; scatter the import by interval shards for >200 WGS samples. Add later samples with `--genomicsdb-update-workspace-path ${COHORT}_gdb` (no `-L`). Genotyping a single sample: `GenotypeGVCFs -V ${SAMPLE}.g.vcf.gz` directly.

DeepVariant gVCFs: GLnexus (`--config DeepVariantWGS` or `DeepVariantWES`), full block in `reference.md`:

```bash
set -euo pipefail
COHORT=cohort
CONFIG=DeepVariantWGS
ls *.dv.g.vcf.gz | sed 's#^#/data/#' > ${COHORT}.gvcf_list.txt
docker run --rm -v "${PWD}":/data -w /data ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
  /usr/local/bin/glnexus_cli --config ${CONFIG} --threads 16 --dir /data/GLnexus.DB --list /data/${COHORT}.gvcf_list.txt \
  | bcftools view -Oz -o ${COHORT}.glnexus.vcf.gz - \
  && bcftools index -t ${COHORT}.glnexus.vcf.gz
```

Validate: `bcftools query -l ${COHORT}.glnexus.vcf.gz | wc -l`.

### 5. Filter

| Cohort | Method | Block |
|---|---|---|
| ≥30 WES / ≥1 WGS (GATK stated minimum; VQSR stabilises above ~30 WGS) | VQSR SNP (truth-sensitivity 99.7) then INDEL (`--max-gaussians 4`, 99.7) | `reference.md` § VQSR |
| Single sample, HaplotypeCaller VCF, GATK ≥4.6 | `NVScoreVariants` (1D default; 2D with `--tensor-type read_tensor -I bam`) → `FilterVariantTranches --info-key CNN_1D` | `reference.md` § NVScoreVariants |
| <30 WES, panels, non-model organisms without resources | Hard filters: SNP `QD<2 ‖ QUAL<30 ‖ SOR>3 ‖ FS>60 ‖ MQ<40 ‖ MQRankSum<-12.5 ‖ ReadPosRankSum<-8`; INDEL `QD<2 ‖ QUAL<30 ‖ FS>200 ‖ ReadPosRankSum<-20` | `reference.md` § Hard filters |
| DRAGEN mode | `VariantFiltration --filter-expression "QUAL < 10.4139" --filter-name DRAGENHardQUAL` | — |
| DeepVariant / GLnexus | Already filtered: keep `FILTER=PASS` (`bcftools view -f PASS`); GLnexus emits `GQ`-based soft filters | — |

RankSum annotations (`MQRankSum`, `ReadPosRankSum`, `BaseQRankSum`) need both ref and alt reads, so they are **absent at hom-alt sites**; GATK `VariantFiltration` treats a missing annotation as pass (with a warning) — do not "fix" this by `--missing-values-evaluate-as-failing`, which would drop every hom-alt call. VQSR on WES: drop `-an DP` (capture depth is not informative).

### 6. Normalise and summarise

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
VCF=cohort.filtered.vcf.gz
OUT=cohort.final.vcf.gz
bcftools norm -f ${REF} -m -both --check-ref w -Oz -o ${OUT} ${VCF} 2> ${OUT%.vcf.gz}.norm.log
bcftools index -t ${OUT}
bcftools stats -f PASS -F ${REF} -s - ${OUT} > ${OUT%.vcf.gz}.stats.txt
```

Validate: `grep -E '^TSTV' ${OUT%.vcf.gz}.stats.txt` and `grep -c 'REF_MISMATCH' ${OUT%.vcf.gz}.norm.log` (must be 0; non-zero ⇒ wrong reference; the warnings go to stderr, not the stats file).

### Accelerated equivalents (same outputs, record version)

```bash
# alternatives: run one (Parabricks germline, Parabricks DeepVariant, Sentieon DNAscope); variables as in steps 1–3
pbrun germline --ref ${REF} --in-fq ${R1} ${R2} "@RG\tID:${FLOWCELL_LANE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${LIB}" --knownSites ${DBSNP} --out-bam ${SAMPLE}.bam --out-variants ${SAMPLE}.g.vcf.gz --gvcf --num-gpus 1
pbrun deepvariant --ref ${REF} --in-bam ${SAMPLE}.dedup.bam --out-variants ${SAMPLE}.dv.vcf.gz --mode shortread --num-gpus 1
sentieon driver -t ${THREADS} -r ${REF} -i ${SAMPLE}.dedup.bam --algo DNAscope --model ${SENTIEON_MODEL} --emit_mode gvcf ${SAMPLE}.g.vcf.gz
```

`pbrun deepvariant` exposes `--haploid-contigs` from Parabricks 4.7; check `pbrun deepvariant --help` for `--par-regions-bed` on your version.

## Resource estimates (approximate; 30× human WGS, 16 cores, scale WES ÷ 10–20)

| Step | RAM | Wall time | Disk |
|---|---|---|---|
| `bwa-mem2 index` (once) | ≈70 GB peak for GRCh38 | ~1 h | ~15 GB index |
| bwa-mem2 mem + sort | 20–30 GB (index ~12 GB + sort buffers) | 6–10 h | sorted BAM ≈ 60–80 GB |
| samtools markdup / MarkDuplicatesSpark | 4 GB / 16 GB+ | 1–2 h / 2–4 h | +1× BAM (tmp) |
| BQSR (BaseRecalibrator + ApplyBQSR) | 8 GB | 3–5 h unscattered | +1× BAM |
| HaplotypeCaller GVCF, 24 shards | 4 GB per shard | 1–2 h with 16 parallel shards; 15–25 h single process | gVCF 5–10 GB |
| DeepVariant 1.10 CPU, 16 shards | 16–32 GB | 4–8 h (GPU: <1 h) | gVCF 3–6 GB |
| GenomicsDBImport, 100 WGS, scattered | 16 GB per shard | 2–6 h | ≈ sum of gVCFs |
| VQSR | 16–32 GB | 1–3 h | small |

Scratch rule: reserve ≥ 4× final BAM size during alignment + markdup. Pass `--tmp-dir` on a local disk to every GATK step; network scratch is the common cause of "slow" HaplotypeCaller.

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Mapped reads (`samtools flagstat`) | ≥ 98 % WGS and WES (WARN 95–98 %, FAIL < 95 %) | Lower ⇒ wrong reference/species, adapter read-through, contamination; check `peek.sh` and a BLAST of 100 unmapped reads |
| Duplicate rate (marked, not removed) | < 10 % PCR-free WGS (WARN 10–20 %), ≤ 20 % PCR WGS, < 25 % WES (WARN 25–40 %, FAIL > 40 %) | Higher ⇒ low input/overamplification; high-depth (>60×) WES legitimately higher; NovaSeq X patterned flowcells add optical duplicates |
| Mean coverage (`mosdepth -n --fast-mode`) | ≥ 30× WGS; ≥ 80× WES mean on-target | Clinical WGS often ≥ 40×; trios tolerate 25× |
| Coverage uniformity | ≥ 90 % of genome ≥ 10× (WGS); ≥ 95 % of targets ≥ 20× (WES; WARN 90–95 %) | GC-biased libraries fail this with acceptable mean; report `mosdepth --thresholds` |
| On-target fraction (WES) | ≥ 60 % | Kit-dependent; <40 % ⇒ wrong BED or failed capture |
| Insert size median | 300–500 bp, unimodal | Bimodal ⇒ adapter dimers or mixed libraries; <200 bp ⇒ overlapping mates, lower effective depth |
| Ti/Tv (PASS, `bcftools stats`) | 2.0–2.1 WGS; 2.8–3.3 WES | Lower ⇒ false positives (typical for unfiltered/hard-filtered); recompute on targets only for WES |
| Het/hom-alt ratio | ≈ 1.5–2.0 (WGS, European ancestry); varies by ancestry | <1.2 ⇒ consanguinity/LOH; >2.5 ⇒ contamination or caller issue |
| Contamination FREEMIX (VerifyBamID2) | < 0.01 (WARN 0.01–0.03, FAIL > 0.03 germline) | Somatic normals often enforced < 0.02; low depth (<5×) or wrong panel build inflates it; > 0.03 exclude from cohort filtering |
| Sex (somalier) | matches metadata | Mismatch ⇒ sample swap; treat as blocker before joint calling |
| Variant count (PASS, WGS GRCh38) | 4.0–5.0 M SNVs, 0.8–1.0 M indels; WES ≈ 20–40 k | Outside ⇒ wrong intervals, wrong ref, or filtering failure |
| dbSNP concordance (`bcftools stats` with `-D`/ID) | ≥ 97 % known WGS | Lower ⇒ FP burden or non-human-like ancestry representation |

Canonical PASS/WARN/FAIL values live in `.claude/skills/genomics-qc/SKILL.md`; this table restates them for the calling context.

Raw-read QC (FastQC/fastp, adapter content, per-base quality) is `.claude/skills/genomics-qc/SKILL.md`.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `Input files reference and reads have incompatible contigs` | `chr` vs no-`chr`, or different assembly | `ref_sanity.sh`; realign or `bcftools annotate --rename-chrs` for VCF resources (NCBI dbSNP `NC_000001.11` style) |
| HaplotypeCaller 10× slower than expected | Network `--tmp-dir`, `-Xmx` > container memory causing swap, or no scatter | Local tmp, `-Xmx` ≤ RAM − 2 GB, `SplitIntervals` scatter |
| GenomicsDBImport "workspace already exists" | Re-run without cleaning | Delete workspace or use `--genomicsdb-update-workspace-path` |
| GenomicsDBImport OOM / too many open files | `--batch-size 0` with hundreds of samples | `--batch-size 50 --reader-threads 4`, raise `ulimit -n` |
| VQSR `No data found` / `Bad input: Values for … annotation not detected` | Too few variants (small panel, single WES) or annotation missing | Switch to hard filters or NVScoreVariants; ensure `-G StandardAnnotation` at calling time |
| VQSR INDEL fails to converge | Default 8 Gaussians on sparse indel set | `--max-gaussians 4` (already in `reference.md`) |
| Many hom-alt variants filtered by hard filters | `--missing-values-evaluate-as-failing` enabled | Remove the flag; RankSum annotations are undefined at hom sites |
| DeepVariant spurious hets on chrX in males | No `--haploid_contigs` | Add `--haploid_contigs=chrX,chrY --par_regions_bed=…` |
| DeepVariant accuracy below published | BQSR BAM used, wrong `--model_type`, or `--regions` omitted for WES | Rerun on dedup BAM with correct model |
| Ti/Tv WES ≈ 2.1 | Stats computed genome-wide including off-target | `bcftools stats -T targets.bed` or call with `-L` |
| Ti/Tv WGS < 1.9 after filtering | FP burden, often from contamination or low depth | Check FREEMIX, coverage; tighten tranche (99.5) |
| `MarkDuplicatesSpark` crashes on sorted-by-name input or ARM Mac | Spark local mode limits | Use `samtools markdup` pipeline (step 1) |
| CRAM unreadable later | Reference moved/changed | Store `${REF}` md5 with the CRAM; set `REF_PATH`/`REF_CACHE` |
| GLnexus `unsupported gVCF` | Non-DeepVariant gVCF or truncated file | Use `--config gatk` for GATK gVCFs only if required; prefer GenomicsDB for GATK inputs |
| sarek `--wes` run calls off-target / fails | `--intervals` missing | `--wes --intervals targets.bed` always together |

## nf-core

`nf-core/sarek` 3.12.0 (exists: https://nf-co.re/pipelines.json). Germline WGS, joint genotyping, both callers:

```bash
nextflow run nf-core/sarek -r 3.12.0 -profile docker \
  --input samplesheet.csv --outdir results \
  --genome GATK.GRCh38 \
  --tools haplotypecaller,deepvariant \
  --joint_germline
```

WES: add `--wes --intervals targets.bed` (both required together; `--wes` alone is rejected). Skip BQSR: `--skip_tools baserecalibrator`. Use `--genome GATK.GRCh38` consistently (iGenomes key; `GRCh38` alone resolves to a different, non-GATK bundle). Samplesheet:

```csv
patient,sample,lane,fastq_1,fastq_2
P1,S1,L001,/data/S1_L001_R1.fastq.gz,/data/S1_L001_R2.fastq.gz
P1,S1,L002,/data/S1_L002_R1.fastq.gz,/data/S1_L002_R2.fastq.gz
```

Add `sex` (`XX`/`XY`) column for ploidy-aware steps. Start from BAM/CRAM with header `patient,sample,bam,bai` (or `cram,crai`) plus `--step markduplicates|prepare_recalibration|variant_calling`. Profiles, `-resume`, SLURM, config overrides: `.claude/skills/genomics-workflow/SKILL.md`. Rare-disease trios with annotation and ranking: `nf-core/raredisease` (exists).

## Validation

- Truth sets: run a GIAB sample (HG001/NA12878, HG002–HG007; IDs and URLs in `reference.md`) through the identical pipeline and benchmark with `hap.py`/`vcfdist` per `.claude/skills/genomics-benchmark/SKILL.md`. Expected 30× WGS DeepVariant 1.10 or GATK+VQSR: SNV F1 ≥ 0.995, indel F1 ≥ 0.99 in GIAB v4.2.1 high-confidence regions (approximate; stratify by difficult regions).
- Without truth: `bcftools stats -f PASS -s -` must show Ti/Tv per table above, per-sample variant counts within 10 % of cohort median, singletons < 5 % of variants per sample, `nMissing` low and even across samples; `bcftools +smpl-stats`/`somalier relate` confirms sex and expected relatedness.
- Mendelian check for trios: `bcftools +mendelian2 -p ped.txt` error rate < 1 % PASS sites (notes in `reference.md`).
- Reference sanity after `bcftools norm --check-ref w`: zero `REF_MISMATCH` lines in the captured stderr (`${OUT%.vcf.gz}.norm.log`).

## Additional resources

- `.claude/skills/genomics-wgs/reference.md` — full VQSR SNP+INDEL blocks, hard-filter blocks, NVScoreVariants + FilterVariantTranches, GLnexus, `SplitIntervals` scatter/gather and haploid X/Y intervals, trio/pedigree (`CalculateGenotypePosteriors`, `VariantAnnotator -A PossibleDeNovo`), GIAB sample IDs. Load when executing step 3–5 or a trio.
- `.claude/skills/genomics/reference.md` — reference FASTA, Broad hg38 bundle, dbSNP/Mills/HapMap/Omni/1000G/Axiom downloads, contig renaming. Load before step 1.
- `.claude/skills/genomics-qc/SKILL.md` — raw-read QC/trimming and MultiQC aggregation.
- `.claude/skills/genomics-benchmark/SKILL.md` — hap.py/vcfdist/stratifications against GIAB truth.
- `.claude/skills/genomics-annotation/SKILL.md` — VEP/ClinVar/gnomAD annotation and variant interpretation of the final VCF.
- `.claude/skills/genomics-somatic/SKILL.md` — tumour/normal and tumour-only SNV/indel calling.
- `.claude/skills/genomics-longread/SKILL.md` — PacBio HiFi/ONT germline calling.
- `.claude/skills/genomics-cnv/SKILL.md` — CNV/SV from the same BAM/CRAM.
- `.claude/skills/genomics-workflow/SKILL.md` — Nextflow/Snakemake execution, profiles, HPC.
