---
name: genomics-somatic
description: Somatic SNV/indel calling from tumour-normal or tumour-only short-read data (Mutect2 full chain, Strelka2, VarDict, DeepSomatic, ensemble), tumour purity/ploidy and allele-specific CN (ASCAT/FACETS/Sequenza/PURPLE), MSI/TMB/HRD/signatures/clonality, and somatic interpretation (AMP/ASCO/CAP tiers, OncoKB, PCGR). Use when inputs are tumour BAM/CRAM pairs, `--panel-of-normals`, `f1r2.tar.gz`, `FilterMutectCalls`, `candidateSmallIndels`, UMI/ctDNA consensus, or the user says "somatic", "tumour-only", "purity", "TMB", "MSI", "actionable". Not for germline SNV/indel (`genomics-wgs`), somatic CNV/SV segmentation mechanics (`genomics-cnv`), long-read somatic callers (`genomics-longread`), or germline ACMG (`genomics-annotation`).
argument-hint: "[call|purity|interpret|umi] [tumour.bam normal.bam]"
paths:
  - "**/*f1r2*.tar.gz"
  - "**/*pon*.vcf.gz"
  - "**/*candidateSmallIndels*.vcf.gz"
  - "**/*.seqz.gz"
  - "**/*contamination*.table"
---

# Somatic variant calling and interpretation

Expert in somatic SNV/indel discovery and tumour genome characterisation from short reads. Outputs are runnable commands, decision tables, QC thresholds, and validation steps; never read data files into context.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" gatk samtools bcftools configureStrelkaSomaticWorkflow.py configManta.py VarDict run_deepsomatic somaticseq alleleCounter snp-pileup sequenza-utils msisensor-pro fgbio somalier pcgr cpsr Rscript python3 nextflow || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| Read | Command | Infer |
|---|---|---|
| BAM/CRAM `@RG SM` | `samtools view -H T.bam \| grep '^@RG'` | Tumour/normal sample names for `-normal`, `--tumor-name`; multiple SM in one file = merge error |
| `@SQ` names/lengths | `ref_sanity.sh ref.fa T.bam N.bam gnomad.vcf.gz pon.vcf.gz` | GRCh38 vs GRCh37, `chr` prefix; `af-only-gnomad.hg38` requires `chr` |
| `@PG` lines | `samtools view -H T.bam \| grep '^@PG'` | Aligner, MarkDuplicates/BQSR done? UMI tools (fgbio) → consensus BAM already |
| Read tags | `samtools view T.bam \| head -2 \| tr '\t' '\n' \| grep -E '^(RX|MI|cD|cM)'` | `RX` raw UMI, `MI`/`cD` consensus → UMI workflow |
| Depth | `samtools coverage T.bam \| head` or mosdepth summary | WES vs WGS vs panel; whether VAF floor is credible |
| Existing VCF header | `bcftools view -h x.vcf.gz \| grep -E '^##(source|GATKCommandLine|tumor_sample|normal_sample)'` | Caller, filtering already applied, sample roles |
| Relatedness | somalier (see `.claude/skills/genomics-qc/SKILL.md`) | Confirm T/N from same individual before any calling |

Ask only if still ambiguous: FFPE vs fresh-frozen; panel target BED and its callable size (TMB); whether a cohort PoN exists; tumour site (for OncoKB/PCGR tiering); consent/approval for OncoKB API use.

## Decision table

| Situation | Choice | Why |
|---|---|---|
| T/N pair, WES/WGS, default | Mutect2 full chain (below) + Strelka2; report intersection or SomaticSeq ensemble | Highest precision; orthogonal error models |
| T/N pair, GPU available | DeepSomatic 1.10.0 `--model_type WGS|WES|FFPE_WGS|FFPE_WES` | Best SEQC2 F1 among single callers; needs ≥30x N |
| Tumour-only | Mutect2 without `-normal` + PoN + gnomAD AF filter; or DeepSomatic `*_TUMOR_ONLY`; PCGR `--tumor_only` | Germline leakage unavoidable: report as "somatic-candidate", not somatic |
| FFPE | Mutect2 with `--f1r2-tar-gz` + `LearnReadOrientationModel`; DeepSomatic FFPE model | C>T/G>A deamination artefacts |
| Panel / ctDNA with UMIs | fgbio consensus → VarDict `-f 0.001`–`0.005` or Mutect2 `--min-allele-fraction`; duplex if available | Raw reads cannot support VAF <1% |
| Low purity (<20%) | Lower VAF threshold only with depth ≥500x; otherwise report limited sensitivity | LOD scales with purity×depth |
| Allele-specific CN + purity/ploidy | ASCAT (default), FACETS (WES), Sequenza (WGS/WES), PURPLE (WGS; via oncoanalyser) | ASCAT needs het SNP BAF; PURPLE integrates SV/SNV |
| Multi-caller ensemble | SomaticSeq `paired` (≥2-of-N consensus or trained classifier) | bcbio-style ≥2-of-3 reduces caller-specific FPs |
| Licensed/commercial speed | Sentieon TNscope; NVIDIA Parabricks `mutectcaller`/`deepsomatic` | Same models, 10–50× faster |
| Full WGS tumour characterisation (CN, SV, drivers, HRD, signatures) | nf-core/oncoanalyser (hmftools WiGiTS) | One pipeline, PURPLE/LINX/CHORD/SIGS |

## Workflow

Full runnable scripts: `reference.md` in this directory. Canonical happy path below.

### 1. Sample-swap and reference checks

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" ${REF} ${TBAM} ${NBAM} af-only-gnomad.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz
somalier extract -d somalier -s sites.hg38.vcf.gz -f ${REF} ${TBAM} ${NBAM}
somalier relate -o tn somalier/*.somalier
```

Validate: `awk -F'\t' 'NR>1 && $3<0.9' tn.pairs.tsv` prints nothing (T/N relatedness ≈1; <0.9 = swap or contamination; LOH in tumour can lower it slightly, not below ~0.8).

### 2. Mutect2 best practice (tumour-normal)

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
NSM=$(samtools view -H ${NBAM} | awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){sub("SM:","",$i); print $i; exit}}')
GNOMAD=af-only-gnomad.hg38.vcf.gz
PON=1000g_pon.hg38.vcf.gz
COMMON=small_exac_common_3.hg38.vcf.gz
INTERVALS=wgs_calling_regions.hg38.interval_list
OUT=somatic
gatk Mutect2 -R ${REF} -I ${TBAM} -I ${NBAM} -normal ${NSM} \
  --germline-resource ${GNOMAD} --panel-of-normals ${PON} \
  --f1r2-tar-gz ${OUT}.f1r2.tar.gz -L ${INTERVALS} \
  --native-pair-hmm-threads 4 -O ${OUT}.unfiltered.vcf.gz
gatk LearnReadOrientationModel -I ${OUT}.f1r2.tar.gz -O ${OUT}.read-orientation-model.tar.gz
gatk GetPileupSummaries -I ${TBAM} -V ${COMMON} -L ${COMMON} -O ${OUT}.tumour.pileups.table
gatk GetPileupSummaries -I ${NBAM} -V ${COMMON} -L ${COMMON} -O ${OUT}.normal.pileups.table
gatk CalculateContamination -I ${OUT}.tumour.pileups.table -matched ${OUT}.normal.pileups.table \
  --tumor-segmentation ${OUT}.segments.table -O ${OUT}.contamination.table
gatk FilterMutectCalls -R ${REF} -V ${OUT}.unfiltered.vcf.gz \
  --contamination-table ${OUT}.contamination.table --tumor-segmentation ${OUT}.segments.table \
  --ob-priors ${OUT}.read-orientation-model.tar.gz -O ${OUT}.filtered.vcf.gz
bcftools view -f PASS -Oz -o ${OUT}.pass.vcf.gz ${OUT}.filtered.vcf.gz && bcftools index -t ${OUT}.pass.vcf.gz
```

Validate: `awk 'NR==2{print "contamination",$2}' ${OUT}.contamination.table` (<0.02–0.03) and `bcftools +counts ${OUT}.pass.vcf.gz`; `${OUT}.filtered.vcf.gz.filteringStats.tsv` lists FP/FN estimates per filter.

Notes: `-L` same file as `-V` in `GetPileupSummaries` (required); `--f1r2-tar-gz` + `--ob-priors` replaces legacy `FilterByOrientationBias`; scatter `Mutect2` by interval with `SplitIntervals` then merge `f1r2` tarballs (`LearnReadOrientationModel -I a.tar.gz -I b.tar.gz`), VCFs (`MergeVcfs`) and stats (`MergeMutectStats`) before filtering. Optional `gatk FilterAlignmentArtifacts -R ${REF} -V ${OUT}.filtered.vcf.gz -I ${TBAM} --bwa-mem-index-image ref.img -O ${OUT}.realigned.vcf.gz` (needs `BwaMemIndexImageCreator` image; costly, use for panels/hotspots).

Tumour-only: drop `-I ${NBAM} -normal`, keep `--germline-resource --panel-of-normals`, filter the same way; expect residual germline (esp. rare, non-gnomAD alleles, populations under-represented in gnomAD); never call a tumour-only variant "somatic" without PoN + gnomAD AF >0.001 exclusion + VAF/CN consistency; `FilterMutectCalls` loses the `normal_artifact` filter.

### 3. Panel of normals

Per normal: `gatk Mutect2 -R ${REF} -I n1.bam --max-mnp-distance 0 -O n1.vcf.gz` → `GenomicsDBImport` → `CreateSomaticPanelOfNormals --germline-resource ${GNOMAD}`. ≥40 normals from the same assay/aligner; never include the matched normal of a case in its own PoN. Script: `reference.md` §PoN.

Validate: `bcftools view -H pon.vcf.gz | wc -l` in expected range (WES 10⁴–10⁵; WGS 10⁶); `bcftools view -h pon.vcf.gz | grep -c '^##source=CreateSomaticPanelOfNormals'` = 1.

### 4. Strelka2 (+ Manta candidate indels)

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
THREADS=16
EXOME_FLAG=""
configManta.py --normalBam ${NBAM} --tumorBam ${TBAM} --referenceFasta ${REF} ${EXOME_FLAG} --runDir manta
manta/runWorkflow.py -m local -j ${THREADS}
configureStrelkaSomaticWorkflow.py --normalBam ${NBAM} --tumorBam ${TBAM} --referenceFasta ${REF} \
  --indelCandidates manta/results/variants/candidateSmallIndels.vcf.gz ${EXOME_FLAG} --runDir strelka
strelka/runWorkflow.py -m local -j ${THREADS}
```

Validate: `bcftools view -f PASS -H strelka/results/variants/somatic.snvs.vcf.gz | wc -l`; set `EXOME_FLAG=--exome` for WES/panel (disables depth filter); add `--callRegions targets.bed.gz` (bgzip+tabix) for panels. Strelka2 writes no VAF field: compute tumour VAF from `AU/CU/GU/TU` tier-1 counts (SNV) or `TAR/TIR` (indel).

### 5. VarDict (low VAF / panels)

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
TSM=tumour
NSM=normal
BED=targets.bed
AF=0.01
VarDict -G ${REF} -f ${AF} -N ${TSM} -b "${TBAM}|${NBAM}" -c 1 -S 2 -E 3 -g 4 -th 8 ${BED} \
  | testsomatic.R | var2vcf_paired.pl -N "${TSM}|${NSM}" -f ${AF} > vardict.vcf
bgzip -f vardict.vcf && tabix -p vcf vardict.vcf.gz
```

Validate: `bcftools view -i 'INFO/STATUS="StrongSomatic"' -H vardict.vcf.gz | wc -l`. VarDict `-c/-S/-E/-g` are BED column indices (chrom/start/end/name), not legacy flags. Keep `StrongSomatic`/`LikelySomatic`; drop `Germline`, `LOH`, `Deletion`.

### 6. DeepSomatic 1.10.0

```bash
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
MODEL=WGS
docker run --rm -v "$PWD":"$PWD" -w "$PWD" google/deepsomatic:1.10.0 run_deepsomatic \
  --model_type=${MODEL} --ref=${REF} --reads_tumor=${TBAM} --reads_normal=${NBAM} \
  --sample_name_tumor=tumour --sample_name_normal=normal \
  --output_vcf=deepsomatic.vcf.gz --num_shards=16 --logging_dir=ds_logs
```

Validate: `bcftools view -f PASS -H deepsomatic.vcf.gz | wc -l`. Model types (v1.10.0): `WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES, WGS_TUMOR_ONLY, PACBIO_TUMOR_ONLY, ONT_TUMOR_ONLY`; tumour-only drops `--reads_normal` and uses `--use_default_pon_filtering=true`. Long-read modes: see `.claude/skills/genomics-longread/SKILL.md`. GPU image tag `1.10.0-gpu`.

### 7. Commercial / GPU one-liners

```bash
# alternatives: run one (Sentieon TNscope or Parabricks Mutect2)
REF=Homo_sapiens_assembly38.fasta
sentieon driver -t 32 -r ${REF} -i tumour.bam -i normal.bam --algo TNscope --tumor_sample tumour --normal_sample normal --dbsnp dbsnp.vcf.gz --pon pon.vcf.gz tnscope.vcf.gz
pbrun mutectcaller --ref ${REF} --in-tumor-bam tumour.bam --tumor-name tumour --in-normal-bam normal.bam --normal-name normal --out-vcf pb_mutect2.vcf
```

Validate: sample names must equal `@RG SM` tags; `bcftools query -l` on output shows both.

### 8. Ensemble

```bash
set -euo pipefail
somaticseq --output-directory somaticseq --genome-reference ${REF} --inclusion-region ${INTERVALS_BED} --threads 16 \
paired --tumor-bam-file ${TBAM} --normal-bam-file ${NBAM} \
  --mutect2-vcf somatic.filtered.vcf.gz --vardict-vcf vardict.vcf.gz \
  --strelka-snv strelka/results/variants/somatic.snvs.vcf.gz --strelka-indel strelka/results/variants/somatic.indels.vcf.gz
```

Validate: `somaticseq/Consensus.sSNV.vcf` + `Consensus.sINDEL.vcf` exist; default consensus = majority vote (`NumCallers` ≥ ceil(N/2)); bcbio-style rule is ≥2 of 3 callers. Train a classifier (`--classifier-model`) only with SEQC2 truth data from the same assay.

### 9. Purity, ploidy, allele-specific CN

| Tool | Input prep | Key call | Output used downstream |
|---|---|---|---|
| ASCAT ≥3.1 | `ascat.prepareHTS()` runs alleleCounter on T/N BAM at 1000G loci (`loci.prefix`, `alleles.prefix`, `genomeVersion="hg38"`) | `ascat.aspcf(penalty=70)` → `ascat.runAscat(gamma=1)` | `aberrantcellfraction`, `ploidy`, `segments` (nMajor/nMinor) |
| FACETS | `snp-pileup -g -q15 -Q20 -P100 -r25,0 dbsnp.vcf.gz out.csv.gz N.bam T.bam` | `preProcSample` → `procSample(cval=150 WES / 300 WGS)` → `emcncf` | `purity`, `ploidy`, `cncf` (tcn.em/lcn.em) |
| Sequenza | `sequenza-utils gc_wiggle -w 50 -f ref.fa -o ref.gc50.wig.gz`; `bam2seqz -n N.bam -t T.bam --fasta ref.fa -gc ref.gc50.wig.gz -o out.seqz.gz`; `seqz_binning -w 50` | `sequenza.extract` → `sequenza.fit` → `sequenza.results` | `*_confints_CP.txt` (cellularity, ploidy), `*_segments.txt` |
| PURPLE (hmftools) | AMBER (BAF) + COBALT (depth ratio) | `java -jar purple.jar -reference N -tumor T -amber_dir amber -cobalt_dir cobalt -somatic_vcf pass.vcf.gz -ref_genome ref.fa -ref_genome_version 38 -gc_profile GC_profile.1000bp.38.cnp -ensembl_data_dir ensembl -output_dir purple` | `*.purple.purity.tsv`, `*.purple.cnv.somatic.tsv`, MSI/TMB/TML in `purity.tsv` |

R scripts: `reference.md` §ASCAT, §FACETS. Validate: purity estimates from two tools agree within ±0.1 and ploidy within ±0.3; disagreement = whole-genome-doubling ambiguity → inspect BAF/logR plots, consider `min_ploidy/max_ploidy` bounds. Segmentation/CNV-caller mechanics live in `.claude/skills/genomics-cnv/SKILL.md`; oncoanalyser runs AMBER/COBALT/PURPLE for you.

### 10. Biomarkers

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
msisensor-pro scan -d ${REF} -o ref.microsat.list
msisensor-pro msi -d ref.microsat.list -n normal.bam -t tumour.bam -o msi/sample -c 20 -b 8
msisensor-pro pro -d ref.microsat.list -t tumour.bam -o msi_to/sample
```

Validate: `head -2 msi/sample` → `%` unstable sites; MSI-H ≥20% (MSIsensor-pro default, also PCGR `--input_msi_threshold`), MSS <10% (10–20% = indeterminate, confirm orthogonally); `-c 20` for WES, 15 for WGS; `pro` (tumour-only) needs a `baseline` built from ≥20 assay-matched normals. MANTIS alternative: `python mantis.py --bedfile loci.bed --genome ${REF} -n normal.bam -t tumour.bam -o mantis.txt` (MSI-H if step-wise difference >0.4).

| Biomarker | Definition / tool | Default threshold |
|---|---|---|
| TMB | non-synonymous (missense+nonsense+frameshift+in-frame+splice) somatic PASS SNV+indel count ÷ callable coding Mb (bases ≥ TMB depth floor inside the target BED, measured with mosdepth, not nominal panel size); exclude germline, VAF <5%, hotspot-only panel bias | ≥10 mut/Mb = TMB-H (FDA/OncoKB; PCGR `--tmb_high_threshold 10`, `--effective_target_size_mb`) |
| HRD scar | scarHRD (`scar_score(seg, reference="grch38", seqz=TRUE)` = HRD-LOH+LST+TAI), CHORD (hmftools, SNV/indel/SV signatures), HRDetect (`signature.tools.lib`); PCGR `--estimate_hrd` | HRD score ≥42 (Myriad, ovarian only); CHORD P(HRD) ≥0.5; HRDetect ≥0.7 |
| Mutational signatures | SigProfilerAssignment `cosmic_fit` (default COSMIC v3.6; pin `--cosmic_version 3.4` to match published work) or MutationalPatterns `fit_to_signatures_strict` with an external COSMIC v3.4 matrix (`get_known_signatures` ships ≤v3.2) | ≥50 SNVs per sample; report cosine similarity of reconstruction ≥0.9 |
| Clonality | PyClone-VI: TSV `mutation_id sample_id ref_counts alt_counts major_cn minor_cn normal_cn tumour_content` from CN + purity above; `pyclone-vi fit -i in.tsv -o out.h5 -c 40 -d beta-binomial -r 10 && pyclone-vi write-results-file -i out.h5 -o clusters.tsv` | Clonal cluster CCF ≈1; report only with purity ≥0.3 and ≥100 SNVs |

### 11. Interpretation

| AMP/ASCO/CAP 2017 tier | Meaning | Evidence |
|---|---|---|
| I | Strong clinical significance | FDA-approved therapy or professional-guideline biomarker for this tumour type (OncoKB Level 1/2, R1; CIViC A) |
| II | Potential clinical significance | FDA-approved therapy in a different tumour type, or investigational/clinical-trial evidence (OncoKB 3A/3B/4, R2; CIViC B/C) |
| III | Unknown clinical significance | Cancer gene, not in population databases, no therapeutic evidence |
| IV | Benign / likely benign | gnomAD AF >1% or known benign |

```bash
set -euo pipefail
pcgr --input_vcf somatic.pass.vcf.gz --input_cna cna_segments.tsv --vep_dir ${VEP_DIR} --refdata_dir ${PCGR_DATA} \
  --output_dir pcgr --genome_assembly grch38 --sample_id ${SAMPLE} --assay WES --tumor_site 15 \
  --tumor_purity 0.6 --tumor_ploidy 2.1 --tumor_dp_tag TDP --tumor_af_tag TAF --control_dp_tag NDP --control_af_tag NAF \
  --estimate_msi --estimate_tmb --estimate_signatures --estimate_hrd --input_msi msi/sample --input_msi_caller msisensor-pro
```

Validate: `ls pcgr/${SAMPLE}.pcgr.grch38.html pcgr/${SAMPLE}.pcgr.grch38.snv_indel.tsv.gz`; PCGR 2.3.2 (2026-08) needs its matching refdata bundle and VEP cache (`<check current: https://sigven.github.io/pcgr/articles/installation.html>`). Germline side of the pair: `cpsr --input_vcf germline.vcf.gz --vep_dir ${VEP_DIR} --refdata_dir ${PCGR_DATA} --output_dir cpsr --genome_assembly grch38 --sample_id ${SAMPLE} --panel_id 0`; ACMG mechanics in `.claude/skills/genomics-annotation/SKILL.md`.

Actionability caveats: OncoKB levels are tumour-type specific (pass `--oncokb_api_token` + `--oncokb_oncotree_code`; commercial/clinical use needs an MSK licence); CGI biomarker set is frozen; CIViC evidence levels A–E carry curator trust ratings, read them; tiers depend on tumour site, so `--tumor_site 0` demotes every site-specific biomarker to Tier II/III; a Tier I variant at VAF 2% in a 30% purity tumour is subclonal, say so; copy-number and fusion evidence need their own calls (see `genomics-cnv`, `nf-core/rnafusion`).

### 12. UMI / liquid biopsy

fgbio chain (full script `reference.md` §UMI): `FastqToBam --read-structures 8M+T +T` → `bwa mem -K 150000000 -Y` + `ZipperBams` → `samtools sort --template-coordinate` → `GroupReadsByUmi --strategy Adjacency --edits 1` (`Paired` strategy for duplex) → `CallMolecularConsensusReads --min-reads 1 --min-input-base-quality 20` (`CallDuplexConsensusReads` for duplex) → realign with `--tags-to-reverse Consensus --tags-to-revcomp Consensus` → `FilterConsensusReads --min-reads 3 --min-base-quality 45 --max-base-error-rate 0.2`. nf-core/fastquorum 2.0.0 implements exactly this. Call on the consensus BAM with VarDict `-f 0.001` or Mutect2 `--min-allele-fraction`; never MarkDuplicates after consensus.

| Depth (dedup / consensus families) | Reliable VAF floor (≥3 alt reads, binomial LOD) | Use |
|---|---|---|
| 30x WGS | ~10% | Clonal drivers only |
| 80–100x WES | ~5% | Clonal + major subclones |
| 500x panel (raw) | ~1% (error floor ~0.1–0.5%) | Subclones; CHIP confounds |
| ≥3000x raw → ≥1000 consensus families | 0.1–0.5% single-strand; 0.01–0.05% duplex | ctDNA, MRD |

Always sequence matched buffy coat/PBMC for ctDNA to exclude CHIP (DNMT3A, TET2, ASXL1, TP53, JAK2 at low VAF in plasma and blood).

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Tumour mean depth | ≥80–100x WES; ≥60x WGS; ≥500x panel | LOD for VAF ≥5% (WES) / ≥10% (WGS); raise for purity <30% |
| Normal mean depth | ≥30x (WGS/WES), ≥100x panel | Below this, germline hets leak as "somatic" (Strelka/Mutect2 normal evidence weak) |
| T/N relatedness (somalier) | ≈1.0 (fail <0.9) | Swap/contamination; tumour LOH lowers slightly |
| Cross-individual contamination | <2% (fail >3%) | `CalculateContamination`; FilterMutectCalls down-weights but cannot rescue |
| Tumour-in-normal (TiN) | <2% | Haematological or adjacent-normal samples; high TiN removes true variants — rescue via DeTiN or relax `normal_artifact` with caution |
| Purity | ≥20% for discovery; ≥30% for CN/clonality | Two-tool agreement ±0.1 |
| FFPE orientation artefact | `orientation` filter rate reported; high C>T at VAF <10% flags deamination | Require `--ob-priors`; treat SBS30/SBS1-like low-VAF excess as artefact |
| Ti/Tv somatic PASS | 1.0–3.0 WES (tissue-dependent; UV/APOBEC shift) | Outside = artefact signature or wrong normal |
| Caller concordance (PASS, VAF ≥10%) | ≥80% Mutect2∩Strelka2 | Lower → alignment/contamination problem |
| Indel fraction | 5–15% of PASS | >25% = homopolymer noise or MSI-H (confirm with MSIsensor-pro) |
| Mutations per Mb sanity | WES 0.5–50; >100 = POLE/MSI/UV or FP flood | Check against tissue priors (TCGA) |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `Sample X is not in BAM header` | `-normal`/`--normal-name` ≠ `@RG SM` | Read SM from header (step 2); sarek requires `<patient>_<sample>` |
| Thousands of low-VAF C>T (G>A) PASS calls | FFPE deamination, `--f1r2-tar-gz` omitted | Rerun with f1r2 + `LearnReadOrientationModel` + `--ob-priors` |
| `GetPileupSummaries` error about intervals | `-L` missing or different file from `-V` | `-V ${COMMON} -L ${COMMON}` |
| Contamination 0.00 with obvious germline leak | Pileup sites not polymorphic in this population / wrong build | Use `small_exac_common_3` matching build; check contig names |
| PASS count 10× expected in tumour-only | No PoN or PoN from another assay | Build PoN from ≥40 same-assay normals; exclude gnomAD AF >0.001 |
| Strelka misses all indels >~20bp | No Manta candidates | Pass `--indelCandidates candidateSmallIndels.vcf.gz` |
| Strelka PASS count near zero on WES | Depth filter | `--exome` |
| VarDict reports `Germline` as majority | Normal has tumour-in-normal or low depth | Check TiN; require normal ≥30x |
| Mutect2 `germline` filter removes a known driver | Driver also a gnomAD allele (e.g. common TP53 polymorphisms) | Inspect `GERMQ`; whitelist hotspots with evidence, document |
| ASCAT fails "no solution" | Low purity or wrong `gender`/`genomeVersion` | Set purity/ploidy bounds, verify sex from coverage (`genomics-qc`) |
| FACETS purity NA | `cval` too low, noisy logR | Raise `cval` to 300–500; filter `snp-pileup -r25,0` |
| TMB 3× published value | Divided by nominal panel size or included synonymous/germline | Use callable coding Mb from mosdepth; non-synonymous PASS only |
| MSIsensor-pro tumour-only calls everyone MSI-H | No assay-matched baseline | Build `baseline` from ≥20 normals; else require paired `msi` |
| PCGR crashes on tumour-only WGS | Unfiltered multi-million-variant VCF | Pre-filter PASS, gnomAD AF, exonic; `--vep_buffer_size` lower |
| somalier relatedness 0.5 | Different individuals or sample swap in LIMS | Stop; resolve identity before calling |
| DeepSomatic runs but all RefCall | Model type mismatch (WES data with WGS model, FFPE with fresh model) | Choose matching `--model_type` |
| Sentieon/Parabricks output sample order swapped | Names given in wrong flags | `bcftools query -l` and fix `--tumor_sample/--normal_sample` |

## nf-core

nf-core/sarek 3.12.0 somatic (pairs by `patient`; `status` 0 = normal, 1 = tumour; default 0; extra `status=1` rows = relapse/second tumour; `sex` optional, used by ASCAT):

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
P1,XX,0,P1_N,L1,P1_N_R1.fastq.gz,P1_N_R2.fastq.gz
P1,XX,1,P1_T,L1,P1_T_R1.fastq.gz,P1_T_R2.fastq.gz
```

```bash
nextflow run nf-core/sarek -r 3.12.0 -profile docker --input samplesheet.csv --outdir results \
  --genome GATK.GRCh38 --tools mutect2,strelka,ascat,cnvkit,manta,msisensorpro,vep
```

Valid somatic `--tools` tokens (3.12.0 schema): `mutect2, strelka, deepsomatic, parabricks_deepsomatic, parabricks_mutectcaller, sentieon_tnscope, freebayes, ascat, cnvkit, controlfreec, manta, tiddit, msisensorpro` (+ `msisensor2`, `lofreq`, `mpileup` tumour-only). Add `--wes --intervals targets.bed` for exomes; `--pon`/`--pon_tbi` for a custom PoN; `--joint_mutect2` for multi-tumour patients. Sarek sets `SM=<patient>_<sample>`; BAM/CRAM inputs must already carry that.

nf-core/oncoanalyser 3.0.0 (hmftools WiGiTS: SAGE, PAVE, AMBER, COBALT, PURPLE, LINX, CHORD, SIGS, CUPPA):

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1,P1,P1_N,normal,dna,bam,/data/P1_N.bam
P1,P1,P1_T,tumor,dna,bam,/data/P1_T.bam
```

```bash
nextflow run nf-core/oncoanalyser -r 3.0.0 -profile docker --mode wgts --genome GRCh38_hmf --input samplesheet.csv --outdir results
```

Note `GRCh38_hmf` is a distinct reference build (no alt contigs, hmf decoys); BAMs aligned to another GRCh38 flavour need realignment or `fastq` input. nf-core/fastquorum 2.0.0 for UMI consensus. Profiles, `-resume`, cluster config: `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

- Truth sets: SEQC2 HCC1395/HCC1395BL (WGS/WES, high-confidence somatic SNV/indel + callable BED; also the DeepSomatic quick-start subset) and COLO829/COLO829BL (hmftools). Scoring with `som.py`/`hap.py`/`truvari` and nf-core/variantbenchmarking: `.claude/skills/genomics-benchmark/SKILL.md`. Expected PASS F1 on SEQC2 WGS ≥0.9 SNV, ≥0.8 indel for Mutect2/Strelka2/DeepSomatic.
- In-silico spike-ins (BAMSurgeon) for assay-specific LOD; dilution series of HCC1395 into HCC1395BL for purity sensitivity curves.
- Sanity stats per run: PASS count vs tissue prior, Ti/Tv, indel fraction, VAF histogram (peak at purity/2 for clonal heterozygous), caller concordance, contamination, relatedness.
- Positive controls: known drivers in the SEQC2 truth VCF (look them up from the truth set, not memory) present at expected VAF; any hotspot whitelist must be versioned.
- Orthogonal confirmation (ddPCR/panel) for reported actionable variants with VAF <10% or depth <100.

## Additional resources

- `reference.md` (this dir): full Mutect2 chain script with scatter/merge, PoN script, ASCAT and FACETS R scripts, Sequenza commands, PCGR/CPSR invocations, fgbio UMI consensus chain; load when generating runnable scripts.
- `.claude/skills/genomics-cnv/SKILL.md`: somatic CNV/SV callers (CNVkit, Control-FREEC, Manta/GRIDSS SV), segmentation plots; load for copy-number mechanics beyond purity/ploidy.
- `.claude/skills/genomics-qc/SKILL.md`: somalier, coverage/duplication QC, FFPE/insert-size diagnostics; load before calling.
- `.claude/skills/genomics-benchmark/SKILL.md`: SEQC2/`som.py`/`truvari` benchmarking mechanics.
- `.claude/skills/genomics-annotation/SKILL.md`: VEP, germline ACMG, CPSR context for the normal.
- `.claude/skills/genomics-longread/SKILL.md`: ClairS/DeepSomatic PACBIO/ONT somatic modes.
- `.claude/skills/genomics-workflow/SKILL.md`: Nextflow/nf-core execution mechanics.
- `.claude/skills/genomics/reference.md`: GATK resource bundle (gnomAD af-only, PoN, `small_exac_common_3`), reference downloads.
