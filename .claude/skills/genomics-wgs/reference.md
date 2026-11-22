# genomics-wgs reference blocks

Supporting commands for `.claude/skills/genomics-wgs/SKILL.md`. All placeholders are `${UPPER_CASE}` declared at the top of each block. Resource VCFs (HapMap 3.3, Omni 2.5, 1000G phase1 high-confidence SNPs, Mills, Axiom Exome Plus, dbSNP) come from the Broad hg38 bundle; acquisition is in `.claude/skills/genomics/reference.md`.

## Scatter/gather intervals (`SplitIntervals`)

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SCATTER=24
OUTDIR=intervals
WGS_CALLING=/refs/GRCh38/wgs_calling_regions.hg38.interval_list
gatk SplitIntervals -R ${REF} -L ${WGS_CALLING} \
  --scatter-count ${SCATTER} --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
  -O ${OUTDIR}
ls ${OUTDIR}/*.interval_list | wc -l
```

`wgs_calling_regions.hg38.interval_list` (Broad bundle) excludes centromeres/N-gaps; without it pass `-L chr1 … chrY` or omit `-L` (then `--subdivision-mode INTERVAL_SUBDIVISION` splits chromosomes). WES: `-L targets.bed --interval-padding 100` and `--scatter-count 8`.

Per-shard HaplotypeCaller then gather:

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
DBSNP=/refs/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf
OUTDIR=intervals
printf '%s\n' ${OUTDIR}/*.interval_list | xargs -P 8 -I{} bash -c 'IL=$1; SHARD=$(basename "${IL}" .interval_list); gatk --java-options "-Xmx4g" HaplotypeCaller -R "$2" -I "$3".bqsr.bam -D "$4" -L "${IL}" -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation --native-pair-hmm-threads 2 -O "$3.${SHARD}.g.vcf.gz"' _ {} ${REF} ${SAMPLE} ${DBSNP}
gatk GatherVcfs $(for IL in ${OUTDIR}/*.interval_list; do printf -- '-I %s.%s.g.vcf.gz ' ${SAMPLE} $(basename ${IL} .interval_list); done) -O ${SAMPLE}.g.vcf.gz
tabix -f -p vcf ${SAMPLE}.g.vcf.gz
bcftools view -H ${SAMPLE}.g.vcf.gz | cut -f1 | uniq | wc -l
```

`xargs -P` exits non-zero if any shard fails (unlike a bare `&`/`wait` loop, which discards exit codes); keep `-P` ≤ RAM / 4 GB. `GatherVcfs` requires shards in genomic order; `SplitIntervals` names (`0000-scattered…`) sort correctly. Replace the `xargs` fan-out with your scheduler array (see `.claude/skills/genomics-workflow/SKILL.md`). `GatherVcfs` does not write an index — always `tabix` after.

### Haploid X/Y for XY samples (HaplotypeCaller)

GRCh38 PAR: `chrX:10001-2781479`, `chrX:155701383-156030895`, `chrY:10001-2781479`, `chrY:56887903-57217415`. Call autosomes + PAR with default ploidy, non-PAR X/Y with `--ploidy 1`, gather.

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
DBSNP=/refs/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf
printf 'chrX\t2781479\t155701382\nchrY\t2781479\t56887902\n' > nonpar_xy.bed
gatk --java-options "-Xmx4g" HaplotypeCaller -R ${REF} -I ${SAMPLE}.bqsr.bam -D ${DBSNP} \
  -L nonpar_xy.bed --ploidy 1 -ERC GVCF -O ${SAMPLE}.xy_haploid.g.vcf.gz
bcftools query -f '[%GT\n]' ${SAMPLE}.xy_haploid.g.vcf.gz | grep -c '/' || true
```

Exclude `nonpar_xy.bed` from the diploid run (`-XL nonpar_xy.bed`). Validation prints 0 when all X/Y non-PAR genotypes are haploid. Joint genotyping mixed-ploidy gVCFs works in GenotypeGVCFs; GLnexus for DeepVariant handles haploid contigs when DV was run with `--haploid_contigs`.

## VQSR (cohort ≥30 WES or ≥1 WGS; GATK 4.7.0.0)

Run SNP first, apply, then INDEL on the SNP-recalibrated VCF. Drop `-an DP` for WES. Use `-an AS_*` + `--use-allele-specific-annotations` only if the gVCFs were produced with `-G AS_StandardAnnotation` on every sample.

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
COHORT=cohort
B=/refs/GRCh38
HAPMAP=${B}/hapmap_3.3.hg38.vcf.gz
OMNI=${B}/1000G_omni2.5.hg38.vcf.gz
KG=${B}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
MILLS=${B}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
AXIOM=${B}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
DBSNP=${B}/Homo_sapiens_assembly38.dbsnp138.vcf
SNP_SENS=99.7
INDEL_SENS=99.7
gatk --java-options "-Xmx24g" VariantRecalibrator -R ${REF} -V ${COHORT}.raw.vcf.gz \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  -mode SNP --max-gaussians 6 \
  --resource:hapmap,known=false,training=true,truth=true,prior=15 ${HAPMAP} \
  --resource:omni,known=false,training=true,truth=true,prior=12 ${OMNI} \
  --resource:1000G,known=false,training=true,truth=false,prior=10 ${KG} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${DBSNP} \
  -O ${COHORT}.snp.recal --tranches-file ${COHORT}.snp.tranches \
  --rscript-file ${COHORT}.snp.plots.R
gatk --java-options "-Xmx8g" ApplyVQSR -R ${REF} -V ${COHORT}.raw.vcf.gz \
  --recal-file ${COHORT}.snp.recal --tranches-file ${COHORT}.snp.tranches \
  --truth-sensitivity-filter-level ${SNP_SENS} --create-output-variant-index true \
  -mode SNP -O ${COHORT}.snp_recal.vcf.gz
gatk --java-options "-Xmx24g" VariantRecalibrator -R ${REF} -V ${COHORT}.snp_recal.vcf.gz \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
  -mode INDEL --max-gaussians 4 \
  --resource:mills,known=false,training=true,truth=true,prior=12 ${MILLS} \
  --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${AXIOM} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${DBSNP} \
  -O ${COHORT}.indel.recal --tranches-file ${COHORT}.indel.tranches
gatk --java-options "-Xmx8g" ApplyVQSR -R ${REF} -V ${COHORT}.snp_recal.vcf.gz \
  --recal-file ${COHORT}.indel.recal --tranches-file ${COHORT}.indel.tranches \
  --truth-sensitivity-filter-level ${INDEL_SENS} --create-output-variant-index true \
  -mode INDEL -O ${COHORT}.filtered.vcf.gz
bcftools view -H -f PASS ${COHORT}.filtered.vcf.gz | wc -l
```

Validation: `grep -v '^#' ${COHORT}.snp.tranches | column -t -s,` — the chosen tranche should keep ≥ 99 % of known (dbSNP) sites with a novel Ti/Tv ≥ 1.9 (WGS); `bcftools view -H ${COHORT}.filtered.vcf.gz | cut -f7 | sort | uniq -c` shows filter distribution. Priors/annotation lists follow the GATK VariantRecalibrator documentation (https://gatk.broadinstitute.org/hc/en-us/articles/54698787665179-VariantRecalibrator) and the WARP JointGenotyping defaults; `--max-gaussians 4` for INDEL prevents the "no data in some Gaussians" failure on sparse indel sets. Large cohorts (>1000 WGS): scatter `ApplyVQSR` per interval and consider `--sample-every-Nth-variant 10` at SNP model-building; GATK's VETS (`ExtractVariantAnnotations`/`TrainVariantAnnotationsModel`/`ScoreVariantAnnotations`) is the newer alternative at biobank scale.

## Hard filters (small cohorts, panels, no truth resources)

Values from GATK "Hard-filtering germline short variants" (https://gatk.broadinstitute.org/hc/en-us/articles/360035890471). `QUAL < 30` is a floor; raise to 50 for noisy panels. Missing annotations (RankSum at hom-alt sites) pass by design.

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
VCF=cohort.raw.vcf.gz
OUT=cohort.hardfiltered.vcf.gz
gatk SelectVariants -R ${REF} -V ${VCF} --select-type-to-include SNP -O snps.vcf.gz
gatk VariantFiltration -R ${REF} -V snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O snps.filtered.vcf.gz
gatk SelectVariants -R ${REF} -V ${VCF} --select-type-to-include INDEL --select-type-to-include MIXED -O indels.vcf.gz
gatk VariantFiltration -R ${REF} -V indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O indels.filtered.vcf.gz
gatk MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O ${OUT}
bcftools view -H ${OUT} | cut -f7 | sort | uniq -c | sort -k1,1nr
```

Validation: filter counts above; `QD2` and `QUAL30` should dominate; if `MQRankSum`/`ReadPosRankSum` filters remove >5 % of SNPs, check for mapping artefacts (wrong reference build or missing decoys). GATK logs `JEXL engine: undefined variable` once per missing annotation — expected. Genotype-level filtering (`--genotype-filter-expression "DP < 10 || GQ < 20" --genotype-filter-name lowDPGQ --set-filtered-genotype-to-no-call`) is a separate decision; record it if applied. `MIXED` records (SNP+indel at one site) are filtered with the INDEL rules.

## NVScoreVariants + FilterVariantTranches (single-sample HaplotypeCaller VCF; GATK ≥4.6, verified in 4.7.0.0 source)

Replaces `CNNScoreVariants` (deprecated). Requires the GATK conda/PyTorch environment (`broadinstitute/gatk:4.7.0.0` Docker image includes it; not in `gatk-lite`). Annotation key written: `CNN_1D` (default `--tensor-type reference`) or `CNN_2D` (`--tensor-type read_tensor`, needs `-I` BAM). Score the genotyped VCF, not the gVCF. GPU via `--accelerator cuda` (default `auto`).

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
SAMPLE=NA12878
B=/refs/GRCh38
HAPMAP=${B}/hapmap_3.3.hg38.vcf.gz
MILLS=${B}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gatk --java-options "-Xmx8g" GenotypeGVCFs -R ${REF} -V ${SAMPLE}.g.vcf.gz -O ${SAMPLE}.vcf.gz
gatk NVScoreVariants -R ${REF} -V ${SAMPLE}.vcf.gz -O ${SAMPLE}.cnn1d.vcf.gz
gatk FilterVariantTranches -V ${SAMPLE}.cnn1d.vcf.gz \
  --resource ${HAPMAP} --resource ${MILLS} \
  --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 \
  --invalidate-previous-filters -O ${SAMPLE}.filtered.vcf.gz
bcftools view -H ${SAMPLE}.filtered.vcf.gz | cut -f7 | sort | uniq -c
```

2D model (reads-aware, slower, higher accuracy): `gatk NVScoreVariants -R ${REF} -V ${SAMPLE}.vcf.gz --tensor-type read_tensor -I ${SAMPLE}.bqsr.bam -O ${SAMPLE}.cnn2d.vcf.gz` then `--info-key CNN_2D`. Flags verified against `NVScoreVariants.java` on GATK master (2026-11-22): `-V -R -O`, `-I`, `--tensor-type {reference,read_tensor}`, `--batch-size`, `--random-seed`, `--tmp-file`, `--accelerator`. Tranche values are the GATK defaults for `FilterVariantTranches`; tighten (`99.9`/`99.0`) for precision-first use. Scatter by chromosome for WGS (`-L chr1`) and `GatherVcfs`.

## GLnexus (DeepVariant gVCF joint calling; v1.4.1)

```bash
set -euo pipefail
COHORT=cohort
CONFIG=DeepVariantWGS
THREADS=16
MEM_GB=32
ls *.dv.g.vcf.gz > ${COHORT}.gvcf_list.txt
docker run --rm -v "${PWD}":/data -w /data ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
  /usr/local/bin/glnexus_cli --config ${CONFIG} --threads ${THREADS} --mem-gbytes ${MEM_GB} \
  --dir /data/GLnexus.DB --list /data/${COHORT}.gvcf_list.txt \
  > ${COHORT}.glnexus.bcf
bcftools view --threads ${THREADS} -Oz -o ${COHORT}.glnexus.vcf.gz ${COHORT}.glnexus.bcf
bcftools index -t ${COHORT}.glnexus.vcf.gz
bcftools query -l ${COHORT}.glnexus.vcf.gz | wc -l
```

Config presets: `DeepVariantWGS`, `DeepVariantWES`, `DeepVariant_unfiltered` (no QC filters, for downstream custom filtering), `gatk`, `gatk_unfiltered`. WES: `--config DeepVariantWES --bed targets.bed`. `--dir` must not exist beforehand (`rm -rf GLnexus.DB` on rerun). Memory: `--mem-gbytes` ≈ 1 GB per 2–3 WGS samples up to node limit; disk ≈ 2× input gVCF total. Output genotypes carry `RNC` (reason-no-call) FORMAT; `GQ` recalibrated in v1.4.1 for 1/1 rare alleles. For >1000 samples, run per chromosome (`--bed chr1.bed`) and `bcftools concat`. Multiallelic sites are joint-genotyped as split records; run `bcftools norm -m -both` before annotation anyway for consistency with the GATK path.

## Trio / pedigree

Pedigree file `trio.ped` (tab-separated, no header): `FAM  CHILD  FATHER  MOTHER  SEX(1=M,2=F)  PHENO(1=unaffected,2=affected)`; one row per individual (parents with `0 0`). Sample IDs must match VCF `SM`.

1. Joint-genotype the trio as a cohort (GenomicsDB or GLnexus); do not merge single-sample VCFs.
2. Refine genotypes and flag de novos (GATK path):

```bash
set -euo pipefail
REF=/refs/GRCh38/Homo_sapiens_assembly38.fasta
TRIO=trio
PED=trio.ped
SUPPORTING=/refs/GRCh38/gnomad.genomes.v4.1.sites.autosomes.subset.vcf.gz
gatk CalculateGenotypePosteriors -R ${REF} -V ${TRIO}.filtered.vcf.gz -ped ${PED} \
  --supporting-callsets ${SUPPORTING} -O ${TRIO}.cgp.vcf.gz
gatk VariantFiltration -R ${REF} -V ${TRIO}.cgp.vcf.gz \
  --genotype-filter-expression "GQ < 20" --genotype-filter-name "lowGQ" -O ${TRIO}.cgp.gq.vcf.gz
gatk VariantAnnotator -R ${REF} -V ${TRIO}.cgp.gq.vcf.gz -ped ${PED} -A PossibleDeNovo -O ${TRIO}.denovo.vcf.gz
bcftools view -H -i 'INFO/hiConfDeNovo!=""' ${TRIO}.denovo.vcf.gz | wc -l
```

`--supporting-callsets` is optional (population priors; use a sites-only gnomAD subset, see `.claude/skills/genomics-annotation/SKILL.md` for the file). `PossibleDeNovo` writes `hiConfDeNovo` (all three GQ ≥ 20) and `loConfDeNovo` (GQ ≥ 10) INFO lists of child sample IDs. Expect ≈ 60–100 hiConf de novo SNVs genome-wide per 30× trio; hundreds ⇒ parent contamination/swap or low parental depth. Confirm relatedness first: `somalier relate --ped trio.ped somalier/*.somalier` (expected relatedness ≈ 0.5 parent–child, ≈ 0 between parents).

3. Mendelian error rate: `bcftools +mendelian2 -p ${PED} -m c ${TRIO}.filtered.vcf.gz` (bcftools ≥1.20; `+mendelian` in older) — < 1 % at PASS sites; higher ⇒ sample swap or over-permissive filtering. Count only; do not remove Mendelian-inconsistent sites before de novo review.

4. DeepVariant trios: DeepTrio (`google/deepvariant:deeptrio-1.10.0`, `run_deeptrio --model_type WGS --reads_child … --reads_parent1 … --reads_parent2 …`) then GLnexus `--config DeepVariant`-family preset; `PossibleDeNovo` via `VariantAnnotator -ped` still works on the GLnexus output after `bcftools norm`.

5. Sex-chromosome ploidy in pedigrees: `CalculateGenotypePosteriors` assumes diploid; skip non-PAR X/Y for male probands or run the autosomal and X-linked analyses separately (`-L` / `-XL nonpar_xy.bed`).

## GIAB samples for validation

| GIAB ID | Coriell | Ancestry / family | Use |
|---|---|---|---|
| HG001 | NA12878 | CEU female, pilot genome | Legacy comparisons; sarek `test_full_germline` uses it |
| HG002 | NA24385 | Ashkenazi son | Primary benchmark; richest truth set (v4.2.1 small variants, T2T-based Q100 draft) |
| HG003 / HG004 | NA24149 / NA24143 | Ashkenazi father / mother | Trio with HG002; de novo and Mendelian validation |
| HG005 | NA24631 | Han Chinese son | Second ancestry; trio with HG006/HG007 |
| HG006 / HG007 | NA24694 / NA24695 | Han Chinese father / mother | Trio with HG005 |

Truth VCF + high-confidence BED per sample: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/<AshkenazimTrio|ChineseTrio|NA12878_HG001>/<HGxxx>_…/NISTv4.2.1/GRCh38/` <check current: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/>. Public 30× Illumina FASTQ/BAM for HG001–HG007 exist in the GIAB FTP `data/` tree and in `s3://giab/`; acquisition and stratification BED downloads are in `.claude/skills/genomics/reference.md`; hap.py/vcfdist invocation and expected F1 ranges are in `.claude/skills/genomics-benchmark/SKILL.md`. Restrict comparison to the high-confidence BED intersected with your capture targets for WES.
