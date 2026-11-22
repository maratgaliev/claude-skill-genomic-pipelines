# genomics-somatic reference scripts

Runnable templates for `.claude/skills/genomics-somatic/SKILL.md`. Every script declares placeholders at the top; replace them, keep `set -euo pipefail`, record tool versions (`gatk --version`, etc.) in the run log. Resource files (`af-only-gnomad.hg38.vcf.gz`, `1000g_pon.hg38.vcf.gz`, `small_exac_common_3.hg38.vcf.gz`, interval lists) come from the GATK resource bundle; acquisition is covered in `.claude/skills/genomics/reference.md`.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## Mutect2 full chain (tumour-normal, scattered)

```bash
#!/usr/bin/env bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
GNOMAD=af-only-gnomad.hg38.vcf.gz
PON=1000g_pon.hg38.vcf.gz
COMMON=small_exac_common_3.hg38.vcf.gz
INTERVALS=wgs_calling_regions.hg38.interval_list
SAMPLE=P1
SCATTER=24
THREADS=4
OUT=${SAMPLE}.mutect2
mkdir -p ${OUT}

NSM=$(samtools view -H ${NBAM} | awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){sub("SM:","",$i); print $i; exit}}')
TSM=$(samtools view -H ${TBAM} | awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){sub("SM:","",$i); print $i; exit}}')
test -n "${NSM}" && test -n "${TSM}" && test "${NSM}" != "${TSM}"
bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" ${REF} ${TBAM} ${NBAM} ${GNOMAD} ${PON} ${COMMON}

gatk SplitIntervals -R ${REF} -L ${INTERVALS} --scatter-count ${SCATTER} -O ${OUT}/intervals
ls ${OUT}/intervals/*.interval_list | wc -l

for IL in ${OUT}/intervals/*.interval_list; do
  TAG=$(basename ${IL} .interval_list)
  gatk --java-options "-Xmx8g" Mutect2 -R ${REF} -I ${TBAM} -I ${NBAM} -normal ${NSM} \
    --germline-resource ${GNOMAD} --panel-of-normals ${PON} \
    --f1r2-tar-gz ${OUT}/${TAG}.f1r2.tar.gz -L ${IL} \
    --native-pair-hmm-threads ${THREADS} -O ${OUT}/${TAG}.unfiltered.vcf.gz &
done
wait
test "$(ls ${OUT}/*.unfiltered.vcf.gz | wc -l)" -eq "${SCATTER}"

F1R2_ARGS=$(for f in ${OUT}/*.f1r2.tar.gz; do printf -- "-I %s " "$f"; done)
gatk LearnReadOrientationModel ${F1R2_ARGS} -O ${OUT}/read-orientation-model.tar.gz
VCF_ARGS=$(for f in ${OUT}/*.unfiltered.vcf.gz; do printf -- "-I %s " "$f"; done)
gatk MergeVcfs ${VCF_ARGS} -O ${OUT}/${SAMPLE}.unfiltered.vcf.gz
STATS_ARGS=$(for f in ${OUT}/*.unfiltered.vcf.gz.stats; do printf -- "-stats %s " "$f"; done)
gatk MergeMutectStats ${STATS_ARGS} -O ${OUT}/${SAMPLE}.unfiltered.vcf.gz.stats
bcftools view -h ${OUT}/${SAMPLE}.unfiltered.vcf.gz | grep -q '^##tumor_sample='

gatk GetPileupSummaries -I ${TBAM} -V ${COMMON} -L ${COMMON} -O ${OUT}/tumour.pileups.table
gatk GetPileupSummaries -I ${NBAM} -V ${COMMON} -L ${COMMON} -O ${OUT}/normal.pileups.table
gatk CalculateContamination -I ${OUT}/tumour.pileups.table -matched ${OUT}/normal.pileups.table \
  --tumor-segmentation ${OUT}/segments.table -O ${OUT}/contamination.table
awk 'NR==2{printf "contamination=%s error=%s\n",$2,$3; if($2>0.03) exit 1}' ${OUT}/contamination.table

gatk FilterMutectCalls -R ${REF} -V ${OUT}/${SAMPLE}.unfiltered.vcf.gz \
  --stats ${OUT}/${SAMPLE}.unfiltered.vcf.gz.stats \
  --contamination-table ${OUT}/contamination.table --tumor-segmentation ${OUT}/segments.table \
  --ob-priors ${OUT}/read-orientation-model.tar.gz \
  --filtering-stats ${OUT}/${SAMPLE}.filteringStats.tsv -O ${OUT}/${SAMPLE}.filtered.vcf.gz
bcftools view -f PASS -Oz -o ${OUT}/${SAMPLE}.pass.vcf.gz ${OUT}/${SAMPLE}.filtered.vcf.gz
bcftools index -t ${OUT}/${SAMPLE}.pass.vcf.gz
bcftools +counts ${OUT}/${SAMPLE}.pass.vcf.gz
bcftools query -f '%FILTER\n' ${OUT}/${SAMPLE}.filtered.vcf.gz | tr ';' '\n' | sort | uniq -c | sort -rn
```

Optional realignment filter (slow; panels/hotspots):

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
IN=P1.mutect2/P1.filtered.vcf.gz
test -f ${REF%.fasta}.img || gatk BwaMemIndexImageCreator -I ${REF} -O ${REF%.fasta}.img
gatk FilterAlignmentArtifacts -R ${REF} -V ${IN} -I ${TBAM} --bwa-mem-index-image ${REF%.fasta}.img -O ${IN%.vcf.gz}.realigned.vcf.gz
bcftools view -H -i 'FILTER~"alignment"' ${IN%.vcf.gz}.realigned.vcf.gz | wc -l
```

Tumour-only variant of the chain: remove `-I ${NBAM} -normal ${NSM}` from Mutect2, `-matched` from CalculateContamination and the normal `GetPileupSummaries`; keep everything else. Post-filter additionally with `bcftools view -f PASS -i 'INFO/POPAF>=3'` (POPAF is −log10 gnomAD AF, so ≥3 keeps alleles with AF ≤0.001).

Per-sample tumour VAF/depth INFO tags for PCGR (`--tumor_dp_tag TDP --tumor_af_tag TAF --control_dp_tag NDP --control_af_tag NAF`):

```bash
set -euo pipefail
IN=P1.mutect2/P1.pass.vcf.gz
TSM=tumour
NSM=normal
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t%AF\t]\n' -s ${TSM},${NSM} ${IN} | cut -f1-8 | bgzip -c > tags.tsv.gz
tabix -s1 -b2 -e2 tags.tsv.gz
printf '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor depth">\n##INFO=<ID=TAF,Number=1,Type=Float,Description="Tumor AF">\n##INFO=<ID=NDP,Number=1,Type=Integer,Description="Normal depth">\n##INFO=<ID=NAF,Number=1,Type=Float,Description="Normal AF">\n' > tags.hdr
bcftools annotate -a tags.tsv.gz -h tags.hdr -c CHROM,POS,REF,ALT,TDP,TAF,NDP,NAF -Oz -o ${IN%.vcf.gz}.pcgr.vcf.gz ${IN}
bcftools index -t ${IN%.vcf.gz}.pcgr.vcf.gz
bcftools view -h ${IN%.vcf.gz}.pcgr.vcf.gz | grep -c 'ID=TAF'
```

## Panel of normals

```bash
#!/usr/bin/env bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
GNOMAD=af-only-gnomad.hg38.vcf.gz
INTERVALS=wgs_calling_regions.hg38.interval_list
NORMALS=normals.txt
PON_DB=pon_db
OUT=pon.hg38.vcf.gz
THREADS=4
test "$(wc -l < ${NORMALS})" -ge 40
mkdir -p pon_vcfs
while read -r BAM; do
  SM=$(samtools view -H ${BAM} | awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){sub("SM:","",$i); print $i; exit}}')
  gatk Mutect2 -R ${REF} -I ${BAM} --max-mnp-distance 0 -L ${INTERVALS} \
    --native-pair-hmm-threads ${THREADS} -O pon_vcfs/${SM}.vcf.gz
  test -s pon_vcfs/${SM}.vcf.gz.tbi
done < ${NORMALS}
V_ARGS=$(for f in pon_vcfs/*.vcf.gz; do printf -- "-V %s " "$f"; done)
gatk --java-options "-Xmx16g" GenomicsDBImport -R ${REF} -L ${INTERVALS} \
  --genomicsdb-workspace-path ${PON_DB} --merge-input-intervals true --reader-threads ${THREADS} ${V_ARGS}
gatk CreateSomaticPanelOfNormals -R ${REF} --germline-resource ${GNOMAD} -V gendb://${PON_DB} -O ${OUT}
bcftools view -h ${OUT} | grep -q '^##source=CreateSomaticPanelOfNormals'
bcftools view -H ${OUT} | wc -l
```

`normals.txt` = one BAM/CRAM path per line, same assay, aligner and reference as the cases, no tumour-bearing samples, no case-matched normals of the cohort being called. Use `--intervals targets.interval_list` for exomes. Default `CreateSomaticPanelOfNormals --min-sample-count 2` keeps sites seen in ≥2 normals.

## ASCAT (R, HTS input)

```r
suppressPackageStartupMessages(library(ASCAT))
tumour_bam  <- "tumour.bam"
normal_bam  <- "normal.bam"
tumour_name <- "P1_T"
normal_name <- "P1_N"
sex         <- "XX"
allelecounter <- Sys.which("alleleCounter")
stopifnot(nzchar(allelecounter))
loci_prefix    <- "G1000_loci_hg38_chr"
alleles_prefix <- "G1000_alleles_hg38_chr"
gc_file  <- "GC_G1000_hg38.txt"
rt_file  <- "RT_G1000_hg38.txt"
nthreads <- 12
ascat.prepareHTS(
  tumourseqfile = tumour_bam, normalseqfile = normal_bam,
  tumourname = tumour_name, normalname = normal_name,
  allelecounter_exe = allelecounter,
  alleles.prefix = alleles_prefix, loci.prefix = loci_prefix,
  gender = sex, genomeVersion = "hg38", nthreads = nthreads,
  tumourLogR_file = paste0(tumour_name, "_LogR.txt"),
  tumourBAF_file  = paste0(tumour_name, "_BAF.txt"),
  normalLogR_file = paste0(normal_name, "_LogR.txt"),
  normalBAF_file  = paste0(normal_name, "_BAF.txt"),
  min_base_qual = 10, additional_allelecounter_flags = "-f 0")
ascat.bc <- ascat.loadData(
  Tumor_LogR_file = paste0(tumour_name, "_LogR.txt"), Tumor_BAF_file = paste0(tumour_name, "_BAF.txt"),
  Germline_LogR_file = paste0(normal_name, "_LogR.txt"), Germline_BAF_file = paste0(normal_name, "_BAF.txt"),
  gender = sex, genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = gc_file, replictimingfile = rt_file)
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc <- ascat.aspcf(ascat.bc, penalty = 70)
ascat.plotSegmentedData(ascat.bc)
ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, write_segments = TRUE)
qc <- ascat.metrics(ascat.bc, ascat.output)
stopifnot(length(ascat.output$aberrantcellfraction) == 1)
write.table(data.frame(sample = tumour_name,
                       purity = ascat.output$aberrantcellfraction,
                       ploidy = ascat.output$ploidy,
                       goodness_of_fit = ascat.output$goodnessOfFit),
            paste0(tumour_name, ".ascat.purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ascat.output$segments, paste0(tumour_name, ".ascat.segments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, paste0(tumour_name, ".ascat.metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
```

Run: `Rscript ascat.R`. Validate: `goodness_of_fit ≥ 90`, purity in (0.1, 1], ploidy 1.5–5; inspect `*.sunrise.png` for a single sharp optimum. Loci/alleles/GC/RT files per build from the ASCAT `ReferenceFiles` documentation (`<check current: https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles>`). WES: add `BED_file =` to `ascat.prepareHTS`; targeted panels: run `ascat.prepareTargetedSeq` on a normal batch first and use `penalty = 25`. Tumour-only: `ascat.predictGermlineGenotypes` instead of normal BAF (lower confidence). Segments go to PCGR `--input_cna` as `Chromosome Start End nMajor nMinor`.

## FACETS (snp-pileup + R)

```bash
set -euo pipefail
DBSNP=dbsnp_146.hg38.vcf.gz
TBAM=tumour.bam
NBAM=normal.bam
SAMPLE=P1
snp-pileup -g -q15 -Q20 -P100 -r25,0 ${DBSNP} ${SAMPLE}.snp_pileup.csv.gz ${NBAM} ${TBAM}
test -s ${SAMPLE}.snp_pileup.csv.gz && zcat ${SAMPLE}.snp_pileup.csv.gz | head -2
```

```r
suppressPackageStartupMessages(library(facets))
sample <- "P1"
cval_pre <- 25
cval <- 150
set.seed(1234)
rcmat <- readSnpMatrix(paste0(sample, ".snp_pileup.csv.gz"))
xx <- preProcSample(rcmat, gbuild = "hg38", cval = cval_pre, snp.nbhd = 250)
oo <- procSample(xx, cval = cval, min.nhet = 15)
fit <- emcncf(oo)
stopifnot(!is.na(fit$purity))
png(paste0(sample, ".facets.png"), width = 1400, height = 1000)
plotSample(x = oo, emfit = fit)
dev.off()
write.table(data.frame(sample = sample, purity = fit$purity, ploidy = fit$ploidy, dipLogR = oo$dipLogR, loglik = fit$loglik),
            paste0(sample, ".facets.purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fit$cncf, paste0(sample, ".facets.cncf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
```

`cval` 150 for WES, 300–500 for WGS (higher = fewer, longer segments); `emcncf` columns `tcn.em`/`lcn.em` are total/minor copy number (PCGR: nMajor = tcn.em − lcn.em, nMinor = lcn.em). `lcn.em = NA` means minor CN unresolved, not 0. Order of BAMs to `snp-pileup` is normal first, then tumour.

## Sequenza

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
SAMPLE=P1
test -s ${REF%.fasta}.gc50.wig.gz || sequenza-utils gc_wiggle -w 50 -f ${REF} -o ${REF%.fasta}.gc50.wig.gz
sequenza-utils bam2seqz -n ${NBAM} -t ${TBAM} --fasta ${REF} -gc ${REF%.fasta}.gc50.wig.gz -o ${SAMPLE}.seqz.gz
sequenza-utils seqz_binning -w 50 -s ${SAMPLE}.seqz.gz -o ${SAMPLE}.small.seqz.gz
tabix -f -s1 -b2 -e2 -S1 ${SAMPLE}.small.seqz.gz && zcat ${SAMPLE}.small.seqz.gz | head -2
```

```r
suppressPackageStartupMessages(library(sequenza))
sample <- "P1"
ext <- sequenza.extract(paste0(sample, ".small.seqz.gz"), assembly = "hg38", chromosome.list = paste0("chr", c(1:22, "X")))
cp <- sequenza.fit(ext)
sequenza.results(sequenza.extract = ext, cp.table = cp, sample.id = sample, out.dir = paste0(sample, "_sequenza"))
```

Outputs `P1_sequenza/P1_confints_CP.txt` (cellularity, ploidy with CIs) and `P1_segments.txt` (`A`/`B` allele CN). Parallelise `bam2seqz` with `-C chr1 chr2 ...` per chromosome then concatenate.

## PURPLE via hmftools (manual; oncoanalyser does this end-to-end)

```bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
TBAM=tumour.bam
NBAM=normal.bam
TSM=P1_T
NSM=P1_N
HMF=hmf_ref_38
THREADS=8
java -jar amber.jar -reference ${NSM} -reference_bam ${NBAM} -tumor ${TSM} -tumor_bam ${TBAM} \
  -ref_genome_version 38 -loci ${HMF}/AmberGermlineSites.38.tsv.gz -output_dir amber -threads ${THREADS}
test -s amber/${TSM}.amber.baf.tsv.gz
java -jar cobalt.jar -reference ${NSM} -reference_bam ${NBAM} -tumor ${TSM} -tumor_bam ${TBAM} \
  -ref_genome ${REF} -gc_profile ${HMF}/GC_profile.1000bp.38.cnp -output_dir cobalt -threads ${THREADS}
test -s cobalt/${TSM}.cobalt.ratio.tsv.gz
java -jar purple.jar -reference ${NSM} -tumor ${TSM} -amber_dir amber -cobalt_dir cobalt \
  -somatic_vcf P1.mutect2/P1.pass.vcf.gz -ref_genome ${REF} -ref_genome_version 38 \
  -gc_profile ${HMF}/GC_profile.1000bp.38.cnp -ensembl_data_dir ${HMF}/ensembl_data_cache -output_dir purple -threads ${THREADS}
cut -f1-6 purple/${TSM}.purple.purity.tsv
```

hmftools expects SAGE/PAVE-annotated VCFs for full feature use; with Mutect2 input, PURPLE still fits purity/ploidy but driver/TMB fields are approximate. Reference bundle (`AmberGermlineSites`, `GC_profile`, ensembl cache): `<check current: https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md>`.

## PCGR / CPSR

```bash
set -euo pipefail
SAMPLE=P1
VEP_DIR=${HOME}/.vep
PCGR_DATA=/refdata/pcgr
TUMOR_SITE=15
pcgr --input_vcf ${SAMPLE}.mutect2/${SAMPLE}.pass.pcgr.vcf.gz --input_cna ${SAMPLE}.cna_segments.tsv \
  --vep_dir ${VEP_DIR} --refdata_dir ${PCGR_DATA} --output_dir pcgr_${SAMPLE} \
  --genome_assembly grch38 --sample_id ${SAMPLE} --assay WES --effective_target_size_mb 34 \
  --tumor_site ${TUMOR_SITE} --tumor_purity 0.6 --tumor_ploidy 2.1 \
  --tumor_dp_tag TDP --tumor_af_tag TAF --control_dp_tag NDP --control_af_tag NAF \
  --tumor_dp_min 20 --tumor_af_min 0.05 --control_dp_min 10 --control_af_max 0.02 \
  --estimate_tmb --tmb_high_threshold 10 --estimate_msi --estimate_signatures --estimate_hrd \
  --input_msi msi/${SAMPLE} --input_msi_caller msisensor-pro --vep_n_forks 8 --vcfanno_n_proc 8
ls pcgr_${SAMPLE}/${SAMPLE}.pcgr.grch38.html pcgr_${SAMPLE}/${SAMPLE}.pcgr.grch38.snv_indel.tsv.gz
cpsr --input_vcf ${SAMPLE}.germline.vcf.gz --vep_dir ${VEP_DIR} --refdata_dir ${PCGR_DATA} \
  --output_dir cpsr_${SAMPLE} --genome_assembly grch38 --sample_id ${SAMPLE} --panel_id 0 --vep_n_forks 8
ls cpsr_${SAMPLE}/${SAMPLE}.cpsr.grch38.html
```

PCGR 2.3.2 (2026-08-10). `--input_cna` columns: `Chromosome Start End nMajor nMinor`. `--effective_target_size_mb` must be the callable coding size measured with mosdepth for panels (`--assay TARGETED`). OncoKB: add `--oncokb_api_token ${ONCOKB_TOKEN} --oncokb_oncotree_code LUAD` (internet access; MSK licence for commercial/clinical use). Tumour-only: `--tumor_only` plus `--exclude_likely_het_germline --exclude_likely_hom_germline --exclude_dbsnp_nonsomatic`, and pre-filter WGS inputs to PASS + exonic. CPSR `--panel_id 0` = CPSR superpanel; use a Genomics England PanelApp id for a focused panel. Both need a refdata bundle matched to the PCGR version and the VEP cache release it was built against (`<check current: https://sigven.github.io/pcgr/articles/installation.html>`).

## fgbio UMI consensus chain (single-strand UMI, inline 8 bp on R1)

```bash
#!/usr/bin/env bash
set -euo pipefail
REF=Homo_sapiens_assembly38.fasta
R1=sample_R1.fastq.gz
R2=sample_R2.fastq.gz
SAMPLE=P1_cfDNA
LIBRARY=lib1
PU=FLOWCELL.1
READ_STRUCTURE="8M+T +T"
THREADS=16
MIN_READS=3
fgbio -Xmx1g --compression 1 --async-io FastqToBam --input ${R1} ${R2} --read-structures ${READ_STRUCTURE} \
  --sample ${SAMPLE} --library ${LIBRARY} --platform-unit ${PU} --output ${SAMPLE}.unmapped.bam
samtools view ${SAMPLE}.unmapped.bam | head -1 | grep -q 'RX:Z:'
samtools fastq ${SAMPLE}.unmapped.bam \
  | bwa mem -t ${THREADS} -p -K 150000000 -Y ${REF} - \
  | fgbio -Xmx4g --compression 1 --async-io ZipperBams --unmapped ${SAMPLE}.unmapped.bam --ref ${REF} --output ${SAMPLE}.mapped.bam
samtools quickcheck ${SAMPLE}.mapped.bam
samtools sort --template-coordinate --threads 8 -o ${SAMPLE}.mapped.sorted.bam ${SAMPLE}.mapped.bam
fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi --input ${SAMPLE}.mapped.sorted.bam \
  --strategy Adjacency --edits 1 --output ${SAMPLE}.grouped.bam --family-size-histogram ${SAMPLE}.family_sizes.txt
awk 'NR>1 && $1>=3 {s+=$3} END {printf "fraction_reads_in_families_ge3=%.3f\n", s}' ${SAMPLE}.family_sizes.txt
fgbio -Xmx4g --compression 1 CallMolecularConsensusReads --input ${SAMPLE}.grouped.bam \
  --output ${SAMPLE}.cons.unmapped.bam --min-reads 1 --min-input-base-quality 20 --threads 4
samtools fastq ${SAMPLE}.cons.unmapped.bam \
  | bwa mem -t ${THREADS} -p -K 150000000 -Y ${REF} - \
  | fgbio -Xmx4g --compression 1 --async-io ZipperBams --unmapped ${SAMPLE}.cons.unmapped.bam --ref ${REF} \
      --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${SAMPLE}.cons.mapped.bam
samtools quickcheck ${SAMPLE}.cons.mapped.bam
fgbio -Xmx8g --compression 0 FilterConsensusReads --input ${SAMPLE}.cons.mapped.bam --output /dev/stdout --ref ${REF} \
  --min-reads ${MIN_READS} --min-base-quality 45 --max-base-error-rate 0.2 \
  | samtools sort --threads 8 -o ${SAMPLE}.cons.filtered.bam --write-index
samtools quickcheck ${SAMPLE}.cons.filtered.bam
samtools view -c ${SAMPLE}.cons.filtered.bam
```

Duplex (e.g. IDT xGen duplex, `3M2S+T 3M2S+T`): `GroupReadsByUmi --strategy Paired`, `CallDuplexConsensusReads --min-reads 1 --min-input-base-quality 20`, then `FilterConsensusReads --min-reads 3 1 1 --min-base-quality 45 --max-base-error-rate 0.2 0.1 0.1`. Read-name UMIs (Illumina `:UMI` in the header) use `--extract-umis-from-read-names` in `FastqToBam`. Calling on `${SAMPLE}.cons.filtered.bam`: VarDict `-f 0.001`, or Mutect2 with `--min-allele-fraction 0.001 --max-reads-per-alignment-start 0` and a UMI-assay PoN; `cD`/`cM` tags give per-read family depth and error counts for post-filtering. Raw-read duplicate rate ≥80% at 3000x is expected; consensus depth, not raw depth, sets the LOD. nf-core/fastquorum 2.0.0 runs this chain with the same defaults.

## Consensus VAF / LOD helper

```python
from math import comb
def min_depth_for_vaf(vaf, min_alt=3, power=0.95):
    d = min_alt
    while True:
        p_detect = 1 - sum(comb(d, k) * vaf**k * (1 - vaf)**(d - k) for k in range(min_alt))
        if p_detect >= power:
            return d
        d += 1
for v in (0.10, 0.05, 0.01, 0.005, 0.001):
    print(f"VAF {v:.3f}: >= {min_depth_for_vaf(v)}x for 95% power with >=3 alt reads")
```

Use this to justify depth targets in the QC table; it ignores the error floor, which for raw Illumina reads is ~0.1–0.5% per base and must be measured on the matched normal / background panel.
