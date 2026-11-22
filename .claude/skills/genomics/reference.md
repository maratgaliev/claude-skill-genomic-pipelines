# References, annotations, databases and data retrieval

All URLs below returned HTTP 200 on 2026-11-22 unless marked `<check current: …>`. Re-check before quoting a release number: GENCODE, Ensembl, ClinVar and dbNSFP move every few months. Record the exact file name and checksum in `versions.txt`.

## 1. Choose the build before anything else

| Situation | Use | Why |
|---|---|---|
| New human short-read project | GRCh38 analysis set, **no alt contigs** | Avoids MAPQ 0 on alt-duplicated genes with non-alt-aware tools; what nf-core/GATK bundles assume |
| GATK/Broad-compatible resources (BQSR, VQSR, Mutect2 PoN) | Broad `Homo_sapiens_assembly38.fasta` (has alts + HLA + decoys) | All Broad VCF resources are on these contig names; use with alt-aware BWA (`bwa-mem2 mem` + `.alt` file) or restrict calling to primary contigs |
| Long reads, repeats, centromeres, population-specific | T2T-CHM13v2.0 | Complete; but ClinVar/gnomAD/VEP support is partial; lift results back to GRCh38 for interpretation |
| Legacy cohort already on GRCh37/hg19 | Stay on GRCh37 for that cohort; lift over only summary results | Realigning mid-project creates batch effects |
| Mouse | GRCm39 (= mm39) | Current since 2020; `mm10` = GRCm38 is legacy |

Never mix `chr1` (UCSC/Broad) and `1` (Ensembl/NCBI) styles inside one analysis. Pick one; convert at download time with `bcftools annotate --rename-chrs`, `samtools reheader`, or `sed` on BED/GTF; verify with `scripts/ref_sanity.sh`.

## 2. Human GRCh38

### Reference FASTA

```bash
set -euo pipefail
# NCBI analysis set, no alts, UCSC names (chr1..chrM) — default for new projects
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip -k GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38_no_alt.fa
samtools faidx GRCh38_no_alt.fa
gatk CreateSequenceDictionary -R GRCh38_no_alt.fa
```

Alternatives: Broad `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta` (+ `.fai`, `.dict`, `.64.alt` etc. in the same prefix; 3.2 GB uncompressed); GENCODE `GRCh38.primary_assembly.genome.fa.gz` (release-matched with the GTF, `chr` names, no alts); UCSC `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz`.

### Indexes

| Tool | Command | RAM / time / size |
|---|---|---|
| bwa-mem2 | `bwa-mem2 index -p GRCh38_no_alt GRCh38_no_alt.fa` | ≈70 GB / 1–2 h / 16 GB; index is **not** interchangeable with `bwa index` |
| bwa | `bwa index GRCh38_no_alt.fa` | 5 GB / 1 h / 5 GB |
| minimap2 | `minimap2 -d GRCh38_no_alt.mmi GRCh38_no_alt.fa` | 10 GB / 5 min |
| STAR | `STAR --runMode genomeGenerate --genomeDir star_150 --genomeFastaFiles GRCh38_no_alt.fa --sjdbGTFfile gencode.gtf --sjdbOverhang 149 --runThreadN 16` | 32 GB / 1 h / 30 GB; `sjdbOverhang = read length − 1` |
| Salmon | see `genomics-rnaseq/bulk.md` (decoy-aware) | 20 GB / 30 min |
| Bowtie2 | `bowtie2-build --threads 16 GRCh38_no_alt.fa GRCh38_no_alt` | 8 GB / 1–2 h |

### GATK resource bundle (Broad contig names == UCSC names; works with the NCBI no-alt FASTA for primary contigs)

```bash
set -euo pipefail
B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx \
         Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
         Homo_sapiens_assembly38.known_indels.vcf.gz Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
         hapmap_3.3.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz.tbi \
         1000G_omni2.5.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz.tbi \
         1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi; do
  wget -c "${B}/${f}"
done
# Somatic resources live in a different bucket
S=https://storage.googleapis.com/gatk-best-practices/somatic-hg38
for f in af-only-gnomad.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz.tbi 1000g_pon.hg38.vcf.gz 1000g_pon.hg38.vcf.gz.tbi \
         small_exac_common_3.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz.tbi; do
  wget -c "${S}/${f}"
done
```

`dbsnp138` is old but is what BQSR/VQSR were tuned with; for annotation use current dbSNP (below). The whole hg38/v0 prefix is >100 GB; download files, not the prefix.

### Annotation (GTF/transcripts)

| Source | File | Notes |
|---|---|---|
| GENCODE (preferred for human/mouse) | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz`; transcripts `gencode.v50.transcripts.fa.gz` | v50 is latest as of 2026-11-22 (`latest_release/` symlink). `chr` names. Strip `.N` version suffixes consistently or keep them everywhere |
| Ensembl | `https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` | Release 116 current; `1`-style names; must match VEP cache release |
| RefSeq | NCBI `GCF_000001405.40_GRCh38.p14_genomic.gtf.gz` | `NC_0000xx` names; rename with the assembly report (below) |
| MANE | `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/` | Select transcript set for clinical reporting |

Rule: GTF, transcript FASTA, Salmon/STAR index, VEP cache and `tx2gene` all come from the **same release**. Write the release into the output directory name (`salmon_gencode50/`).

### Variant databases

| Database | URL | Contig style | Notes |
|---|---|---|---|
| dbSNP (current) | `https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz` (+`.tbi`) | `NC_000001.11` | Rename: build a two-column map from `GCF_000001405.40_GRCh38.p14_assembly_report.txt` (columns 7→10: RefSeq accession → UCSC name) then `bcftools annotate --rename-chrs map.txt` |
| ClinVar | `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` (+`.tbi`, `.md5`) | `1`-style | Weekly; record the `##fileDate`; add `chr` with `--rename-chrs` for UCSC-style projects |
| gnomAD v4.1 | `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{1..22,X,Y}.vcf.bgz`; `.../genomes/gnomad.genomes.v4.1.sites.chrN.vcf.bgz` | `chr` | 730,947 exomes + 76,215 genomes; exomes ≈ 100 GB total — subset with `bcftools view -R targets.bed` or `bcftools annotate -x` to AF fields on download |
| gnomAD SV v4 / CNV v4 | `https://gnomad.broadinstitute.org/downloads#v4-structural-variants` | `chr` | For SV AF annotation (SVAFotate/AnnotSV) |
| VEP cache 116 | `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz` | n/a | ≈25 GB; `vep_install -a cf -s homo_sapiens -y GRCh38 --CACHE_VERSION 116` |
| ENCODE blacklist v2 | `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz`; `mm10-blacklist.v2.bed.gz` | `chr` | No official GRCm39 list from Boyle lab (404 on `mm39`); lift mm10 with `CrossMap.py bed` and say so |
| somalier sites | `https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz` | `chr` | For relatedness/sex/ancestry; see `genomics-qc` |
| GIAB truth sets | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/` | `chr` | Details in `genomics-benchmark` |
| Cytobands | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz` | `chr` | Plots, circos |

dbSNP rename map:

```bash
set -euo pipefail
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
grep -v '^#' GCF_000001405.40_GRCh38.p14_assembly_report.txt | awk -F'\t' '$10!="na"{print $7"\t"$10}' > refseq2ucsc.map
bcftools annotate --rename-chrs refseq2ucsc.map -Oz -o dbsnp_GRCh38_ucsc.vcf.gz GCF_000001405.40.gz
tabix -p vcf dbsnp_GRCh38_ucsc.vcf.gz
```

## 3. T2T-CHM13v2.0

```bash
set -euo pipefail
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz
```

Chain files for liftover to/from GRCh38 are under `.../T2T/CHM13/assemblies/chain/`. Expect partial support in VEP (use `--gff` + `--fasta`), none in ANNOVAR, and growing support in gnomAD (v4 has a CHM13 liftover subset).

## 4. Mouse GRCm39

| Item | URL |
|---|---|
| FASTA (Ensembl, `1`-style) | `https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz` |
| FASTA (UCSC, `chr`-style) | `https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz` |
| FASTA (NCBI) | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz` |
| GENCODE M39 GTF (`chr`-style) | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.primary_assembly.annotation.gtf.gz` (latest as of 2026-11-22; M40 not yet published) |
| dbSNP mouse | NCBI dbSNP dropped non-human; use Ensembl variation or MGP (Sanger Mouse Genomes Project) VCFs |

## 5. Pre-built bundles

| Bundle | Command | Caveat |
|---|---|---|
| AWS iGenomes (nf-core `--genome`) | `aws s3 --no-sign-request sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./GATK.GRCh38/` | Hundreds of GB; GTF/annotations are old (Ensembl 75-era for some keys). Use `--genome GATK.GRCh38` for sarek, `--genome GRCh38` for rnaseq, and prefer `--fasta/--gtf` with a pinned GENCODE for anything expression-related |
| refgenie | `refgenie init -c g.yaml; refgenie pull hg38/fasta hg38/bwa_index -c g.yaml` | Versioned asset digests; good for shared HPC |
| genomepy | `genomepy install GRCh38.p14 NCBI --annotation` | Quick, less control over exact files |
| Ensembl rapid/others | `https://ftp.ensembl.org/pub/release-116/` | Non-model organisms |

## 6. Public data retrieval

| Source | Command | Notes |
|---|---|---|
| SRA (fast path) | `fastq-dl --accession SRR000000 --provider ena --cpus 8` | Pulls ENA mirrored FASTQ directly; no SRA-format conversion |
| SRA (toolkit) | `prefetch --max-size 200G SRR000000 && fasterq-dump --split-files --threads 8 --outdir fastq SRR000000 && pigz fastq/*.fastq` | Needs `vdb-config --prefetch-to-cwd`; dbGaP data requires `vdb-config --import prj_NNNN.ngc` and is **controlled access** |
| ENA | `enaDataGet -f fastq -d fastq SRR000000` or `ffq SRR000000 \| jq` for URLs | ENA has md5 for every file: verify |
| GEO | `pysradb gse-to-srp GSE000000` → SRP → `pysradb srp-to-srr` → above | GEO supplementary files (counts) via `https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE000000/suppl/` |
| 1000 Genomes 30x GRCh38 | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/` | CRAMs on GRCh38 full analysis set (needs matching `.fa` for CRAM decoding) |
| GIAB raw data | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/` | Use for pipeline validation; see `genomics-benchmark` |
| ENCODE | `https://www.encodeproject.org/` with `xargs -L 1 curl -O -J -L < files.txt` | Batch download via the portal's `files.txt` |
| 10x datasets | `https://www.10xgenomics.com/datasets` | Chemistry listed per dataset; use for pipeline smoke tests |

Always: `md5sum -c` or compare to the provider's checksum; record accession, date and checksum in `versions.txt`.

## 7. Spike-in genomes (ChIP-Rx / CUT&RUN normalisation)

| Spike-in | FASTA | Notes |
|---|---|---|
| *D. melanogaster* dm6 (Active Motif ChIP-Rx) | `https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz` (+ `dm6.chrom.sizes`) | Build a combined index (`cat GRCh38_no_alt.fa dm6.fa`) with prefixed contig names, or align twice and count reads with MAPQ ≥ 30 per genome |
| *E. coli* K-12 MG1655 (CUT&RUN/CUT&Tag carry-over) | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz` (`NC_000913.3`) | Assembly `GCF_000005845.2`; scale factor = constant / spike-in reads; see `genomics-chipseq` |

## 8. Mappability, segmental duplications, common SNP sites (GATK CNV)

| Track | URL | Use |
|---|---|---|
| Umap k100 single-read mappability, hg38 | `https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz` (index: `https://bismap.hoffmanlab.org/`; Bismap k100 for bisulfite: `.../raw/hg38/k100.bismap.bed.gz`) | `gatk AnnotateIntervals --mappability-track` — bedGraph with score 0–1 in column 4; sort, bgzip, `tabix -p bed` |
| UCSC segmental duplications, hg38 | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz` | `--segmental-duplication-track` wants BED: `zcat genomicSuperDups.txt.gz \| cut -f2-4,27 \| sort -k1,1 -k2,2n \| bgzip > segdups.hg38.bed.gz && tabix -p bed segdups.hg38.bed.gz` (col 27 = fracMatch) |
| Common biallelic SNPs for `CollectAllelicCounts` / `GetPileupSummaries` | `small_exac_common_3.hg38.vcf.gz` (section 2, somatic bucket) | ExAC-derived, ≈60 k sites, AF-only; for denser allelic counts use `af-only-gnomad.hg38.vcf.gz` filtered with `bcftools view -m2 -M2 -v snps -i 'AF>0.05'` |

Both GATK tracks must be on the same contig names as the reference `.dict`; GATK ignores unindexed `.bed.gz`.

## 9. Region BEDs: PAR, tandem repeats, exclusion lists

GRCh38 PAR (NCBI/Ensembl, 1-based): PAR1 chrX:10,001–2,781,479 = chrY:10,001–2,781,479; PAR2 chrX:155,701,383–156,030,895 = chrY:56,887,903–57,217,415. Write as 0-based BED:

```bash
printf 'chrX\t10000\t2781479\tPAR1\nchrX\t155701382\t156030895\tPAR2\nchrY\t10000\t2781479\tPAR1\nchrY\t56887902\t57217415\tPAR2\n' > GRCh38_PAR.bed
```

Use it to exclude PAR from haploid-ploidy calls on XY samples (`-XL GRCh38_PAR.bed`) or to mask chrY PAR (Broad `Homo_sapiens_assembly38` already hard-masks chrY PAR to N).

| BED | URL | Use |
|---|---|---|
| Sniffles2 tandem-repeat annotations | `https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations` (`human_GRCh38_no_alt_analysis_set.trf.bed`) | `sniffles --tandem-repeats`; improves SV calls in repeats; `chr` names |
| DELLY exclude templates | `https://github.com/dellytools/delly/tree/main/excludeTemplates` (`human.hg38.excl.tsv`) | `delly call -x human.hg38.excl.tsv`; drops centromeres/telomeres/decoys |
| ENCODE blacklist v2 | section 2 (`hg38-blacklist.v2.bed.gz`) | Peak filtering (ChIP/ATAC) and CNV/SV FP removal |

## 10. Dosage sensitivity and SV/CNV population frequency

| Resource | URL | Notes |
|---|---|---|
| ClinGen dosage map (HI/TS scores) | `https://ftp.clinicalgenome.org/` — `ClinGen_gene_curation_list_GRCh38.tsv`, `ClinGen_region_curation_list_GRCh38.tsv`, `ClinGen_haploinsufficiency_gene_GRCh38.bed`, `ClinGen_triplosensitivity_gene_GRCh38.bed`, `ClinGen_recurrent_CNV_V2.1-hg38.bed` | Updated nightly; score 3 = sufficient evidence; ACMG/ClinGen CNV scoring (`genomics-cnv`) |
| gnomAD SV v4.1 (genomes) | `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz` (+ `.bed.gz`) | `chr` names; AF per SVTYPE; use with SVAFotate/AnnotSV |
| gnomAD CNV v4.1 (exomes) | `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz` | Rare coding CNVs (GATK-gCNV); no `.bed.gz` companion |
| Browser download page | `https://gnomad.broadinstitute.org/downloads` | Lists current SV/CNV file names and checksums |

## 11. Pathogenicity scores and liftover chains

| Resource | URL | Licence / notes |
|---|---|---|
| CADD v1.7 GRCh38 | `https://cadd.bihealth.org/download` → `https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz` (≈80 GB, + `.tbi`) | Free for academic/non-commercial; commercial via licence. VEP plugin `CADD` |
| dbNSFP 5.x | `https://www.dbnsfp.org/download` | Academic build free (registration); commercial build licensed. Aggregates REVEL, AlphaMissense, SIFT, PolyPhen etc.; prefer it over separate files when licence allows |
| SpliceAI precomputed scores | `https://basespace.illumina.com/s/otSPW8hnhaZR` (BaseSpace account; `spliceai_scores.raw.snv.hg38.vcf.gz`) | Illumina non-commercial licence; VEP plugin `SpliceAI` |
| AlphaMissense | `https://zenodo.org/records/8208688` (`AlphaMissense_hg38.tsv.gz`) | CC BY-NC-SA 4.0 — non-commercial; VEP plugin `AlphaMissense` |
| REVEL | `https://zenodo.org/records/7072866` (`revel-v1.3_all_chromosomes.zip`, GRCh37+GRCh38 positions) | Free for research; VEP plugin `REVEL` |
| UCSC liftover chains | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz`; `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz`; T2T: `.../hs1/liftOver/hs1ToHg38.over.chain.gz`, `.../hg38/liftOver/hg38ToHs1.over.chain.gz` (also section 3) | `CrossMap.py vcf chain.gz in.vcf ref.fa out.vcf` or `gatk LiftoverVcf`; lift results, not reads |

## 12. Metagenomics databases

| Database | URL | Notes |
|---|---|---|
| Kraken2/Bracken pre-built indexes | `https://benlangmead.github.io/aws-indexes/k2` | Standard (≈80 GB RAM), Standard-8/16 (capped), PlusPF, nt; Bracken `kmer_distrib` included; record the dated tarball name |
| GTDB-Tk reference data | `https://ecogenomics.github.io/GTDBTk/installing/index.html` → `https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz` | ≈110 GB unpacked; set `GTDBTK_DATA_PATH`; release must match the GTDB-Tk version |
| SILVA SSU/LSU rRNA | `<check current: https://www.arb-silva.de/download/>` (`SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`; HEAD timed out on 2026-11-22) | 16S/18S amplicons; DADA2/QIIME2-formatted trainsets are separate downloads |
| MetaPhlAn 4 | `metaphlan --install --bowtie2db ${DB}` pulls from `http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/` | HUMAnN requires the matching ChocoPhlAn/UniRef release |

## 13. Bisulfite / EM-seq conversion controls

| Control | Accession | Retrieval |
|---|---|---|
| Lambda phage (unmethylated) | `NC_001416.1` | `efetch -db nuccore -id NC_001416.1 -format fasta > lambda.fa` or `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=fasta&retmode=text` |
| pUC19 (CpG-methylated) | `L09137` | `efetch -db nuccore -id L09137 -format fasta > pUC19.fa` (returns `L09137.2`) |

Append both to the reference before `bismark_genome_preparation`; conversion rate = 1 − (methylated C / total C on lambda); pUC19 CpG methylation should be ≈ 95–100 %. efetch endpoints answer GET only (HEAD returns 405).

## 14. Provenance snippet

```bash
set -euo pipefail
OUT=results
mkdir -p "${OUT}"
{
  echo "date: $(date -u +%FT%TZ)"
  echo "reference: $(readlink -f "${REF}") md5=$(md5sum "${REF}" | cut -d' ' -f1)"
  echo "gtf: $(readlink -f "${GTF}") md5=$(md5sum "${GTF}" | cut -d' ' -f1)"
  for t in samtools bcftools gatk bwa-mem2 STAR salmon vep nextflow; do
    command -v "${t}" >/dev/null 2>&1 && bash "${CLAUDE_SKILL_DIR}/scripts/check_env.sh" "${t}"
  done
} >> "${OUT}/versions.txt"
```
