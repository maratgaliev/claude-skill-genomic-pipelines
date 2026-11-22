---
name: genomics-metagenomics
description: Microbiome and metagenome analysis — shotgun taxonomic profiling (Kraken2/Bracken, MetaPhlAn 4, sylph), host depletion, assembly/binning/MAG QC (MEGAHIT, metaSPAdes, MetaBAT2, SemiBin2, DAS Tool, CheckM2, GTDB-Tk), functional and AMR profiling (HUMAnN 3, eggNOG-mapper, AMRFinderPlus, RGI), 16S/ITS amplicons (QIIME 2, DADA2, SILVA) and compositional diversity statistics. Use when inputs are metagenomic or amplicon FASTQs, `*.kreport`, `*.bracken`, MetaPhlAn profiles, bins/MAGs, `*.qza`/`*.qzv`, feature tables, or the user says microbiome, 16S, MAG, Kraken, taxprofiler, mag, ampliseq. Not for host genomes (`genomics-wgs`), long-read microbial assembly polishing (`genomics-longread`) or RNA-seq of a single organism (`genomics-rnaseq`).
argument-hint: "[profile|assemble|amplicon|amr|diversity] [path]"
paths:
  - "**/*.kreport*"
  - "**/*.bracken"
  - "**/*.qza"
  - "**/*.qzv"
  - "**/*_profile.tsv"
---

# Genomics: metagenomics and microbiome

Expert in shotgun and amplicon microbiome pipelines: routing by data type, host read removal, read-based taxonomic profiling with explicit database provenance, assembly-based MAG recovery and QC, functional/AMR annotation, amplicon denoising, and compositional statistics with negative-control contamination handling. Outputs are commands, decision tables and validation checks.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" kraken2 bracken metaphlan sylph hostile bowtie2 megahit spades.py metabat2 SemiBin2 DAS_Tool checkm2 gtdbtk humann emapper.py amrfinder rgi qiime nextflow Rscript || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

Human reads removed from microbiome data are still controlled-access material: keep host-depleted FASTQs, never the host-matched reads, and record the depletion reference.

## Inspect first

| Check | Command | Infers |
| --- | --- | --- |
| Amplicon vs shotgun | `peek.sh R1.fastq.gz 400 \| awk 'NR%4==2{print substr($1,1,20)}' \| sort \| uniq -c \| sort -rn \| head -3` | One dominant 5′ k-mer (primer, e.g. `GTGYCAGCMGCCGCGGTAA` 515F, `CCTACGGGNGGCWGCAG` 341F) = amplicon; diverse = shotgun. |
| Read length / platform | `peek.sh R1.fastq.gz 8`; header `@M0` (MiSeq), `@A0`/`@LH` (NovaSeq/X), ONT UUIDs | 2×250/2×300 MiSeq → 16S V3–V4; 2×150 NovaSeq → shotgun; ONT → long-read route. |
| Host fraction (quick) | `kraken2 --db ${DB} --report /dev/stdout --output /dev/null --paired R1 R2 \| awk '$4=="S" && $6~/Homo sapiens/'` on a 100k-read subsample | % human decides depletion stringency and whether data are low-biomass. |
| Existing reports | `peek.sh x.kreport 15` | Kraken2 report (6 cols, `U/R/D/K/P…`) vs Bracken (`name taxid lvl kraken_assigned added new_est frac`) vs MetaPhlAn (`#mpa_v…` header line names the database). |
| Database provenance | `cat ${DB}/*.md5 ${DB}/library_report.tsv 2>/dev/null \| head`; `metaphlan --version`; `ls ${GTDBTK_DATA_PATH}` | Which index/date was used; never merge profiles from different database builds. |
| Negative controls present | sample names matching `blank|neg|ntc|extraction`; metadata column | Required for low-biomass (skin, tissue, placenta, air, cfDNA); absence must be reported as a limitation. |
| Assembly / bins | `grep -c '>' assembly.fa; seqkit stats bins/*.fa` | Contig count, N50, bin sizes before any QC. |
| QIIME artifacts | `qiime tools peek x.qza` | Semantic type (`FeatureTable[Frequency]`, `SampleData[PairedEndSequencesWithQuality]`), provenance. |

Ask only if still ambiguous: primer pair/region for amplicons, host species if not human, sample type (biomass), study design (paired/longitudinal), and whether MAGs or read-level profiles are the deliverable.

## Decision table

| Situation | Choose | Why |
| --- | --- | --- |
| Short-read shotgun, taxonomy, fast | Kraken2 + Bracken (Standard or PlusPF) | k-mer; high recall; needs Bracken for abundance; database RAM-bound. |
| Shotgun, precise relative abundance, human gut/body sites | MetaPhlAn 4 (SGB database) | Marker genes; low false positives; profiles sum to 100 % of known; misses novel taxa outside SGB set. |
| Shotgun, species-level with ANI/containment, huge DBs, low RAM | sylph (`sylph profile`) | Sketch-based; GTDB r232 DB; reports coverage and ANI; no reads-per-taxon output. |
| Database | Human-associated: Standard/PlusPF (RefSeq) or MetaPhlAn SGB; environmental/novel: GTDB (Kraken2 GTDB index or sylph GTDB); viruses: add PlusPF/Viral; eukaryotic pathogens: PlusPF/PlusPFP | Database determines what can be named; record date string in every output. |
| Kraken2 `--confidence` | 0.1 default for Standard on 2×150; 0.05 for 2×75; ≥0.2 when FP species matter (clinical/low-biomass) | Raises precision, lowers classified %; always re-run Bracken after changing it. |
| Host removal | `hostile clean` (default index T2T-CHM13v2.0 + HLA) or Bowtie2 `--very-sensitive` against GRCh38 **and** T2T-CHM13 | GRCh38 alone leaves ~0.1–0.5 % human reads (centromeres, HLA); for non-human hosts build the host index. |
| Assembly | MEGAHIT (large/complex, memory-bound); metaSPAdes (higher N50, ≤~100 M read pairs); co-assembly only with `group` logic | metaSPAdes RAM ≈ 1 GB per 1 M read pairs rough upper bound; MEGAHIT ≈ 1/4 of that. |
| Binning | MetaBAT2 + SemiBin2 + MaxBin2 → DAS Tool (or MAGScoT) refinement; CheckM2 for quality; GTDB-Tk `classify_wf` for taxonomy | Ensemble refinement recovers 20–40 % more MQ/HQ bins than a single binner. |
| MAG quality classes | HQ: ≥90 % completeness, ≤5 % contamination, 5S/16S/23S + ≥18 tRNAs (MIMAG); MQ: ≥50 %/≤10 % | Report MIMAG class per bin; never aggregate LQ bins into abundance tables. |
| Functional profiling | HUMAnN 3.9 (UniRef90, with MetaPhlAn 4 taxonomy) for pathways; eggNOG-mapper on predicted proteins (Prodigal/Pyrodigal) for MAGs | HUMAnN 4 is alpha (EC-filtered UniRef90 only; incompatible with MetaPhlAn ≥4.2); use 3.9 for production. |
| AMR | AMRFinderPlus (`amrfinder -u` then `--plus`) on contigs/proteins; RGI (CARD) on contigs; abricate for quick multi-DB screens | AMRFinderPlus curated, point mutations for named organisms (`--organism`); RGI adds CARD ontology; reads-level AMR needs ≥ several × coverage — report as presence only. |
| 16S/ITS amplicons | QIIME 2 amplicon distribution (2026.7 or later) with DADA2 denoising; SILVA 138.2 or 144 classifier matching the release; nf-core/ampliseq 3.1.0 | ASVs over OTUs; classifier must be trained for the primer region or use full-length weighted classifier. |
| Diversity statistics | `vegan`/`phyloseq`/`mia` (TreeSummarizedExperiment); alpha on rarefied **or** coverage-based; beta on Aitchison (CLR) or Bray–Curtis with rarefaction; differential abundance: ANCOM-BC2, LinDA, ALDEx2 | Compositionality: relative abundances are not independent; rarefaction debate — rarefy for presence/absence and alpha, use CLR/compositional models for abundance. |
| Contamination | `decontam` (frequency + prevalence with negative controls); remove taxa prevalent in blanks; report kit-ome genera (`Ralstonia`, `Bradyrhizobium`, `Cutibacterium` in low-biomass) | Low-biomass samples are dominated by reagent contaminants; filtering without blanks is guesswork. |
| Long-read metagenomes (ONT/PacBio HiFi) | Taxonomy: Kraken2 (`--minimum-hit-groups 3`), sylph, MetaPhlAn ≥4.2 long-read mode; assembly: metaFlye/hifiasm-meta | Short-read binning cut-offs change; see `.claude/skills/genomics-longread/SKILL.md` for basecalling and polishing. |

## Workflow

Variables used throughout: `SAMPLE`, `THREADS`, `K2DB` (Kraken2 index dir), `DBVER` (date string recorded from the index).

1. Host depletion (human; for other hosts build the Bowtie2 index from that genome).

```bash
set -euo pipefail
SAMPLE=S1
THREADS=16
mkdir -p clean
hostile clean --fastq1 "raw/${SAMPLE}_R1.fastq.gz" --fastq2 "raw/${SAMPLE}_R2.fastq.gz" --index human-t2t-hla --threads "${THREADS}" --out-dir clean > "clean/${SAMPLE}.hostile.json"
```
Validate: `jq '.[0] | {reads_in, reads_out, reads_removed_proportion}' clean/${SAMPLE}.hostile.json` (proportion removed consistent with sample type: stool <1 %, saliva 50–90 %, tissue >95 %). Bowtie2 alternative: `bowtie2 --very-sensitive -p ${THREADS} -x ref/GRCh38_T2T -1 R1 -2 R2 --un-conc-gz clean/${SAMPLE}_R%.fastq.gz -S /dev/null`, then validate with `zcat clean/${SAMPLE}_R1.fastq.gz | awk 'END{print NR/4" pairs kept"}'`. Raw-read QC/trimming (fastp) lives in `.claude/skills/genomics-qc/SKILL.md`.

2. Kraken2 + Bracken profiling.

```bash
set -euo pipefail
SAMPLE=S1
THREADS=16
K2DB=db/k2_standard_20261022
READLEN=150
DBVER="$(basename "${K2DB}")"
kraken2 --db "${K2DB}" --threads "${THREADS}" --paired --confidence 0.1 --minimum-hit-groups 2 --report "tax/${SAMPLE}.kreport" --output "tax/${SAMPLE}.kraken2.tsv" --gzip-compressed "clean/${SAMPLE}_R1.fastq.gz" "clean/${SAMPLE}_R2.fastq.gz"
bracken -d "${K2DB}" -i "tax/${SAMPLE}.kreport" -o "tax/${SAMPLE}.S.bracken" -w "tax/${SAMPLE}.S.breport" -r "${READLEN}" -l S -t 10
echo -e "${SAMPLE}\tkraken2\t${DBVER}\tconfidence=0.1" >> tax/provenance.tsv
```
Validate: `awk '$4=="U"{u=$1} $4=="R"{r=$1} END{printf "classified=%.1f%%\n", 100*r/(r+u)}' tax/${SAMPLE}.kreport`; Bracken file has `fraction_total_reads` summing to ≈1 (`awk 'NR>1{s+=$7} END{print s}' tax/${SAMPLE}.S.bracken`). Memory: index size from the Langmead table is the RAM floor (Standard 2026-10 ≈109 GB; Standard-16 ≈15 GB; PlusPF ≈116 GB); `--memory-mapping` trades speed for RAM.

3. MetaPhlAn 4 (marker-based) and sylph (sketch-based), either as the primary profile or as cross-check.

```bash
set -euo pipefail
SAMPLE=S1
THREADS=16
MPADB=db/metaphlan
metaphlan "clean/${SAMPLE}_R1.fastq.gz,clean/${SAMPLE}_R2.fastq.gz" --input_type fastq --nproc "${THREADS}" --bowtie2db "${MPADB}" --bowtie2out "tax/${SAMPLE}.mpa.bz2" --unclassified_estimation -o "tax/${SAMPLE}.mpa_profile.tsv"
head -1 "tax/${SAMPLE}.mpa_profile.tsv" >> tax/provenance.tsv
```
Validate: first line is `#mpa_vJan26_CHOCOPhlAnSGB_...` (or the version you installed with `metaphlan --install --index <name>`); `grep -c 't__SGB' tax/${SAMPLE}.mpa_profile.tsv` > 0. SGB names: `t__SGB1234` = known species-level genome bin; `_group` suffix = merged SGBs; `UNCLASSIFIED` row only with `--unclassified_estimation`.

```bash
set -euo pipefail
SAMPLE=S1
THREADS=16
SYLDB=db/gtdb-r232-c200-dbv1.syldb
sylph sketch -1 "clean/${SAMPLE}_R1.fastq.gz" -2 "clean/${SAMPLE}_R2.fastq.gz" -d sketches -t "${THREADS}"
sylph profile "${SYLDB}" "sketches/${SAMPLE}_R1.fastq.gz.paired.sylsp" -t "${THREADS}" > "tax/${SAMPLE}.sylph.tsv"
```
Validate: `cut -f1-5,8 tax/${SAMPLE}.sylph.tsv | head` (columns include `Taxonomic_abundance`, `Sequence_abundance`, `Adjusted_ANI`, `Eff_cov`); convert to taxonomy strings with `sylph-tax taxprof`.

4. Merge profiles across samples (one tool, one database build per table).

```bash
set -euo pipefail
combine_bracken_outputs.py --files tax/*.S.bracken -o tax/bracken_species_merged.tsv
merge_metaphlan_tables.py tax/*.mpa_profile.tsv > tax/metaphlan_merged.tsv
```
Validate: `awk -F'\t' 'NR==1{print NF-1" samples"}' tax/bracken_species_merged.tsv`; confirm every input had the same `DBVER` in `tax/provenance.tsv` (`cut -f3 tax/provenance.tsv | sort -u | wc -l` → 1 per tool).

5. Assembly and binning (per sample or per `group` co-assembly).

```bash
set -euo pipefail
SAMPLE=S1
THREADS=32
megahit -1 "clean/${SAMPLE}_R1.fastq.gz" -2 "clean/${SAMPLE}_R2.fastq.gz" -t "${THREADS}" --min-contig-len 1000 -o "asm/${SAMPLE}"
bowtie2-build --threads "${THREADS}" "asm/${SAMPLE}/final.contigs.fa" "asm/${SAMPLE}/contigs"
bowtie2 -p "${THREADS}" -x "asm/${SAMPLE}/contigs" -1 "clean/${SAMPLE}_R1.fastq.gz" -2 "clean/${SAMPLE}_R2.fastq.gz" | samtools sort -@ 8 -o "asm/${SAMPLE}/reads.bam" - && samtools quickcheck "asm/${SAMPLE}/reads.bam"
samtools index "asm/${SAMPLE}/reads.bam"
jgi_summarize_bam_contig_depths --outputDepth "asm/${SAMPLE}/depth.txt" "asm/${SAMPLE}/reads.bam"
metabat2 -i "asm/${SAMPLE}/final.contigs.fa" -a "asm/${SAMPLE}/depth.txt" -o "bins/${SAMPLE}/metabat2/bin" -t "${THREADS}" -m 1500 --seed 1
SemiBin2 single_easy_bin -i "asm/${SAMPLE}/final.contigs.fa" -b "asm/${SAMPLE}/reads.bam" -o "bins/${SAMPLE}/semibin2" --environment human_gut -t "${THREADS}"
run_MaxBin.pl -contig "asm/${SAMPLE}/final.contigs.fa" -abund "asm/${SAMPLE}/depth.txt" -out "bins/${SAMPLE}/maxbin2/bin" -thread "${THREADS}"
```
Validate: `seqkit stats asm/${SAMPLE}/final.contigs.fa` (N50, total length); `samtools flagstat asm/${SAMPLE}/reads.bam | grep 'mapped ('` (≥70 % reads assembled for gut; environmental often 30–50 %); `ls bins/${SAMPLE}/*/ | wc -l`. Use `--environment` matching the biome or omit it for SemiBin2 self-training; metaSPAdes alternative: `spades.py --meta -1 R1 -2 R2 -t ${THREADS} -m 500 -o asm/${SAMPLE}_spades`.

6. Bin refinement and QC.

```bash
set -euo pipefail
SAMPLE=S1
THREADS=32
Fasta_to_Contig2Bin.sh -i "bins/${SAMPLE}/metabat2" -e fa > "bins/${SAMPLE}/metabat2.c2b.tsv"
Fasta_to_Contig2Bin.sh -i "bins/${SAMPLE}/semibin2/output_bins" -e fa.gz > "bins/${SAMPLE}/semibin2.c2b.tsv"
Fasta_to_Contig2Bin.sh -i "bins/${SAMPLE}/maxbin2" -e fasta > "bins/${SAMPLE}/maxbin2.c2b.tsv"
DAS_Tool -i "bins/${SAMPLE}/metabat2.c2b.tsv,bins/${SAMPLE}/semibin2.c2b.tsv,bins/${SAMPLE}/maxbin2.c2b.tsv" -l metabat2,semibin2,maxbin2 -c "asm/${SAMPLE}/final.contigs.fa" -o "bins/${SAMPLE}/dastool/${SAMPLE}" --write_bins -t "${THREADS}" --score_threshold 0.5
checkm2 predict --threads "${THREADS}" --input "bins/${SAMPLE}/dastool/${SAMPLE}_DASTool_bins" -x fa --output-directory "qc/checkm2_${SAMPLE}"
gtdbtk classify_wf --genome_dir "bins/${SAMPLE}/dastool/${SAMPLE}_DASTool_bins" -x fa --out_dir "qc/gtdbtk_${SAMPLE}" --cpus "${THREADS}" --mash_db db/gtdbtk_mash
```
Validate: `awk -F'\t' 'NR>1 && $2>=90 && $3<=5' qc/checkm2_${SAMPLE}/quality_report.tsv | wc -l` (HQ count); `cut -f1,2 qc/gtdbtk_${SAMPLE}/gtdbtk.bac120.summary.tsv | head`; `gtdbtk check_install`. GTDB-Tk 2.8.0 uses reference data R232 (`GTDBTK_DATA_PATH`; ~110 GB RAM for pplacer unless `--pplacer_cpus 1`/`--skip_ani_screen` trade-offs are accepted). Dereplicate MAGs across samples with `dRep dereplicate -pa 0.9 -sa 0.95` before abundance estimation.

7. Functional and AMR annotation.

```bash
set -euo pipefail
SAMPLE=S1
THREADS=16
humann --input "clean/${SAMPLE}_cat.fastq.gz" --output "func/${SAMPLE}" --threads "${THREADS}" --taxonomic-profile "tax/${SAMPLE}.mpa_profile.tsv" --nucleotide-database db/chocophlan --protein-database db/uniref90
humann_renorm_table --input "func/${SAMPLE}/${SAMPLE}_cat_pathabundance.tsv" --output "func/${SAMPLE}/pathabundance_cpm.tsv" --units cpm
amrfinder -u
amrfinder --nucleotide "asm/${SAMPLE}/final.contigs.fa" --plus --threads "${THREADS}" -o "amr/${SAMPLE}.amrfinder.tsv"
rgi main --input_sequence "asm/${SAMPLE}/final.contigs.fa" --output_file "amr/${SAMPLE}.rgi" --input_type contig --num_threads "${THREADS}" --clean
```
Validate: `grep -c 'UNMAPPED\|UNINTEGRATED' func/${SAMPLE}/*_pathabundance.tsv` (both present); `amrfinder --database_version`; `wc -l amr/${SAMPLE}.amrfinder.tsv amr/${SAMPLE}.rgi.txt`. HUMAnN input is the concatenated pair (`cat R1 R2 > _cat`); for MAG proteins: `pyrodigal -i bin.fa -a bin.faa -p meta` then `emapper.py -i bin.faa -o bin --cpu ${THREADS} --data_dir db/eggnog`.

8. 16S amplicons with QIIME 2 (V4 515F/806R shown; adjust trunc lengths from the demux quality plot, never by default).

```bash
set -euo pipefail
THREADS=16
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.tsv --input-format PairedEndFastqManifestPhred33V2 --output-path amp/demux.qza
qiime demux summarize --i-data amp/demux.qza --o-visualization amp/demux.qzv
qiime cutadapt trim-paired --i-demultiplexed-sequences amp/demux.qza --p-front-f GTGYCAGCMGCCGCGGTAA --p-front-r GGACTACNVGGGTWTCTAAT --p-discard-untrimmed --p-cores "${THREADS}" --o-trimmed-sequences amp/trimmed.qza
qiime dada2 denoise-paired --i-demultiplexed-seqs amp/trimmed.qza --p-trunc-len-f 220 --p-trunc-len-r 180 --p-n-threads "${THREADS}" --o-table amp/table.qza --o-representative-sequences amp/rep-seqs.qza --o-denoising-stats amp/dada2-stats.qza
qiime feature-classifier classify-sklearn --i-classifier db/silva-144-99-515-806-classifier.qza --i-reads amp/rep-seqs.qza --p-n-jobs "${THREADS}" --o-classification amp/taxonomy.qza
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences amp/rep-seqs.qza --o-alignment amp/aln.qza --o-masked-alignment amp/masked.qza --o-tree amp/unrooted.qza --o-rooted-tree amp/rooted.qza
qiime diversity core-metrics-phylogenetic --i-phylogeny amp/rooted.qza --i-table amp/table.qza --p-sampling-depth 10000 --m-metadata-file metadata.tsv --output-dir amp/core
```
Validate: `qiime metadata tabulate --m-input-file amp/dada2-stats.qza --o-visualization amp/dada2-stats.qzv` — `percentage of input non-chimeric` ≥50 % per sample; choose `--p-sampling-depth` from `qiime feature-table summarize` so that ≥90 % of real samples are retained (blanks will drop — expected). Classifier: download the SILVA classifier matching your QIIME 2 release and region from the QIIME 2 data resources page (SILVA 138.2 and 144 available; `<check current: https://library.qiime2.org/data-resources>`); for ITS use UNITE with `--p-trunc-len 0` and ITSxpress.

9. Diversity, contamination and differential abundance (R).

```r
suppressPackageStartupMessages({library(phyloseq); library(decontam); library(vegan); library(ANCOMBC)})
ps <- qiime2R::qza_to_phyloseq(features="amp/table.qza", tree="amp/rooted.qza", taxonomy="amp/taxonomy.qza", metadata="metadata.tsv")
sample_data(ps)$is_neg <- sample_data(ps)$sample_type == "blank"
contam <- isContaminant(ps, method="combined", neg="is_neg", conc="dna_conc_ng_ul", threshold=0.1)
ps_clean <- prune_taxa(!contam$contaminant, ps)
ps_clean <- subset_samples(ps_clean, !is_neg)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
clr <- microbiome::transform(ps_clean, "clr")
ait <- vegan::vegdist(t(otu_table(clr)), method="euclidean")
perm <- vegan::adonis2(ait ~ group + batch, data=data.frame(sample_data(ps_clean)), permutations=999, by="margin")
da <- ancombc2(data=ps_clean, tax_level="Genus", fix_formula="group + batch", p_adj_method="BH", prv_cut=0.10, lib_cut=1000, struc_zero=TRUE)
write.csv(da$res, "stats/ancombc2_genus.csv", row.names=FALSE)
print(perm); print(table(contam$contaminant)); sessionInfo()
```
Validate: `table(contam$contaminant)` non-zero when blanks exist; `betadisper` p>0.05 before trusting PERMANOVA group effect; read-depth not correlated with PC1 (`cor.test(sample_sums(ps_clean), scores)`).

## QC thresholds

| Metric | Default | Rationale / when to deviate |
| --- | --- | --- |
| Reads per shotgun sample after host removal | ≥5 M pairs for profiling; ≥20 M for assembly/MAGs | Below 1 M, species-level profiles are unstable; low-biomass samples can't be rescued by depth. |
| Host reads remaining | <0.1 % human after hostile/T2T | Residual human reads classify as spurious bacteria and are a privacy risk. |
| Kraken2 classified % | Gut 60–90 % (Standard); soil/marine 10–40 % | Low on environmental is expected; low on gut → wrong DB or heavy host. |
| Kraken2 `--confidence` | 0.1 | Lower for short reads/low depth; 0.2–0.5 when precision (clinical) matters; always record. |
| Bracken threshold `-t` | 10 reads | Removes singleton species; raise with depth (e.g. 50 at >20 M). |
| MetaPhlAn unclassified fraction | Gut <30 %; other sites higher | >50 % in gut → non-human-associated DB mismatch or contamination. |
| Assembly N50 / mapped reads | N50 ≥5 kb (gut); ≥70 % reads map back | Complex environmental samples legitimately fall below. |
| Bin quality (CheckM2) | HQ ≥90/≤5; MQ ≥50/≤10; report strain heterogeneity | MIMAG; CheckM2 completeness is model-based — confirm rRNA/tRNA with barrnap/tRNAscan for HQ claims. |
| GTDB-Tk | `classification` present; `warnings` empty; `fastani_ani` ≥95 % for species | Novel species: genus-only call, name as `g__X sp.` — not a failure. |
| DADA2 retained non-chimeric | ≥50 % of input | Lower → trunc too long (quality) or primers not removed. |
| Rarefaction depth | ≥10 000 reads (16S) retaining ≥90 % of real samples | Report dropped samples; for shotgun alpha use ≥1 M reads or coverage-based (iNEXT). |
| Negative controls | ≥1 extraction blank + 1 PCR NTC per batch; blanks ≤5 % of mean sample reads | More reads in blanks than samples = low-biomass study; conclusions limited to taxa absent from blanks. |
| Mock community | Expected taxa recovered, false positives <1 % relative abundance | Zymo/ATCC mock per run validates DB + pipeline. |

## Pitfalls

| Symptom | Cause | Fix |
| --- | --- | --- |
| Species counts differ 2× between runs of the "same" pipeline | Database build drift (quarterly Kraken2 indexes, MetaPhlAn vJan25→vJan26, GTDB r226→r232) | Pin the index path/date in `provenance.tsv`; reprofile all samples when upgrading; never merge tables across builds. |
| Kraken2 reports many species at low abundance | Raw Kraken2 counts used as abundance | Always run Bracken; apply `-t`; use `--confidence`. |
| Relative abundances of taxon A "decrease" when B blooms | Compositional closure | Use CLR/ALR, ANCOM-BC2/ALDEx2; or absolute quantification (spike-ins, qPCR). |
| Skin/tissue/placenta profiles dominated by `Ralstonia`, `Burkholderia`, `Bradyrhizobium`, `Cutibacterium` | Reagent/kit contamination in low-biomass samples | Sequence blanks, run `decontam`, report kit-ome; don't interpret without controls. |
| Taxa present in every sample at ~0.01 % including blanks | Index hopping / cross-talk on patterned flow cells | Use unique dual indexes; subtract blank prevalence; apply per-sample minimum relative abundance (e.g. 0.01–0.1 %). |
| Human DNA still ~1 % after GRCh38 Bowtie2 | Unplaced/centromeric/HLA human sequence absent from GRCh38 | Add T2T-CHM13 (hostile default) and HLA; second pass with Kraken2 `Homo sapiens` extraction. |
| MetaPhlAn fails with database/version mismatch | `--index` from a different MetaPhlAn major version (vJan21 vs vJan25/26 SGB) | `metaphlan --install --index <matching>`; HUMAnN 3.9 requires MetaPhlAn ≤4.1 profiles or `--bypass-nucleotide-index` with mapped taxonomy. |
| Kraken2 killed (OOM) | Index size > RAM | Use Standard-16/PlusPF-16 capped indexes, `--memory-mapping`, or a bigger node; RAM floor = index size column. |
| Many bins with >10 % contamination | Single binner, co-assembly of divergent samples, repeats | DAS Tool/MAGScoT refinement; bin per sample with multi-sample depth (`jgi_summarize_bam_contig_depths` over all BAMs); `--binning_map_mode group` in nf-core/mag. |
| 16S DADA2 yields thousands of singleton ASVs | Primers not trimmed, or trunc lengths cut into low-quality tails | `cutadapt --p-discard-untrimmed`; pick trunc at Q≥30 median; merge requires ≥12 bp overlap after truncation. |
| V3–V4 taxonomy stops at family | Classifier trained on wrong region or full-length un-weighted | Train/download region-specific SILVA classifier for QIIME 2 release; check `--p-confidence 0.7`. |
| PERMANOVA significant but dispersion differs | Heteroscedastic groups | `betadisper` + `permutest`; report both. |
| GTDB-Tk and NCBI names disagree | Different taxonomies (GTDB splits/renames, e.g. `Escherichia flexneri`) | Report GTDB names with release (R232) and map to NCBI via `gtdb_to_ncbi_majority_vote.py` when needed. |

## nf-core

All three exist in `pipelines.json` (checked 2026-11-22). Nextflow/profile/`-resume` mechanics: `.claude/skills/genomics-workflow/SKILL.md`.

**nf-core/taxprofiler 2.0.1** — read-based profiling (Kraken2, Bracken, MetaPhlAn, sylph, Centrifuge, KrakenUniq, KMCP, ganon, mOTUs, MALT, metacache, melon), host removal (`--perform_shortread_hostremoval --hostremoval_reference ref.fa`), standardised tables (`--run_profile_standardisation`).

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
S1,RUN1,ILLUMINA,clean/S1_R1.fastq.gz,clean/S1_R2.fastq.gz,
```

```csv
tool,db_name,db_params,db_type,db_path
kraken2,k2std_20261022,--confidence 0.1,short,db/k2_standard_20261022.tar.gz
bracken,k2std_20261022,--confidence 0.1;-r 150,short,db/k2_standard_20261022.tar.gz
metaphlan,mpa_vJan26,,short,db/metaphlan
sylph,gtdb_r232,,short,db/gtdb-r232-c200-dbv1.syldb
```
`db_params` for bracken is `kraken2_args;bracken_args`. Run: `nextflow run nf-core/taxprofiler -r 2.0.1 -profile singularity --input samplesheet.csv --databases databases.csv --outdir results --perform_shortread_qc --perform_shortread_hostremoval --hostremoval_reference ref/chm13v2.0.fa --run_kraken2 --run_bracken --run_metaphlan --run_sylph --run_profile_standardisation`. Validate: `ls results/taxpasta/*.tsv results/multiqc/multiqc_report.html`.

**nf-core/mag 5.5.0** — assembly (MEGAHIT/metaSPAdes), binning (MetaBAT2, MaxBin2, CONCOCT, SemiBin2 via `--binning_map_mode`), refinement (`--refine_bins_dastool`), QC (`--run_checkm2`, BUSCO), taxonomy (`--gtdb_db`), optional host removal (`--host_fasta`/`--host_genome`).

```csv
sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform
S1,g1,clean/S1_R1.fastq.gz,clean/S1_R2.fastq.gz,,ILLUMINA,
```
Minimal header `sample,group,short_reads_1,short_reads_2` is accepted; optional `run` column enables lane merging. Run: `nextflow run nf-core/mag -r 5.5.0 -profile singularity --input samplesheet.csv --outdir results --host_fasta ref/chm13v2.0.fa --refine_bins_dastool --run_checkm2 --gtdb_db db/gtdbtk_r232.tar.gz`. Validate: `ls results/GenomeBinning/QC/checkm2_summary.tsv results/Taxonomy/GTDB-Tk/*/gtdbtk_summary.tsv`.

**nf-core/ampliseq 3.1.0** — cutadapt + DADA2 + taxonomy (default `--dada_ref_taxonomy sbdi-gtdb=R11-RS232-1`; `silva=144`, `silva=138.2`, `gtdb=R11-RS232`, `greengenes2=2024.09`, `unite-fungi=10.0` available) + QIIME 2 diversity + ANCOM-BC. Samplesheet is TSV (`assets/samplesheet.tsv`):

```text
sample	fastq_1	fastq_2	run
S1	raw/S1_R1.fastq.gz	raw/S1_R2.fastq.gz	1
```
Run: `nextflow run nf-core/ampliseq -r 3.1.0 -profile singularity --input samplesheet.tsv --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT --metadata metadata.tsv --dada_ref_taxonomy silva=144 --outdir results`. Validate: `ls results/dada2/ASV_table.tsv results/qiime2/diversity/`.

Related pipelines that exist: `nf-core/funcscan` 4.0.0 (AMR/BGC/AMP screening on contigs), `nf-core/createtaxdb` 3.2.0 (builds matching Kraken2/Bracken/MetaPhlAn-free DBs from one genome set), `nf-core/metatdenovo` (metatranscriptome assembly).

## Validation

- Mock community (ZymoBIOMICS D6300/D6331, ATCC MSA-1003) per batch: expected 8–20 taxa recovered within 2× of nominal abundance; count false positives above 0.1 %.
- Negative controls: blank read counts and taxa reported in the methods; `decontam` output retained.
- Cross-tool concordance: Kraken2/Bracken vs MetaPhlAn vs sylph species lists overlap (Jaccard on taxa ≥0.5 %) — disagreements flag DB gaps, not biology.
- Provenance table per run: tool version, database name/date/MD5, `--confidence`, read length used for Bracken, host reference, classifier artifact UUID (`qiime tools peek`).
- MAGs: CheckM2 + GTDB-Tk tables, dRep clusters, MIMAG class per bin; CAMI benchmark datasets or synthetic communities for pipeline comparison.
- Statistics: `sessionInfo()`, permutation seeds, rarefaction depth and dropped samples recorded.
- Benchmark mechanics (truth sets, scoring scripts): `.claude/skills/genomics-benchmark/SKILL.md`.

## Additional resources

- `.claude/skills/genomics-qc/SKILL.md`: fastp/FastQC trimming and complexity filtering before profiling.
- `.claude/skills/genomics-longread/SKILL.md`: ONT/HiFi basecalling, metaFlye/hifiasm-meta, polishing for long-read metagenomes.
- `.claude/skills/genomics-workflow/SKILL.md`: Nextflow/Snakemake execution, profiles, resource configs for taxprofiler/mag/ampliseq.
- `.claude/skills/genomics/reference.md`: downloading T2T-CHM13, GRCh38, Kraken2/GTDB/SILVA database archives and recording checksums.
- `.claude/skills/genomics-benchmark/SKILL.md`: scoring profiles/MAGs against CAMI or mock-community truth.
- Database index tables: Kraken2/Bracken <https://benlangmead.github.io/aws-indexes/k2> (Standard/PlusPF/PlusPFP/GTDB, dated 2026-10-22 at time of writing); MetaPhlAn DB index `<check current: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/>`; sylph pre-built DBs <https://sylph-docs.github.io/pre%E2%80%90built-databases/>; GTDB-Tk data (R232) <https://ecogenomics.github.io/GTDBTk/installing/index.html>; AMRFinderPlus DB `amrfinder -u` (latest observed 2026-10-29.1).
