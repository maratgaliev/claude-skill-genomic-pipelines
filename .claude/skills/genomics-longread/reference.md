# genomics-longread — reference

Supporting detail for `SKILL.md`. Sections are referenced as "§n" from the workflow. Versions verified 2026-11-22 are listed in `SKILL.md` → Environment.

## §1 ONT basecalling details (dorado 2.1.x)

### Model naming

`dna_<pore>_<enzyme>_<speed>bps_<fast|hac|sup>@v<X.Y.Z>[_<mod>@v<N>]`, e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v3`. The short "model complex" form `sup,5mCG_5hmCG` resolves to the newest compatible model for the pod5's chemistry/sample rate; pin it explicitly in production (`sup@v5.2.0,5mCG_5hmCG@v3`) and record it from the BAM `@RG DS:` line.

| Data | Current simplex model | Notes |
|---|---|---|
| R10.4.1 E8.2, 5 kHz (MinKNOW ≥ 23.x) | `hac@v6.0.0` (fastest, ≈Q18 median), `sup@v5.2.0` (≈Q20.6 median; no sup@v6 released) | Both ship 5mCG_5hmCG, 5mC_5hmC (all-context), 6mA, 4mC_5mC mod models |
| R10.4.1 E8.2, 4 kHz (2022–mid 2023) | `*@v4.1.0` family (hac/sup) | v4.2+ models reject 4 kHz data; Clair3 has matching `v410` models |
| R9.4.1 | `dna_r9.4.1_e8_sup@v3.6` (legacy; last R9 models) | Clair3 `r941_prom_sup_g5014`; do not mix with R10 models |
| RNA004 (direct RNA) | `rna004_130bps_sup@v5.2.0` (+`m6A`, `pseU`, `inosine_m6A`, `m5C`) | Reads are reverse-complement of the mRNA; dorado emits sense orientation; use `--estimate-poly-a` |
| FLO-PRO114P (PromethION Plus) | models only in dorado 1.3.3 (MinKNOW 26.01 LTS) | Not yet in 2.x at time of writing |

### Useful flags

| Flag | Purpose |
|---|---|
| `--models-directory DIR` / `DORADO_MODELS_DIRECTORY` | Cache models; otherwise re-downloaded and deleted per run |
| `--emit-moves` | Adds `mv:B:c` move table (needed for Clair3 `*_with_mv` + `--enable_dwell_time`, `dorado polish`, signal tools); ~2× BAM size |
| `--min-qscore N` | Drops reads below mean Q (10 is a sane default; use 7 for metagenomics yield) |
| `--kit-name SQK-NBD114-24` | Inline barcode classification into `BC` tag/read groups; then `dorado demux --no-classify` |
| `--barcode-both-ends` | Require both barcodes; lowers mis-assignment for pooled clinical samples |
| `--no-trim` | Keep adapters/barcodes (required if you will `demux` later from the BAM) |
| `--reference ref.mmi` | Basecall + align in one pass (minimap2 inside dorado; `--mm2-opts "-x lr:hq --MD"`) |
| `--estimate-poly-a` | Poly(A)/(T) tail length → `pt:i` tag (cDNA/dRNA) |
| `--resume-from previous.bam` | Resume an interrupted basecall |
| `--device cuda:0,1` / `--batchsize 0` | GPU selection / auto batch; CPU basecalling of sup is impractical for WGS |

### Duplex

```bash
set -euo pipefail
POD5_DIR=pod5/
SAMPLE=S1
dorado duplex sup,5mCG_5hmCG "${POD5_DIR}" --models-directory /opt/dorado-models > "${SAMPLE}.duplex.bam"
samtools view -c -e '[dx]==1' "${SAMPLE}.duplex.bam"
```

Output mixes simplex (`dx:i:0`), simplex parents of duplex pairs (`dx:i:-1`) and duplex reads (`dx:i:1`, ~Q30). For variant/methylation calling drop `dx:i:-1` records to avoid double counting. Duplex rate depends on library prep (ligation kit, duplex-optimised loading); 10–40 % of bases is typical, <5 % means the run was not duplex-oriented and simplex sup is the better use of GPU time.

### Demultiplexing

```bash
set -euo pipefail
SAMPLE_RUN=run1
dorado demux --kit-name SQK-NBD114-24 --emit-summary --output-dir demux_"${SAMPLE_RUN}" "${SAMPLE_RUN}.unaligned.bam"
ls demux_"${SAMPLE_RUN}"/
```

Already-classified BAM (basecalled with `--kit-name`): `dorado demux --no-classify --output-dir out/ in.bam`. Aligned input needs `--no-trim` (optionally `--sort-bam`) or alignments are discarded. Sample sheets: `--sample-sheet sheet.csv` (MinKNOW format: `flow_cell_id,kit,barcode,alias`).

### pod5 / fast5

```bash
set -euo pipefail
pod5 convert fast5 fast5/*.fast5 --output pod5/ --one-to-one fast5/ --threads 8
pod5 inspect summary pod5/*.pod5 | head
pod5 view pod5/ --include "read_id,channel,end_reason" --output reads.tsv
```

`pod5 subset` splits by channel/time for re-basecalling subsets; `pod5 merge` consolidates small files (dorado prefers few large files). Guppy and fast5 are legacy; re-basecalling old fast5 with current dorado models is worth it whenever variant or methylation accuracy matters (R9 data still needs R9 models).

### Other dorado subcommands

`dorado aligner ref.mmi in.bam > aligned.bam` (keeps all tags, uses `lr:hq` defaults); `dorado summary in.bam > summary.tsv` (sequencing-summary-like table for pycoQC-style QC); `dorado trim`; `dorado correct` (HERRO read correction for assembly); `dorado polish` (assembly polishing); `dorado smallvar` (ONT's bundled small-variant caller, gVCF since 2.1.0 — treat as emerging next to Clair3/DeepVariant).

## §2 Small-variant model matrix

| Data (from `@RG DS:` / `@PG`) | Clair3 2.0.x `--platform` / `--model_path` | DeepVariant 1.10 `--model_type` |
|---|---|---|
| ONT R10.4.1 5 kHz, dorado sup v5.2.0 | `ont` / `r1041_e82_400bps_sup_v520` (`_with_mv` + `--enable_dwell_time` if `mv` tags) | `ONT_R104` |
| ONT R10.4.1 5 kHz, dorado hac v6.0.0 | `ont` / `r1041_e82_400bps_hac_v600` (`_with_mv` variant exists) | `ONT_R104` |
| ONT R10.4.1 5 kHz, dorado v5.0.0 / v4.3.0 | `ont` / `r1041_e82_400bps_{sup,hac}_v500` or `_v430` | `ONT_R104` |
| ONT R10.4.1 4 kHz, dorado/guppy v4.1.0 | `ont` / `r1041_e82_400bps_{sup,hac}_v410` | `ONT_R104` |
| ONT R9.4.1 Guppy5/6 sup | `ont` / `r941_prom_sup_g5014` | not supported (R9 model removed) |
| PacBio HiFi Revio | `hifi` / `hifi_revio` | `PACBIO` |
| PacBio HiFi Sequel II | `hifi` / `hifi_sequel2` | `PACBIO` |
| Non-human / draft reference | add `--include_all_ctgs` (Clair3 defaults to chr1–22,X,Y) | `--regions` as needed |
| Haploid organism | `--haploid_precise` (1/1 only) or `--haploid_sensitive` | `--haploid_contigs <all contigs>` |
| XY sample | `--gender male` (≥2.0.3; handles chrX/Y ploidy) | `--haploid_contigs chrX,chrY --par_regions_bed GRCh38_PAR.bed` |
| Long-read RNA | Clair3-RNA (separate repo) | `--model_type` RNA-seq model (1.10, short-read trained; check case study) |

Clair3 outputs `merge_output.vcf.gz` (+ `merge_output.gvcf.gz` with `--gvcf`, `phased_merge_output.vcf.gz` with `--enable_phasing`). Clair3 ≥ 2.0 runs on CUDA/ROCm/Apple GPUs (`--use_gpu`, `--device cuda:0`), ~5× faster than CPU. DeepVariant 1.10 emits natively phased long-read VCF/gVCF (`PS` tags) — WhatsHap is still needed for haplotagging BAMs.

Somatic small variants (ClairS, ClairS-TO, DeepSomatic long-read models): `.claude/skills/genomics-somatic/SKILL.md`.

## §3 Structural variants

### cuteSV with platform parameters (second caller for consensus)

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
PLATFORM=ont
mkdir -p cutesv_"${SAMPLE}"
if [ "${PLATFORM}" = hifi ]; then
  cuteSV --threads 16 --genotype --sample "${SAMPLE}" --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 "${SAMPLE}.bam" "${REF}" "${SAMPLE}.cutesv.vcf" cutesv_"${SAMPLE}"
else
  cuteSV --threads 16 --genotype --sample "${SAMPLE}" --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 "${SAMPLE}.bam" "${REF}" "${SAMPLE}.cutesv.vcf" cutesv_"${SAMPLE}"
fi
bgzip -f "${SAMPLE}.cutesv.vcf" && bcftools index -f "${SAMPLE}.cutesv.vcf.gz"
```

`--min_support` defaults to 10; lower to 3–5 at 20× or for mosaic work. Force-calling/regenotyping moved to `cuteFC` (cuteSV ≥ 2.1.2). Default `--min_mapq 20`, `--min_size 30`.

### pbsv (PacBio legacy/SMRT Link path)

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
TRF=human_GRCh38_no_alt_analysis_set.trf.bed
pbsv discover --hifi --tandem-repeats "${TRF}" "${SAMPLE}.bam" "${SAMPLE}.svsig.gz"
pbsv call --hifi -j 16 "${REF}" "${SAMPLE}.svsig.gz" "${SAMPLE}.pbsv.vcf"
bgzip -f "${SAMPLE}.pbsv.vcf" && bcftools index -f "${SAMPLE}.pbsv.vcf.gz"
```

`--hifi` (alias of `--ccs`) relaxes thresholds to `-A 3 -B 2 -O 3 -S 0 -P 10`. Requires `SM` in `@RG` (set via `pbmm2 --sample`). Joint calling: pass all samples' `.svsig.gz` to one `pbsv call`. pbsv deprecated CNV calling in 2.8; sawfish (SV+CNV from local assembly) is PacBio's current recommendation and HiPhase/pacvar default.

### sawfish joint calling (trio / cohort)

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
for S in HG002 HG003 HG004; do sawfish discover --threads 16 --ref "${REF}" --bam "${S}.bam" --output-dir sawfish_"${S}"; done
sawfish joint-call --threads 16 --sample sawfish_HG002 --sample sawfish_HG003 --sample sawfish_HG004 --output-dir sawfish_trio
bcftools view -H sawfish_trio/genotyped.sv.vcf.gz | wc -l
```

Add `--expected-cn data/expected_cn/expected_cn.hg38.XY.bed --cnv-excluded-regions data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz` to `discover` for usable CNV output (files ship with the sawfish distribution).

### Sniffles2 population / trio merge

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
TRF=human_GRCh38_no_alt_analysis_set.trf.bed
for S in HG002 HG003 HG004; do sniffles --threads 16 --reference "${REF}" --tandem-repeats "${TRF}" --input "${S}.bam" --snf "${S}.snf"; done
sniffles --threads 16 --reference "${REF}" --input HG002.snf HG003.snf HG004.snf --vcf trio.sniffles.vcf.gz
bcftools query -l trio.sniffles.vcf.gz
```

A TSV (`path<TAB>sample_id`) can replace the `.snf` list. `--mosaic` enables low-VAF (≥5 %) calling for somatic/mosaic SVs; `--genotype-vcf known.vcf` force-genotypes a published set. Mendelian-violation rate in a trio should be < 5 % of PASS SVs.

### Somatic SV (emerging; expect method churn)

| Tool | Platform | Command skeleton | Notes |
|---|---|---|---|
| Severus 1.7 | ONT + HiFi, T/N or tumour-only | `severus --target-bam T.haplotagged.bam --control-bam N.haplotagged.bam --phasing-vcf N.phased.vcf.gz --vntr-bed ${TRF} --out-dir severus_out -t 16` | Haplotype-aware; `--PON PoN_1000G_hg38.tsv.gz` for tumour-only; outputs `somatic_SVs/severus_somatic.vcf` + breakpoint graphs |
| SAVANA 1.3.8 | ONT (`--ont`, default) / HiFi (`--pb`) | `savana --tumour T.bam --normal N.bam --ref ${REF} --outdir savana_out --threads 16 --contigs contigs.chr.hg38.txt` | Classifier trained on COLO829; add `--snp_vcf` for allele-specific CN |
| Sniffles2 `--mosaic` | both | `sniffles --mosaic --input T.bam --vcf T.mosaic.vcf.gz` | Not T/N aware; subtract normal with `bcftools isec`/Truvari |

Mark somatic long-read SV calls as research-grade; validate on COLO829 (ONT/HiFi truth set) before using on patient data.

### SV post-processing

Multi-caller merge: Jasmine or `truvari collapse` after `bcftools sort | bcftools norm`; keep `SUPP_VEC`. Annotate with AnnotSV/SVAnnotate — see `.claude/skills/genomics-annotation/SKILL.md`. Strip `BND` on decoys and `chrM`; require `PASS`, `SVLEN ≥ 50`, `RE/DV ≥ 3` (ONT) or `≥ 2` (HiFi, sawfish genotypes already filtered).

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## §4 Methylation detail

### modkit options worth knowing

| Option | Use |
|---|---|
| `--cpg --combine-strands` | CpG-only, strand-collapsed bedMethyl (one row per CpG) |
| `--phased --prefix S` | Haplotype split (`S_hp1.bedmethyl`, `S_hp2.bedmethyl`, `S_combined.bedmethyl`) from `HP`-tagged BAM |
| `--combine-mods` | Report 5mC+5hmC as "modified" (matches bisulfite semantics) |
| `--filter-threshold 0.8` / `--mod-thresholds m:0.8` | Override the auto 10th-percentile confidence filter; record it |
| `--motif GATC 1` / `--motif A 0` | Bacterial 6mA/4mC contexts; `--modified-bases 6mA` with `--organism` presets |
| `--bgzf` | Write bgzipped output (then `tabix -p bed`) — required for `modkit dmr` |
| `--region chr20` / `--include-bed` | Subset for fast QC |
| `modkit extract calls` | Per-read, per-site calls (for allele-specific or single-molecule analyses) |
| `modkit pileup-hemi` | Paired CpG (hemi-methylation) from duplex reads |
| `modkit bedmethyl tobigwig` | Browser tracks |

### DMR

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
REGIONS=cpg_islands.bed
for S in case ctrl; do bgzip -kf "${S}.bedmethyl" && tabix -f -p bed "${S}.bedmethyl.gz"; done
modkit dmr pair -a ctrl.bedmethyl.gz -b case.bedmethyl.gz -r "${REGIONS}" --ref "${REF}" --base C --threads 16 -o case_vs_ctrl.cpgi.dmr.bed --log-filepath dmr.log
modkit dmr pair -a ctrl.bedmethyl.gz -b case.bedmethyl.gz --ref "${REF}" --base C --threads 16 --segment case_vs_ctrl.segments.bed -o case_vs_ctrl.sites.bed
sort -k5,5gr case_vs_ctrl.cpgi.dmr.bed | head
```

`dmr multi` handles > 2 samples; `dmr pair` without `-r` runs single-site mode (add `--segment` for HMM segmentation; 0.6.x reports Cohen's h and effect size). Replicate-aware statistics (DSS, methylKit, NanoMethViz) need bedMethyl → count matrix; treat ONT/HiFi CpG counts like WGBS counts with ≥ 10× per site and follow `.claude/skills/genomics-methylation/SKILL.md` for the differential-methylation design.

### pb-CpG-tools

`aligned_bam_to_cpg_scores --bam S.bam --output-prefix S --threads 16 --pileup-mode model` (default `model`; `count` for raw proportions). Produces `S.combined.bed.gz` (+ `S.hap1/hap2` when BAM carries `HP`) and bigWigs. Requires on-instrument 5mC calling (Revio default) — `MM/ML` must be present and `pbmm2` preserves them; minimap2 needs `-y` after `samtools fastq -T MM,ML`.

### Modifications beyond CpG 5mC

| Target | Model/tool | Note |
|---|---|---|
| 5hmC | `5mCG_5hmCG` model → modkit code `h` | Low in most tissues except brain; validate against oxBS/ACE-seq |
| All-context 5mC (CHG/CHH, plants) | `5mC_5hmC` model; `modkit pileup --motif CHH 0` etc. | Do not use `--cpg` |
| 6mA / 4mC (bacteria, Fiber-seq) | `6mA`, `4mC_5mC` models; HiFi Fiber-seq via `fibertools` | Sniff motif enrichment with `modkit motif search` |
| RNA m6A, pseU, inosine, m5C | RNA004 `sup@v5.2.0` + `m6A,pseU` | Only direct RNA; calibrate per site with IVT control; dRNA m6A proportions are tissue/cell-type dependent |

## §5 Tandem repeats

| Item | ONT (Straglr 1.5.6) | HiFi (TRGT 5.1.0) |
|---|---|---|
| Targeted genotyping | `straglr.py S.bam ref.fa out --loci catalog.bed --genotype_in_size --min_support 2 --nprocs 16` | `trgt genotype --genome ref.fa --reads S.bam --repeats catalog.bed --karyotype XY --threads 16 --output-prefix S` (`--karyotype XX` default) |
| Genome-wide expansion scan | `straglr.py S.bam ref.fa scan --min_str_len 2 --max_str_len 100 --min_ins_size 100 --exclude hg38.exclude.bed --min_support 2 --max_num_clusters 2 --nprocs 32` | Catalog-driven only (Adotto ~1 M loci or Broad genome-wide catalog) |
| Catalog | `catalog.bed` 4-col: chrom start end motif (STRchive provides Straglr format) | STRchive (pathogenic), Adotto (genome-wide), Broad TR catalog; bundled `pathogenic_repeats.hg38.bed` is deprecated |
| Output | TSV (`genotype` col: size(reads);size(reads)) + BED + spanning reads | `S.vcf.gz` (`AL` allele length, `MC` motif counts, `AM` methylation) + `S.spanning.bam`; then `bcftools sort`, `samtools sort` + index |
| Visualise | IGV on spanning reads | `trgt plot --genome ref.fa --repeats catalog.bed --vcf S.sorted.vcf.gz --spanning-reads S.spanning.sorted.bam --repeat-id FMR1 --image FMR1.svg`; `--plot-type waterfall` for mosaicism |
| Cohort | concatenate TSVs; `straglr-compare` scripts | `trgt merge --vcf *.vcf.gz --genome ref.fa -o cohort.vcf.gz` or TDB |
| Targeted panels | adaptive sampling BED → `--loci` same BED | `--preset targeted` (PureTarget) raises `--max-depth` to 10000, uses `cluster` genotyper |

Pathogenic thresholds come from STRchive (`strchive.org/loci`); repeat size in the VCF is in reference coordinates and may differ from clinical PCR sizing by flanking bases — report motif copies and size both. Methylation at the FMR1/C9orf72 expansion is read directly from `AM` (TRGT) or from `modkit pileup --include-bed locus.bed`.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## §6 Assembly, long-read RNA, targeted sequencing

### De novo assembly (brief)

| Input | Command | Expectation (human) |
|---|---|---|
| HiFi 30× | `hifiasm -o S.asm -t 32 S.hifi.fq.gz` | contig N50 > 30 Mb, QV > 50; `gfatools gfa2fa S.asm.bp.p_ctg.gfa` |
| HiFi + ONT UL (≥ 100 kb N50, 30×) | `hifiasm -o S.asm -t 32 --ul S.ul.fq.gz S.hifi.fq.gz` | near-T2T; UL closes rDNA/centromere gaps |
| HiFi + Hi-C (phased) | `hifiasm -o S.asm -t 32 --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz [--ul ul.fq.gz] S.hifi.fq.gz` | `*.hic.hap1/hap2.p_ctg.gfa`; add `--telo-m CCCTAA --dual-scaf` for scaffolds |
| Trio (parental Illumina) | `yak count -k31 -b37 -t16 -o pat.yak pat_R*.fq.gz` ×2; `hifiasm -1 pat.yak -2 mat.yak ...` | best haplotype purity |
| ONT R10 sup only | `hifiasm --ont -o S.asm -t 32 S.ont.fq.gz` (≥ 0.21) or `flye --nano-hq S.ont.fq.gz --genome-size 3.1g --out-dir flye_S --threads 32` then `dorado polish` / medaka | Flye: faster, collapsed; hifiasm `--ont`: phased, needs ≥ 30× sup |
| ONT R9 / noisy | `flye --nano-raw ...` + 2× medaka | legacy |

QC: `quast.py -r chm13v2.0.fa -t 16 -o quast_S S.fa`; `merqury.sh S.meryl S.hap1.fa S.hap2.fa S` (k=21 meryl DB from HiFi or Illumina); `busco -i S.fa -m genome -l primates_odb12 -c 16 -o busco_S`; `paftools.js asmstat`; check mitochondrial contig and telomere counts (`seqtk telo`).

### Long-read RNA

| Platform / library | Pre-processing | Alignment | Quantification / isoform discovery |
|---|---|---|---|
| PacBio Iso-Seq / Kinnex (MAS-Seq) | `skera split movie.hifi_reads.bam mas16_primers.fasta movie.segmented.bam` (Kinnex only) → `lima movie.segmented.bam IsoSeq_v2_primers_12.fasta movie.fl.bam --isoseq --peek-guess` → `isoseq refine movie.fl.*.bam primers.fasta movie.flnc.bam --require-polya` → `isoseq cluster2 flnc.fofn clustered.bam` | `pbmm2 align --preset ISOSEQ --sort ref.fa clustered.bam clustered.aligned.bam` or `minimap2 -ax splice:hq -uf --MD -y` | `isoseq collapse clustered.aligned.bam collapsed.gff` → `pigeon classify`; or IsoQuant/bambu on FLNC alignments |
| ONT cDNA (PCS/PCB kits) | dorado `--kit-name SQK-PCB114-24`, `--estimate-poly-a`; `pychopper` for strand/full-length if not dorado-trimmed | `minimap2 -ax splice:hq -uf --junc-bed anno.bed12 -y --MD` (hac/sup v5+); `-ax splice -uf -k14` only for noisy/R9 | bambu (`bambu(reads, annotations, genome, ncore)`), IsoQuant (`isoquant --bam S.bam --reference ref.fa --genedb anno.gtf --data_type nanopore -o out`), oarfish (`oarfish -a S.transcriptome.bam -o S --filter-group no-filters --model-coverage` on transcriptome alignments or `--reads --annotated tx.fa --seq-tech ont-cdna`), FLAIR 3 (`flair align/correct/collapse/quantify`) |
| ONT direct RNA (RNA004) | `dorado basecaller rna004_130bps_sup@v5.2.0,m6A,pseU --estimate-poly-a pod5/` | `minimap2 -ax splice -uf -k14 --MD -y` (reads are sense); `--secondary=no` | same tools; `--data_type nanopore`; mod calls via `modkit pileup --motif DRACH 2` |
| scRNA/snRNA (10x + ONT) | nf-core/scnanoseq | — | — |

`paftools.js gff2bed anno.gtf > anno.bed12` builds the junction BED. Use GENCODE comprehensive GTF matched to the reference (`.claude/skills/genomics/reference.md`). Novel-isoform claims need ≥ 3 FLNC/full-length reads and SQANTI3/pigeon structural categories; report FSM/ISM/NIC/NNC counts. Differential isoform usage: `bambu` → `DEXSeq`/`satuRn`; long-read fusion detection: `pbfusion` (HiFi), `JAFFAL`/`ctat-LR-fusion` (ONT).

### Targeted and adaptive sampling

- **Adaptive sampling (ONT ReadUntil):** target BED (padded ≥ 10 kb per side) loaded in MinKNOW; expect 3–10× enrichment, ~400–600 bp rejected reads flagged `end_reason=data_service_unblock_mux_change` in `sequencing_summary.txt`. Compute on-target depth with `mosdepth --by targets.bed`; filter rejects (`< 1 kb`) before SV/STR calling; off-target reads still give a low-pass genome for CNV/aneuploidy. Enriched regions of hundreds of Mb reduce enrichment; keep targets ≤ 1–3 % of the genome.
- **Cas9/PCR panels (ONT), PureTarget / HiFi target capture (PacBio):** very high local depth (1 000–10 000×) — set `--max-depth`/`--preset targeted` (TRGT), raise `--min_support` for cuteSV/Sniffles, downsample with `samtools view -s` for Clair3 (`--chunk_num=-1` for amplicons).
- **Mitochondrial genome:** align to full reference (not chrM alone) to divert NUMTs; call haploid at high depth; rotate the reference or use `mutserve`/`mitorsaw` for the D-loop junction.
