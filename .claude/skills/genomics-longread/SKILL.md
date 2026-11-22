---
name: genomics-longread
description: Long-read DNA and RNA analysis for Oxford Nanopore (R10.4.1, dorado) and PacBio HiFi (Revio/Sequel II) — basecalling, alignment, long-read QC, small variants (Clair3/DeepVariant), phasing, SVs (Sniffles2/sawfish/pbsv/cuteSV), methylation (modkit/pb-CpG-tools), STRs (Straglr/TRGT), assembly and Iso-Seq/ONT cDNA. Use when inputs are pod5/fast5, dorado or HiFi uBAM (MM/ML tags), `*.hifi_reads.bam`, sequencing_summary.txt, or the user names minimap2 map-ont/map-hifi, Sniffles, modkit, Clair3, TRGT, hifiasm, Flye, Iso-Seq. Not for short-read WGS (`genomics-wgs`), short-read SV/CNV (`genomics-cnv`), somatic SNV/indel (`genomics-somatic`), or truth-set benchmarking mechanics (`genomics-benchmark`).
argument-hint: "[ont|hifi|rna] [path]"
paths:
  - "**/*.pod5"
  - "**/*.fast5"
  - "**/*.hifi_reads.bam"
  - "**/*.hifi_reads.*.bam"
  - "**/sequencing_summary*.txt"
  - "**/*.snf"
  - "**/*.svsig.gz"
  - "**/*.flnc.bam"
  - "**/*bedmethyl*"
---

# Long-read genomics (ONT + PacBio HiFi)

Expert in Oxford Nanopore and PacBio HiFi DNA/RNA analysis. Outputs are pinned commands, platform × question decision tables, QC thresholds with rationale, and validation plans. Long-read callers live only here; short-read equivalents live in the sibling skills linked at the end.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" dorado pod5 minimap2 pbmm2 samtools bcftools cramino NanoPlot pycoQC run_clair3.sh whatshap longphase hiphase sniffles sawfish pbsv cuteSV modkit aligned_bam_to_cpg_scores straglr.py trgt hifiasm flye isoseq lima oarfish isoquant || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version. Verified versions (2026-11): dorado 2.1.2 (DNA models `hac@v6.0.0`, `sup@v5.2.0`), pod5 0.3.49, minimap2 2.31, pbmm2 26.2.0 (SMRT Link 26.2), samtools ≥1.20, Clair3 2.0.5 (PyTorch; v1 TF models incompatible), DeepVariant 1.10.0, WhatsHap 2.8, longphase 2.1.2, HiPhase 1.7.1, Sniffles 2.8.1, sawfish 2.3.1, pbsv 2.11.0, cuteSV 2.1.4, modkit 0.6.4, pb-CpG-tools 3.0.0, Straglr 1.5.6, TRGT 5.1.0, hifiasm 0.25.0, Flye 2.9.7, isoseq 4.3.0 (bioconda tag v26.2.0), lima 26.2.1, bambu 3.12.1, oarfish 0.13.0, IsoQuant 4.1.0, FLAIR 3.0.1, NanoPlot 1.48.0, cramino 2.0.0, pycoQC 2.5.2 (unmaintained since 2020; summary-file QC only), Severus 1.7, SAVANA 1.3.8, QUAST 5.3.0, merqury 1.4.1, BUSCO 6.1.0.

## Inspect first

| Input | Command | Infer |
|---|---|---|
| Any BAM/uBAM | `samtools view -H x.bam \| grep -E '^@(RG\|PG\|SQ)' \| head -30` | Platform (`PL:ONT`/`PL:PACBIO`), basecaller + model in `@PG`/`@RG DS:`, already aligned? (`@SQ` present), reference contig naming (`chr1` vs `1`), sample `SM` |
| Dorado uBAM | `samtools view x.bam \| head -2 \| cut -f12- \| tr '\t' '\n' \| cut -c1-40` | `MM`/`ML` present (mods called), `mv` (moves; enables Clair3 `--enable_dwell_time`), `dx:i:1` (duplex), `qs:f` (mean qscore), `RG` model string → R9 vs R10, 4 kHz vs 5 kHz, hac vs sup |
| HiFi BAM | same, look for `rq:f`, `np:i`, `fi/fp/ri/rp` kinetics, `MM/ML` (5mC called on-instrument) | `rq ≥0.99` = HiFi; `fail_reads.bam` must be excluded; kinetics tags inflate size (`pbmm2` strips them) |
| pod5/fast5 directory | `ls *.pod5 \| wc -l; du -sh .` ; `pod5 inspect summary <one.pod5>` | Raw signal; fast5 → convert; sample rate 4 kHz (old) vs 5 kHz (current R10.4.1) decides model family |
| `sequencing_summary*.txt` | `head -1 \| tr '\t' '\n'` | MinKNOW run; `pycoQC`-ready; column `end_reason` reveals adaptive-sampling unblocks |
| FASTQ (headers) | `peek.sh x.fq.gz` | ONT header has `runid= model_version_id=`; HiFi names `m84.../ccs`; `MM:Z` in comment means tags survived export |
| Read length | `cramino x.bam` or `samtools stats x.bam \| grep ^RL \| tail -5` | N50 class: HiFi 15–20 kb, ONT ligation 10–40 kb, UL >100 kb, cDNA 0.5–3 kb |
| Reference | `ref_sanity.sh ref.fa x.bam` | GRCh38 no-alt vs T2T-CHM13 (see `.claude/skills/genomics/reference.md`); mismatch aborts |

Ask only if still ambiguous: platform/chemistry when no `@PG`/`@RG` (raw FASTQ with stripped headers); whether duplex was intended; sex karyotype for TRGT/haploid calling; tumour/normal pairing for somatic SV; target BED for adaptive sampling or panels; whether modification calls are needed (decides re-basecalling).

## Decision table

| Question | ONT R10.4.1 (dorado sup/hac v5.2+/v6.0) | PacBio HiFi (Revio/Sequel II) |
|---|---|---|
| Basecall / reads | `dorado basecaller sup,5mCG_5hmCG pod5/` → uBAM with MM/ML; duplex only for Q30 needs | Done on instrument: `*.hifi_reads.bam` (`rq≥0.99`); `lima` for barcodes |
| Align | `minimap2 -ax lr:hq -y --MD` (≥Q20 reads; `map-ont` for R9/noisy) | `pbmm2 align --preset HIFI --sort` (keeps MM/ML, strips kinetics) or `minimap2 -ax map-hifi -y` |
| Small variants | Clair3 `--platform=ont --model_path=r1041_e82_400bps_sup_v520` (or `hac_v600`); DeepVariant `--model_type ONT_R104` | DeepVariant `--model_type PACBIO` (natively phased since 1.10); Clair3 `hifi_revio` |
| Phase + haplotag | `whatshap phase` + `haplotag`, or `longphase phase --ont` (faster) | `hiphase` (joint SNV+SV+TRGT) |
| Germline SV | `sniffles --tandem-repeats <trf.bed>`; cuteSV ONT params as second caller | `sawfish discover/joint-call` (PacBio current) or `pbsv discover --tandem-repeats` + `pbsv call --ccs`; Sniffles2 also fine |
| Population SV | Sniffles2 `--snf` per sample → joint `.snf` call | `sawfish joint-call --sample-csv` or `pbsv call` on all `.svsig.gz` |
| Somatic SV (emerging) | Severus (T/N, phased BAMs) or SAVANA `--ont` | Severus; SAVANA `--pb` |
| Methylation | `modkit pileup` → bedMethyl; `modkit dmr` | `aligned_bam_to_cpg_scores` (pb-CpG-tools) or `modkit pileup` |
| STR / repeat expansion | Straglr (genome scan or `--loci`) | TRGT `genotype -b <catalog.bed>` (STRchive/Adotto) |
| CNV | Spectre (ONT) — see `.claude/skills/genomics-cnv/SKILL.md` | HiFiCNV / sawfish CNV track — same pointer |
| De novo assembly | Flye `--nano-hq`; hifiasm `--ont` (R10 sup) | hifiasm (+`--ul` ONT UL, `--h1/--h2` Hi-C) |
| Long-read RNA | `minimap2 -ax splice -uf -k14` (dRNA) / `splice:hq -uf` (cDNA); bambu, IsoQuant, oarfish | `lima --isoseq` → `isoseq refine --require-polya` → `cluster2`; `pbmm2 --preset ISOSEQ`; IsoQuant/bambu |
| Reference | GRCh38 no-alt default; T2T-CHM13 for SV/STR/centromere work | same; liftover rules in `.claude/skills/genomics/reference.md` |

Caller agreement rule: for clinical-grade SV sets run ≥2 callers and report the intersection separately from the union; record which caller supports each call.

## Workflow

1. **Basecall (ONT only; skip if uBAM exists).** Current input is pod5; convert fast5 first.

```bash
set -euo pipefail
POD5_DIR=pod5/
SAMPLE=S1
MODELS=/opt/dorado-models
if ls fast5/*.fast5 >/dev/null 2>&1; then pod5 convert fast5 fast5/*.fast5 --output "${POD5_DIR}" --one-to-one fast5/; fi
dorado basecaller sup,5mCG_5hmCG "${POD5_DIR}" --models-directory "${MODELS}" --emit-moves --min-qscore 10 > "${SAMPLE}.unaligned.bam"
dorado summary "${SAMPLE}.unaligned.bam" > "${SAMPLE}.summary.tsv"
```

Validate: `samtools view ${SAMPLE}.unaligned.bam | head -1 | grep -o -E 'MM:Z:[^ \t]{0,20}|mv:B' ; samtools view -H ${SAMPLE}.unaligned.bam | grep -o 'dna_r10[^ ]*'` — MM present and model string recorded. Use `hac,5mCG_5hmCG` when GPU-hours are the constraint (hac@v6.0.0 ≈ Q18 median vs sup@v5.2.0 ≈ Q20.6 on reads; sup still wins for assembly). Duplex, demux, R9 models: `reference.md` §1.

2. **Align.** Keep modification tags: dorado/HiFi BAM → FASTQ must carry `MM,ML` in the comment and minimap2 must copy comments with `-y`.

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
THREADS=32
PRESET=lr:hq
samtools fastq -@ 8 -T MM,ML,mv,qs,dx "${SAMPLE}.unaligned.bam" | minimap2 -t "${THREADS}" -ax "${PRESET}" -y --MD -Y --secondary=no -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT" "${REF}" - | samtools sort -@ 8 -m 2G -o "${SAMPLE}.bam" - && samtools index "${SAMPLE}.bam" && samtools quickcheck "${SAMPLE}.bam"
```

HiFi: `pbmm2 align --preset HIFI --sort -j ${THREADS} --sample ${SAMPLE} ${REF} ${SAMPLE}.hifi_reads.bam ${SAMPLE}.bam` (BAM in, sorted BAM out, tags preserved, kinetics stripped). Spliced RNA presets: `reference.md` §6.
Validate: `samtools view ${SAMPLE}.bam | head -1000 | grep -c 'MM:Z'` (>0) and `samtools flagstat ${SAMPLE}.bam | grep -E 'mapped|supplementary'` (mapped ≥95 %, supplementary 1–10 %).

3. **QC.** 

```bash
set -euo pipefail
cramino -t 8 --hist --karyotype "${SAMPLE}.bam" > "${SAMPLE}.cramino.txt"
NanoPlot -t 8 --bam "${SAMPLE}.bam" --N50 -o nanoplot_"${SAMPLE}"
```

ONT run-level: `pycoQC -f sequencing_summary.txt -a ${SAMPLE}.bam -o ${SAMPLE}.pycoQC.html`. Validate against the table in *QC thresholds*; `grep -E 'N50|Yield|identity|coverage' ${SAMPLE}.cramino.txt`.

4. **Small variants.** Model must match chemistry *and* basecaller version; `--include_all_ctgs` for non-human or unplaced contigs; set haploid handling for XY samples.

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
THREADS=32
MODEL=/opt/models/r1041_e82_400bps_sup_v520
run_clair3.sh --bam_fn="${SAMPLE}.bam" --ref_fn="${REF}" --threads="${THREADS}" --platform=ont --model_path="${MODEL}" --output=clair3_"${SAMPLE}" --enable_phasing --longphase_for_phasing --gvcf
```

DeepVariant alternative (ONT or HiFi): `run_deepvariant --model_type ONT_R104|PACBIO --ref ${REF} --reads ${SAMPLE}.bam --output_vcf ${SAMPLE}.dv.vcf.gz --output_gvcf ${SAMPLE}.dv.g.vcf.gz --num_shards ${THREADS} --haploid_contigs chrX,chrY --par_regions_bed GRCh38_PAR.bed` (XY sample; omit for XX). Clair3 move-table models (`*_with_mv` + `--enable_dwell_time`) need `mv` tags from `--emit-moves`.
Validate: `bcftools stats clair3_${SAMPLE}/merge_output.vcf.gz | grep -E '^SN|^TSTV' | head -12` — ~4.0–4.5 M PASS variants per 30× human genome, Ti/Tv 2.0–2.1, het/hom-alt ≈ 1.5–2.0.

5. **Phase and haplotag** (needed for methylation by haplotype, Severus, allele-specific work).

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
VCF=clair3_"${SAMPLE}"/merge_output.vcf.gz
whatshap phase --ignore-read-groups --reference "${REF}" -o "${SAMPLE}.phased.vcf.gz" "${VCF}" "${SAMPLE}.bam"
tabix -f "${SAMPLE}.phased.vcf.gz"
whatshap haplotag --ignore-read-groups --tag-supplementary --reference "${REF}" -o "${SAMPLE}.haplotagged.bam" "${SAMPLE}.phased.vcf.gz" "${SAMPLE}.bam"
samtools index "${SAMPLE}.haplotagged.bam"
whatshap stats --gtf "${SAMPLE}.phaseblocks.gtf" "${SAMPLE}.phased.vcf.gz" > "${SAMPLE}.phase_stats.txt"
```

HiFi: `hiphase --threads 16 --reference ${REF} --bam ${SAMPLE}.bam --output-bam ${SAMPLE}.haplotagged.bam --vcf ${SAMPLE}.dv.vcf.gz --output-vcf ${SAMPLE}.dv.phased.vcf.gz --vcf ${SAMPLE}.sv.vcf.gz --output-vcf ${SAMPLE}.sv.phased.vcf.gz --summary-file ${SAMPLE}.hiphase.tsv`. Faster ONT option: `longphase phase -s ${VCF} -b ${SAMPLE}.bam -r ${REF} -t 16 -o ${SAMPLE}.longphase --ont` then `longphase haplotag`.
Validate: `grep -E 'Phase block N50|phased' ${SAMPLE}.phase_stats.txt` — N50 ≥ 0.5 Mb at 30× ONT 20 kb N50; `samtools view ${SAMPLE}.haplotagged.bam | grep -c 'HP:i'` > 0.

6. **Structural variants.**

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
TRF=human_GRCh38_no_alt_analysis_set.trf.bed
sniffles --threads 16 --reference "${REF}" --tandem-repeats "${TRF}" --input "${SAMPLE}.haplotagged.bam" --vcf "${SAMPLE}.sniffles.vcf.gz" --snf "${SAMPLE}.snf" --output-rnames --minsvlen 50 --phase
bcftools index -f "${SAMPLE}.sniffles.vcf.gz"
```

HiFi primary: `sawfish discover --threads 16 --ref ${REF} --bam ${SAMPLE}.bam --output-dir sawfish_${SAMPLE}` then `sawfish joint-call --threads 16 --sample sawfish_${SAMPLE} --output-dir sawfish_${SAMPLE}_jc` (emits SV + CNV; needs expected-CN/excluded-region tracks for CNV). pbsv, cuteSV platform parameters, population merge, somatic (Severus/SAVANA): `reference.md` §3.
Validate: `bcftools view -H -f PASS ${SAMPLE}.sniffles.vcf.gz | wc -l` — 20–30 k PASS calls per human genome; `bcftools query -f '%SVTYPE\n' ${SAMPLE}.sniffles.vcf.gz | sort | uniq -c` — INS ≈ DEL, INS slightly more on T2T; BND <5 %.

7. **Methylation.**

```bash
set -euo pipefail
REF=GRCh38_no_alt.fa
SAMPLE=S1
modkit pileup --threads 16 --reference "${REF}" --cpg --combine-strands --phased --prefix "${SAMPLE}" "${SAMPLE}.haplotagged.bam" modkit_"${SAMPLE}"/
modkit summary --threads 8 "${SAMPLE}.haplotagged.bam" > "${SAMPLE}.modkit_summary.txt"
```

Unphased: drop `--phased` and give a single `.bedmethyl` output path instead of a directory. HiFi: `aligned_bam_to_cpg_scores --bam ${SAMPLE}.bam --output-prefix ${SAMPLE}.cpg --threads 16` (needs `rq`/MM tags). DMR, 5hmC, 6mA, direct-RNA mods: `reference.md` §4.
Validate: `awk '$10>=10' modkit_${SAMPLE}/${SAMPLE}_combined.bedmethyl | wc -l` (CpGs with ≥10 valid cov; expect ~25 M of 28 M at 30×); global CpG methylation 70–80 % for blood/tissue (`modkit summary` fraction of `m`).

8. **Tandem repeats (optional).** ONT: `straglr.py ${SAMPLE}.bam ${REF} straglr_${SAMPLE} --loci STRchive.hg38.bed --genotype_in_size --min_support 2 --nprocs 16`. HiFi: `trgt genotype --genome ${REF} --reads ${SAMPLE}.bam --repeats catalog.bed --karyotype XY --threads 16 --output-prefix ${SAMPLE}.trgt` then `bcftools sort`/index and `trgt plot` for any expanded locus. Catalog choice, clinical loci, merge: `reference.md` §5.
Validate: `bcftools view -H ${SAMPLE}.trgt.sorted.vcf.gz | wc -l` equals catalog rows; inspect `AL`/`MC` fields for loci above pathogenic thresholds.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Read-level accuracy (median) | ONT sup@v5.2 ≥ Q20 (99 %); hac@v6.0 ≥ Q18; HiFi `rq` ≥ 0.99 (Q20) by definition, median Q30+ | Below → wrong model/chemistry (R9 data with R10 model) or degraded flowcell; HiFi `fail_reads` must not be in the input |
| Read N50 | ONT ligation ≥ 10 kb (target 15–30 kb); UL ≥ 50 kb; HiFi 15–20 kb; cDNA 1–2 kb | Short N50 kills SV/phasing sensitivity; rapid kits and FFPE legitimately give 3–8 kb |
| Depth (mean, aligned) | Germline SNV/SV 30× (min 20×); methylation 20×; somatic SV ≥ 60× T / 30× N; STR targeted ≥ 20 spanning reads | Clair3/DeepVariant F1 plateaus ~30×; phasing N50 grows with depth to ~40× then with read length |
| Mapped fraction | ≥ 95 % (ONT), ≥ 98 % (HiFi) | Lower → contamination, adapter chimeras, wrong reference, or unfiltered adaptive-sampling rejects |
| Supplementary alignments | 1–10 % of records | >15 % → chimeric library, reference mismatch, or `-Y` missing; needed intact for SV callers |
| Gap-compressed identity | ONT ≥ 97 %; HiFi ≥ 99.5 % (cramino) | Reports read errors + SV divergence; drop expectation on T2T-divergent regions |
| Duplex fraction (if run) | 10–40 % of bases | Only duplex (`dx:i:1`) reach Q30; low fraction means the library/kit was not duplex-optimised |
| Adaptive-sampling on-target | ≥ 5× enrichment; rejected reads ≈ 400–600 bp | Count `end_reason=data_service_unblock_mux_change`; filter rejects (<1 kb) before SV calling |
| SNV Ti/Tv (WGS) | 2.0–2.1 | <1.9 → false positives from homopolymer indels miscounted or model mismatch |
| SV count (PASS, human) | 20–30 k; INS:DEL ≈ 1:1 | <15 k → low depth/short reads; >40 k → noisy reads, no TR annotation, or T2T vs GRCh38 mixed |
| Phase block N50 | ≥ 0.5 Mb (30×, 20 kb N50 ONT); ≥ 1 Mb HiFi + HiPhase with SVs | Short blocks → low het density (consanguinity) or low depth |
| CpG sites ≥ 10× | ≥ 85 % of ~28 M (GRCh38) | Lower → low depth or MM/ML lost during alignment |
| Isoform/RNA: FLNC yield | ≥ 70 % of HiFi reads FLNC; mapped spliced ≥ 90 % | Low FLNC → primer/poly(A) issue; use `--require-polya` only for poly(A)-selected libraries |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `modkit pileup` reports 0 mod calls / `aligned_bam_to_cpg_scores` empty | MM/ML lost: FASTQ export without tags, or realignment without `-y` | `samtools fastq -T MM,ML` → `minimap2 -y`; or `dorado aligner`/`pbmm2` on the uBAM directly; never re-align from plain FASTQ |
| Sniffles/cuteSV find almost no SVs, or thousands of tiny BNDs | Supplementary alignments dropped (`samtools view -F 2048`), hard clipping (`-Y` absent), or `--secondary=yes` noise | Keep flags 2048, use `-Y`, drop secondaries; pbmm2 already enforces this |
| Variant QUAL/F1 far below expectation on ONT | R9 data called with R10 model (or 4 kHz vs 5 kHz, hac vs sup mismatch) in Clair3/DeepVariant | Read `@RG DS:basecall_model=` and pick the matching Clair3 model (table in `reference.md` §2); re-basecall legacy data with current dorado if pod5/fast5 exist |
| Clair3 v2 refuses legacy model directory | TensorFlow (v1/Rerio) models incompatible with v2 PyTorch | Use `clair3_models_pytorch` / converted Rerio set; bundled models live under `/opt/models` or `${CONDA_PREFIX}/bin/models` |
| Hundreds of het SNVs/indels on chrM, or SVs anchored in chrM | NUMT reads mis-placed; chrM circular junction; haploid treated as diploid | Call chrM separately (haploid, high depth), mask NUMT loci, ignore `MT` BNDs; use `--haploid_contigs` only for chrX/chrY not chrM |
| Male sample shows het calls across chrX | Diploid genotyping on hemizygous contigs | DeepVariant `--haploid_contigs chrX,chrY --par_regions_bed`; Clair3 `--gender` (≥2.0.3) or split by PAR |
| `pbmm2`/`pbsv` complain about missing `SM` | HiFi FASTQ input without read group | `--rg '@RG\tID:x\tSM:y'` or `--sample`; joint calling keys on `SM` |
| BAM write error "CIGAR > 65535 ops" | Ultra-long reads | `minimap2 -L` (moves CIGAR to `CG` tag) or output CRAM; pbmm2 does this implicitly |
| Duplex reads double-count SV/methylation evidence | Simplex parents emitted together with duplex children | Filter `dx:i:-1` (simplex parent of a duplex pair) before calling: `samtools view -e '![dx] \|\| [dx]!=-1'` |
| Adaptive-sampling run has huge count of ~500 bp reads | Rejected (unblocked) molecules | Filter `< 1 kb` and `end_reason` rejects; compute on-target depth separately |
| Half the CpGs missing in bedMethyl, methylation ~50 % everywhere | Strands not combined or hemi-methylation on one haplotype reported as unphased | `--cpg --combine-strands`; partition by `HP` for allele-specific methylation |
| T2T-CHM13 run yields more INS than DEL and extra calls in centromeres | Expected: GRCh38 gaps closed; satellite arrays align but with low MAPQ | Report by MAPQ bin; use T2T TR annotation; compare with GRCh38 only inside confident regions |
| `isoseq cluster` out of memory on Kinnex data | Legacy clustering | `isoseq cluster2`; run `skera split` first for MAS-Seq arrays |
| Iso-Seq/dRNA alignments mostly antisense or unspliced | Missing `-uf`, wrong preset, or no `--junc-bed` for noisy reads | `-ax splice:hq -uf` (cDNA/Iso-Seq), `-ax splice -uf -k14` (dRNA), add `--junc-bed` from GTF |

## nf-core

Generic Nextflow/profile mechanics: `.claude/skills/genomics-workflow/SKILL.md`. All pipelines below exist in `https://nf-co.re/pipelines.json` (checked 2026-11-22).

| Pipeline | Latest | Scope | Samplesheet header |
|---|---|---|---|
| nf-core/nanoseq | 3.1.0 | ONT DNA + cDNA/dRNA: demux, minimap2/graphmap2, QC, Clair3/medaka SNV, Sniffles/cuteSV SV, modification, isoform quant (bambu/StringTie2, FLAIR) | `group,replicate,barcode,input_file,fasta,gtf` |
| nf-core/longraredisease | 1.0.0 | ONT/PacBio rare disease: multi-caller SV consensus (Sniffles, cuteSV, SVIM + Jasmine), optional Clair3/DeepVariant, Spectre/HiFiCNV, Straglr, modkit, HPO prioritisation; input uBAM/BAM | `sample,file_path,hpo_terms,sex,phenotype,family_id,maternal_id,paternal_id` |
| nf-core/scnanoseq | 1.3.0 | ONT single-cell/nucleus cDNA (10x) → barcode correction, isoform/gene matrices | `sample,fastq,cell_count` |
| nf-core/pacvar | 1.1.0 | PacBio HiFi WGS and PureTarget: pbmm2, DeepVariant, pbsv/sawfish, TRGT, HiPhase | `sample,bam,pbi` |
| nf-core/methylong | 2.0.0 | ONT + PacBio 5mC: align → methylation → SNV → phasing → DMR | `group,sample,path,ref,method` (`method` = `ont` or `pacbio`) |
| nf-core/isoseq | 3.0.1 | PacBio Iso-Seq: lima → refine → cluster → align → TAMA/isoform collapse | `sample,seq_data,pbi,start_from` |

Vendor workflows (not nf-core, but current references): `epi2me-labs/wf-human-variation` 2.8.1 (ONT SNV/SV/CNV/STR/mod), `PacificBiosciences/HiFi-human-WGS-WDL` 4.1.0.

## Validation

- Small variants: benchmark on GIAB HG002 with `hap.py`/`vcfdist` using the **T2T-Q100 (v1.1) benchmark** for whole-genome and HG002 v4.2.1 for GRCh38 legacy comparisons; expected ONT sup 30× SNP F1 ≥ 0.995, indel F1 ≥ 0.90 (homopolymer indels dominate FN/FP); HiFi SNP ≥ 0.998, indel ≥ 0.98. Mechanics: `.claude/skills/genomics-benchmark/SKILL.md`.
- SVs: Truvari (`bench`, `refine`) against **GIAB HG002 SV v1.1** (T2T-based, replaces v0.6) and **CMRG v1.00** for medically relevant genes; expect ≥ 0.9 F1 at 30× for INS/DEL ≥ 50 bp in confident regions; report `--pctseq 0` and default separately.
- Phasing: switch error rate < 0.5 % vs GIAB trio-phased truth (`whatshap compare`).
- Methylation: correlation ≥ 0.9 with WGBS/EM-seq on the same sample at CpGs ≥ 10×; check imprinted loci (H19/IGF2, SNRPN) show ~50 % and allele-specific split in `hp1`/`hp2` bedMethyl after `--phased`.
- STR: known expansions in GIAB HG002/HG001 (e.g. STRchive truth loci) recovered within ±1 motif unit; visualise with `trgt plot` / Straglr spanning reads.
- Assembly: QUAST vs T2T-CHM13, merqury QV ≥ 50 (HiFi) / ≥ 40 (ONT sup+polish), BUSCO ≥ 98 % complete, hifiasm phased contig N50 > 50 Mb.
- Always: `ref_sanity.sh` across BAM/VCF/BED; re-run caller on a 1-chromosome subset after any model/version change and diff counts before full run.

## Additional resources

- `.claude/skills/genomics-longread/reference.md` — dorado details (duplex, demux, R9/4 kHz legacy models, pod5), Clair3/DeepVariant model matrix, cuteSV/pbsv/population/somatic SV commands, modkit DMR + pb-CpG-tools options, STR catalogs, assembly (hifiasm/Flye/QC), long-read RNA (Iso-Seq, ONT cDNA/dRNA, bambu/oarfish/IsoQuant/FLAIR), adaptive sampling. Load when a step above says "`reference.md` §n".
- `.claude/skills/genomics/reference.md` — GRCh38 vs T2T-CHM13 choice, downloads, TR/PAR/exclusion BEDs, liftover. Load before picking a reference.
- `.claude/skills/genomics-benchmark/SKILL.md` — GIAB truth sets, hap.py/Truvari mechanics. Load for any accuracy claim.
- `.claude/skills/genomics-cnv/SKILL.md` — CNV calling incl. long-read Spectre/HiFiCNV context. Load for copy-number questions.
- `.claude/skills/genomics-somatic/SKILL.md` — somatic SNV/indel (ClairS/DeepSomatic pointers). Load for tumour small variants; somatic SV stays here.
- `.claude/skills/genomics-annotation/SKILL.md` — VEP/SV annotation and ACMG/ClinGen rules once VCFs exist.
- `.claude/skills/genomics-rnaseq/SKILL.md` — short-read RNA; load when comparing long-read isoform results with Illumina quantification.
- `.claude/skills/genomics-qc/SKILL.md` — generic read QC/trimming (short-read); long-read QC stays here.
- `.claude/skills/genomics-workflow/SKILL.md` — nf-core/Nextflow/SLURM mechanics.
