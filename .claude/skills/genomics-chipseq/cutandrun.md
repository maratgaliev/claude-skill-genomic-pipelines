# CUT&RUN and CUT&Tag — deviations from the ChIP-seq core

Load after `SKILL.md`. Shared steps (dedup, blacklist, deepTools, DiffBind, annotation, motifs) are unchanged. What differs: antibody-tethered pA-MNase (CUT&RUN) or pA-Tn5 (CUT&Tag) releases short, low-background fragments; there is no input, the control is IgG; carry-over E. coli DNA from the enzyme prep is the natural spike-in; and duplicates in target samples are often genuine (identical cut positions at a tethered site), so dedup is applied to IgG by default and to targets only when PCR bottlenecking is evident.

## Inspect first (CUT&RUN-specific)

| Read | Infer |
|---|---|
| Fragment-size histogram (`samtools view -f 66 \| cut -f9`) | TF/CUT&RUN: dominant < 120 bp; histone: 150–300 bp nucleosomal ladder; CUT&Tag: 10 bp sawtooth on top |
| Sample names / design | which samples are IgG, whether one IgG per condition or a shared one |
| Reads mapping to E. coli (`bowtie2` against K12 MG1655) | spike fraction 0.01–10 %; > 10 % means little target chromatin was recovered |
| Read length | 25–50 bp PE is normal for CUT&Tag (fragments are short); adapters read-through is universal → trim per `genomics-qc` |
| Kit | commercial kits with purified pA-Tn5 may carry almost no E. coli; then spike-in normalisation is impossible unless exogenous spike DNA was added |

Ask only if still ambiguous: CUT&RUN vs CUT&Tag (Tn5 duplicates and the 9 bp insertion offset apply only to CUT&Tag), whether exogenous spike-in (yeast/Drosophila) was added, and whether cell numbers were equal across samples.

## Decision table

| Situation | Choice |
|---|---|
| Alignment to target genome | Bowtie2 `--local --very-sensitive-local --no-mixed --no-discordant --dovetail -I 10 -X 700` (CUT&RUN; nf-core/cutandrun default); CUT&Tag with ≥ 50 bp reads may use `--end-to-end --very-sensitive` instead |
| Alignment to spike genome | same settings plus `--no-overlap --no-dovetail` to reject cross-mapping fragments |
| Normalisation, equal cell input, measurable E. coli | spike-in scale factor (`C / spike_reads`) |
| Normalisation, kit with no E. coli, or unequal cells | CPM via `bamCoverage --normalizeUsing CPM --binSize 1`; say so |
| Peak calling, sharp TF or H3K4me3 | SEACR `stringent` with IgG; MACS3 `-f BAMPE -q 0.01` as cross-check |
| Peak calling, broad mark (H3K27me3, H3K9me3) | SEACR `relaxed` with IgG; MACS3 `--broad` cross-check |
| No IgG available | SEACR numeric threshold `0.01` (top 1 % of regions) `non stringent`; MACS3 without `-c` — flag as unreplicated-control |
| Dedup target reads? | only if PCR duplication is high (`NRF < 0.8`) or the library was over-amplified; always dedup IgG |
| Differential binding | DiffBind with `dba.normalize(spikein = <E. coli BAMs>)` when spike is valid, else `DBA_NORM_LIB` + background |

## Workflow

### 1. Align to target and spike-in

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
REF_IDX=/refs/GRCh38/bowtie2/GRCh38
SPIKE_IDX=/refs/ecoli_K12_MG1655/bowtie2/K12
THREADS=8
bowtie2 -p "${THREADS}" --local --very-sensitive-local --no-mixed --no-discordant --dovetail \
  -I 10 -X 700 -x "${REF_IDX}" -1 "${SAMPLE}_R1.trim.fq.gz" -2 "${SAMPLE}_R2.trim.fq.gz" \
  2> "${SAMPLE}.target.bowtie2.log" \
  | samtools view -b -q 20 -F 1804 -f 2 - \
  | samtools sort -@ 4 -o "${SAMPLE}.target.bam" - && samtools quickcheck "${SAMPLE}.target.bam"
bowtie2 -p "${THREADS}" --local --very-sensitive-local --no-mixed --no-discordant --no-overlap --no-dovetail \
  -I 10 -X 700 -x "${SPIKE_IDX}" -1 "${SAMPLE}_R1.trim.fq.gz" -2 "${SAMPLE}_R2.trim.fq.gz" \
  2> "${SAMPLE}.spike.bowtie2.log" \
  | samtools view -b -q 20 -f 2 - \
  | samtools sort -@ 4 -o "${SAMPLE}.spike.bam" - && samtools quickcheck "${SAMPLE}.spike.bam"
samtools index "${SAMPLE}.target.bam"
samtools index "${SAMPLE}.spike.bam"
```

`--dovetail` allows mates that extend past each other (short fragments with read-through); `--local` tolerates untrimmed adapter tails. `-X 700` caps fragments at tri-nucleosome. MAPQ 20 is the nf-core default; raise to 30 for TF targets. Validation: `grep overall "${SAMPLE}".*.bowtie2.log`.

### 2. Spike-in scale factor

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
CONST=10000
SPIKE_FRAGS=$(samtools view -c -f 66 "${SAMPLE}.spike.bam")
TARGET_FRAGS=$(samtools view -c -f 66 "${SAMPLE}.target.bam")
awk -v s="${SPIKE_FRAGS}" -v t="${TARGET_FRAGS}" -v c="${CONST}" \
  'BEGIN{printf "spike_frags\t%d\ntarget_frags\t%d\nspike_fraction\t%.5f\nscale_factor\t%.6f\n", s, t, s/(s+t), c/s}' \
  | tee "${SAMPLE}.spike.tsv"
```

Scale factor = `CONST / spike fragments` (Henikoff lab convention `CONST=10000`; any constant works as long as it is identical across the comparison). Valid only when cell number and enzyme amount are constant per sample; if spike fragments < 1000 the factor is noise — fall back to CPM. Validation: spike fraction similar (within ~3×) across replicates of the same antibody.

### 3. Dedup policy and fragment BED

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
IGG=IgG_rep1
gatk MarkDuplicates -I "${IGG}.target.bam" -O "${IGG}.dedup.bam" -M "${IGG}.dup_metrics.txt" --REMOVE_DUPLICATES false --CREATE_INDEX true
gatk MarkDuplicates -I "${SAMPLE}.target.bam" -O "${SAMPLE}.marked.bam" -M "${SAMPLE}.dup_metrics.txt" --REMOVE_DUPLICATES false --CREATE_INDEX true
samtools sort -n -@ 4 "${SAMPLE}.marked.bam" \
  | bedtools bamtobed -bedpe -i - \
  | awk 'BEGIN{OFS="\t"} $1==$4 && $6-$2<1000 {print $1,$2,$6}' \
  | sort -k1,1 -k2,2n > "${SAMPLE}.fragments.bed"
```

Targets are marked, not removed (decide after reading `PERCENT_DUPLICATION` and NRF); IgG duplicates are also only marked and are excluded by flag (`samtools view -F 1024 "${IGG}.dedup.bam"`) when IgG fragments are extracted or counted. Physical removal (`--REMOVE_DUPLICATES true`) on either only when the user asks for it. Blacklist filtering as `SKILL.md` §3 (`bedtools intersect -v` on the fragment BED with `-a`/`-b` works directly on BED). Validation: `awk '{print $3-$2}' "${SAMPLE}.fragments.bed" | sort -n | uniq -c | sort -k1,1nr | head` shows the expected fragment mode.

### 4. Fragment-size expectations

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
awk '{l=$3-$2; if(l<=700) h[l]++} END{for(i in h) print i"\t"h[i]}' "${SAMPLE}.fragments.bed" | sort -k1,1n > "${SAMPLE}.fraglen.tsv"
awk '$1<120{s+=$2} $1>=150&&$1<=300{n+=$2} {t+=$2} END{printf "sub120\t%.3f\nnucleosomal\t%.3f\n", s/t, n/t}' "${SAMPLE}.fraglen.tsv"
```

| Target | Expected | Deviation meaning |
|---|---|---|
| TF (CTCF, SOX2) CUT&RUN | majority < 120 bp, small ~180 bp shoulder | nucleosomal-dominated = over-digestion or antibody bound to chromatin non-specifically |
| Histone marks | ladder at ~180, 360, 540 bp; minor sub-120 | sub-120 dominant = MNase over-digestion |
| CUT&Tag (any) | nucleosomal ladder plus 10 bp periodicity (Tn5 helical preference) | absence of ladder with huge sub-100 = over-tagmentation |
| IgG | flat, low-count, broad distribution | an IgG with a crisp ladder is carry-over from the wash, not noise; still usable as background |

Validation: histogram saved per sample (`bamPEFragmentSize` on `target.bam` also works).

### 5. Scaled bedGraph and SEACR 1.3

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
IGG=IgG_rep1
GENOME_SIZES=/refs/GRCh38/GRCh38.fa.fai
SCALE=$(awk '$1=="scale_factor"{print $2}' "${SAMPLE}.spike.tsv")
IGG_SCALE=$(awk '$1=="scale_factor"{print $2}' "${IGG}.spike.tsv")
bedtools genomecov -bg -scale "${SCALE}" -i "${SAMPLE}.fragments.bed" -g "${GENOME_SIZES}" > "${SAMPLE}.scaled.bedgraph"
bedtools genomecov -bg -scale "${IGG_SCALE}" -i "${IGG}.fragments.bed" -g "${GENOME_SIZES}" > "${IGG}.scaled.bedgraph"
bash SEACR_1.3.sh "${SAMPLE}.scaled.bedgraph" "${IGG}.scaled.bedgraph" non stringent "seacr/${SAMPLE}.igg"
bash SEACR_1.3.sh "${SAMPLE}.scaled.bedgraph" 0.01 non stringent "seacr/${SAMPLE}.top01"
```

SEACR arguments: `<target.bedgraph> <IgG.bedgraph | numeric threshold> <norm|non> <stringent|relaxed> <prefix>`. Use `non` when both bedGraphs are already spike-scaled (SEACR `norm` would rescale IgG a second time); use `norm` with unscaled bedGraphs. `stringent` = threshold at the maximum of the IgG curve; `relaxed` = at the knee; broad marks → `relaxed`. Output `prefix.stringent.bed` has columns chrom, start, end, total signal, max signal, max-signal region. Validation: `wc -l seacr/${SAMPLE}.igg.stringent.bed` and fraction of fragments in peaks (`bedtools intersect -u -a fragments.bed -b peaks | wc -l`).

### 6. MACS3 cross-check

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
IGG=IgG_rep1
GSIZE=hs
macs3 callpeak -t "${SAMPLE}.target.bam" -c "${IGG}.dedup.bam" -f BAMPE -g "${GSIZE}" -n "${SAMPLE}" \
  --outdir macs3 -q 0.01 --keep-dup all --broad --broad-cutoff 0.1
```

Drop `--broad --broad-cutoff` for TF/H3K4me3. `--keep-dup all` because CUT&RUN duplicates at a tethered site are expected; MACS3 cannot use the spike-in factor, so IgG background is estimated from read counts (`--scale-to small` when IgG is much deeper). Validation: SEACR and MACS3 peak sets overlap > 70 % for sharp targets; large disagreement points to IgG scaling.

### 7. Tracks and differential binding

```bash
set -euo pipefail
SAMPLE=H3K27me3_rep1
GENOME_SIZES=/refs/GRCh38/GRCh38.fa.fai
sort -k1,1 -k2,2n "${SAMPLE}.scaled.bedgraph" > "${SAMPLE}.scaled.sorted.bedgraph"
bedGraphToBigWig "${SAMPLE}.scaled.sorted.bedgraph" "${GENOME_SIZES}" "${SAMPLE}.scaled.bw"
```

Without a valid spike-in: `bamCoverage -b "${SAMPLE}.target.bam" --normalizeUsing CPM --binSize 1 --extendReads -o "${SAMPLE}.cpm.bw"`. DiffBind: build the sample sheet with SEACR or MACS3 peaks, `bamControl` = IgG, then `dba.normalize(dba, spikein = c("rep1.spike.bam", ...))` (DiffBind ≥ 3.0 computes spike-based library sizes) or `normalize = DBA_NORM_LIB, background = TRUE` when no spike. Validation: `dba.plotPCA` separates conditions, IgG samples cluster away from targets.

## QC thresholds

No ENCODE standard exists for CUT&RUN/CUT&Tag; defaults derive from Meers 2019 (CUT&RUN), Kaya-Okur 2019 / Zheng 2020 protocol (CUT&Tag) and nf-core/cutandrun reporting. State this when reporting.

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| Sequencing depth | 3–8 M PE reads per histone sample; 5–10 M for TFs; IgG ≥ 2 M | low background means saturation at ≤ 10 M; deeper only for very broad marks |
| Alignment rate (target) | > 80 % | lower = adapter read-through untrimmed, or heavy E. coli fraction |
| Spike-in fraction | 0.1–10 % | > 10–20 %: little target chromatin; < 0.01 %: factor unusable, use CPM |
| Duplication (target, not removed) | < 20 % for histone at ≤ 8 M reads; higher acceptable for TFs with few sites | high dup + high NRF-by-position = tethered-site duplicates, not PCR |
| Duplication (IgG) | expect high; always removed | non-specific low-complexity |
| Fragment size | target-class pattern per §4 table | mismatch = digestion/tagmentation problem |
| FRiP (SEACR stringent) | > 0.2 histone; > 0.05 TF | far above ChIP because background is near zero; depends on caller |
| Peaks (SEACR stringent, IgG) | H3K27me3 20–60k, H3K4me3 10–25k, CTCF 20–40k | counts far outside suggest wrong `stringent/relaxed` or scaling |
| Replicate Pearson (1 kb bins, log CPM) | > 0.8 | lower for very sparse TF data |
| IgG peaks | < 5 % of target peak count | an IgG that "calls" thousands of peaks had carry-over or an unscaled bedGraph |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| SEACR returns zero peaks with IgG | both bedGraphs scaled and `norm` used, or IgG bedGraph far deeper and unscaled | `non` with scaled files; or scale IgG by its own spike factor |
| SEACR peaks cover the whole genome | `relaxed` with a numeric threshold `0.1`+ | threshold `0.01`, `stringent` |
| Very low alignment rate | `--end-to-end` on reads with adapter tails | `--local` or trim fully; see `genomics-qc` |
| Spike factor varies 50× between replicates | kit with trace E. coli, unequal cells, or spike BAM accepted dovetailed/overlapping pairs | check `--no-overlap --no-dovetail`; abandon spike scaling |
| Nucleosomal ladder in a TF CUT&RUN | over-digestion (MNase time/temperature) | call peaks on sub-120 bp fragments only (`awk '$3-$2<120'`); note in methods |
| Huge duplicate rate after `--REMOVE_DUPLICATES true` on targets collapsed peaks | genuine tethered duplicates removed | keep duplicates for targets; mark IgG and exclude with `-F 1024` |
| MACS3 finds 10× fewer peaks than SEACR | IgG used as `-c` is deeper than target; MACS3 scales target down | `--scale-to small` or SEACR primary |
| `bedtools genomecov` errors on chromosome order | fragments BED sorted differently from `.fai` | `sort -k1,1 -k2,2n` and a `.fai` of the same reference (`ref_sanity.sh`) |
| nf-core/cutandrun reports "spike-in alignment 0 %" | `--spikein_genome` left at E. coli while yeast/Drosophila spike was used | set `--spikein_genome R64-1-1` or `BDGP6`, or `--spikein_fasta` |

## nf-core/cutandrun (3.2.2)

```csv
group,replicate,fastq_1,fastq_2,control
h3k27me3,1,h3k27me3_1_R1.fastq.gz,h3k27me3_1_R2.fastq.gz,igg_ctrl
h3k27me3,2,h3k27me3_2_R1.fastq.gz,h3k27me3_2_R2.fastq.gz,igg_ctrl
igg_ctrl,1,igg_1_R1.fastq.gz,igg_1_R2.fastq.gz,
igg_ctrl,2,igg_2_R1.fastq.gz,igg_2_R2.fastq.gz,
```

Verified 2026-11-22 (`https://nf-co.re/cutandrun/usage`): five columns, `control` names a `group`; IgG replicates are matched by replicate number when counts agree, otherwise replicate 1 of the control group is used for all. Defaults: dedup IgG only (`--dedup_target_reads` to dedup targets), `--normalisation_mode Spikein` (alternatives `CPM|RPKM|BPM|None`), `--peakcaller seacr` (`seacr,macs2` runs both; first is primary), `--spikein_genome K12-MG1655`, `--minimum_alignment_q_score 20`, `--remove_mitochondrial_reads false`.

```bash
set -euo pipefail
nextflow run nf-core/cutandrun -r 3.2.2 -profile docker --input samplesheet.csv --outdir results \
  --fasta GRCh38.fa --gtf GRCh38.gtf --blacklist hg38-blacklist.v2.bed \
  --spikein_genome K12-MG1655 --normalisation_mode Spikein --peakcaller seacr,macs2 --use_control true
```

Generic Nextflow mechanics: `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

- Positive-control loci: H3K4me3 at GAPDH/ACTB promoters, H3K27me3 over HOX clusters, CTCF motif (JASPAR 2024 MA0139) centred in ≥ 60 % of top-1000 CTCF peaks.
- IgG track is flat at the same loci on the same y-scale.
- Spike-scaled replicate tracks of the same antibody overlay within ~2×; if they do not, report CPM tracks alongside.
- Peak overlap with a public ChIP-seq set for the same mark/cell type (`bedtools jaccard`); CUT&RUN peaks are typically narrower, so use `-f 0.3` reciprocal overlap rather than Jaccard alone.
- Benchmark mechanics: `.claude/skills/genomics-benchmark/SKILL.md`.
