---
name: genomics-benchmark
description: Benchmarks variant calls (germline SNV/indel, SV, CNV, somatic) against truth sets with hap.py, rtg vcfeval, Truvari and som.py, choosing the right GIAB/SEQC2 truth set, confident regions and stratifications, and interpreting precision/recall/F1. Use when the user says "benchmark", "truth set", "GIAB", "HG002", "hap.py", "vcfeval", "Truvari", "precision/recall", "F1", "validate my caller", or supplies a GIAB/SEQC2 VCF. Not for producing the calls (genomics-wgs, genomics-somatic, genomics-longread, genomics-cnv) or running pipelines (genomics-workflow).
argument-hint: "[small|sv|cnv|somatic] [query.vcf.gz]"
paths:
  - "**/*benchmark*.vcf.gz"
  - "**/HG00[1-7]_*v4.2.1*"
  - "**/*v5.0q*"
  - "**/*.extended.csv"
  - "**/summary.json"
---

# Genomics benchmark

Expert in GA4GH-style variant benchmarking: picks the truth set and confident region that match the assay and reference build, runs the comparison engine correctly, stratifies, and reports numbers that can be defended. Outputs are commands, decision tables and a filled results table; never a bare F1.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" hap.py som.py rtg truvari bcftools bedtools tabix docker singularity || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version. hap.py has no conda build for recent Python; run it from `jmcdani20/hap.py:v0.3.12` (the image used by DeepVariant's own case studies; `pkrusche/hap.py` is the legacy upstream image). Truvari: `pip install truvari==5.5.0` or `quay.io/biocontainers/truvari`.

## Inspect first

| Read | How | Infer |
|---|---|---|
| Query VCF header | `peek.sh query.vcf.gz` → `##reference`, `##contig` lines, `##source`/`##DeepVariant_version`/`##GATKCommandLine` | build (GRCh37/38/CHM13), `chr` prefix, caller + version, single vs multi-sample |
| Query body | `bcftools view -H query.vcf.gz \| head`; `bcftools query -l` | sample name, FILTER usage, gVCF (`<NON_REF>`) vs VCF, symbolic SV ALTs (`<DEL>`) vs sequence-resolved |
| Truth VCF + BED | same; `cut -f1 confident.bed \| sort -u` | truth build, contig naming, whether BED covers X/Y |
| Reference | `ref_sanity.sh ref.fa truth.vcf.gz query.vcf.gz confident.bed` | exits 1 on contig/length mismatch: stop and fix before anything else |
| Target BED (WES/panel) | `peek.sh targets.bed` | whether a `-T` restriction is needed; padding already applied? |
| Sample identity | header sample name, `##DRAGENCommandLine`, run sheet | which HG00x; sex (HG002/HG003/HG005 male; HG001/HG004/HG006/HG007 female) |

Ask only if still ambiguous: which GIAB sample was sequenced; whether the query is pre-filtered (PASS only) or raw; assay target BED for WES/panel; whether the user wants genotype-aware (default) or allele-only matching.

## Decision table

| Situation | Truth set | Confident region | Engine |
|---|---|---|---|
| Germline short-read WGS/WES, GRCh37/38, any HG001–HG007 | GIAB v4.2.1 small variants | `*_v4.2.1_benchmark_noinconsistent.bed` (or `_benchmark.bed` for HG001/5/6/7) | `hap.py --engine vcfeval` |
| Germline HG002, want hardest regions / CHM13 / long-read or assembly-based calls | GIAB HG002 v5.0q (assembly-based, from T2T-Q100 v1.1; still labelled *draft*) | `HG002_<build>_v5.0q_smvar.benchmark.bed` | `hap.py --engine vcfeval --gender male` + v3.6 stratifications |
| Medically relevant genes (HG002) | CMRG v1.00 small variants (+ CMRG SV) | `HG002_GRCh38_CMRG_smallvar_v1.00.bed` | hap.py; report separately from genome-wide |
| SV, short- or long-read, HG002 | GRCh37: `NIST_SV_v0.6` Tier1 (mature, GRCh37 only); GRCh38/CHM13: v5.0q `stvar` | Tier1 BED / `stvar.benchmark.bed` | `truvari bench` (+ `--refine` for v5.0q) |
| CNV / depth-only caller | v5.0q `stvar` DEL/DUP ≥ 1 kb subset, or an orthogonal array/dPCR set | `stvar.benchmark.bed` | `truvari bench --pctseq 0 --pctsize 0.7`; wittyer / `rtg cnveval` via nf-core/variantbenchmarking |
| Somatic SNV/indel, tumor-normal | SEQC2 HCC1395 vs HCC1395BL v1.2.1 (GRCh38) | `High-Confidence_Regions_v1.2.bed` | `som.py` (hap.py image) or `hap.py --somatic`-style allele matching; stratify by VAF |
| Somatic in non-reference tumour | no public truth; use dilution series / replicate concordance / orthogonal validation | assay targets | report concordance, not F1 |
| Non-human / no truth set | Mendelian consistency on trios, Ti/Tv, het/hom, known-site concordance | — | see Validation |
| Two callers, no truth | concordance only | — | `rtg vcfeval` with one as "baseline", label it as such |
| Multi-caller, multi-sample, reproducible report | nf-core/variantbenchmarking 1.6.1 | pipeline params | see nf-core section |

Rules: truth and query must be on the **same build and contig naming**; never liftover the query to match the truth if the truth exists for the query build. Benchmark **PASS calls only** unless the point is to tune filters (then report both). WES/panel: pass the assay BED with `-T`/`--targets_bed` **and** the confident BED with `-f` — both, not one.

## Workflow

1. Fetch truth files (checksum them; record URL + date in the results table).

```bash
set -euo pipefail
BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
SAMPLE_DIR=AshkenazimTrio/HG002_NA24385_son
BUILD=GRCh38
OUT=truth
mkdir -p "${OUT}"
for f in HG002_${BUILD}_1_22_v4.2.1_benchmark.vcf.gz HG002_${BUILD}_1_22_v4.2.1_benchmark.vcf.gz.tbi HG002_${BUILD}_1_22_v4.2.1_benchmark_noinconsistent.bed; do
  curl -fsSL -o "${OUT}/${f}" "${BASE}/${SAMPLE_DIR}/NISTv4.2.1/${BUILD}/${f}"
done
curl -fsSL -o "${OUT}/strat.tar.gz" "${BASE}/genome-stratifications/v3.6/genome-stratifications-${BUILD}@all.tar.gz"
tar -xzf "${OUT}/strat.tar.gz" -C "${OUT}"
```

Validate: `bcftools view -H truth/*.vcf.gz | wc -l` (v4.2.1 HG002 GRCh38 ≈ 4.0 M records) and `ls truth/${BUILD}@all/*.tsv`.

Other samples: `ChineseTrio/HG005_NA24631_son`, `ChineseTrio/HG006_NA24694_father`, `ChineseTrio/HG007_NA24695_mother`, `NA12878_HG001`, `AshkenazimTrio/HG003_NA24149_father`, `AshkenazimTrio/HG004_NA24143_mother`; same `NISTv4.2.1/<build>/` layout. v5.0q: `${SAMPLE_DIR}/v5.0q/HG002_${BUILD}_v5.0q_{smvar,stvar}.{vcf.gz,benchmark.bed}` (builds GRCh37, GRCh38, CHM13v2.0). CMRG: `${SAMPLE_DIR}/CMRG_v1.00/GRCh38/{SmallVariant,StructuralVariant}/`. SV v0.6: `${SAMPLE_DIR}/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.{vcf.gz,bed}` (GRCh37 only). SEQC2: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/` (`high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`, `high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz`, `High-Confidence_Regions_v1.2.bed`). Download mechanics for references live in `.claude/skills/genomics/reference.md`.

2. Normalise the query the same way the truth is (split multiallelics, left-align, drop non-PASS, keep sample of interest). Do **not** decompose MNPs differently between the two files.

```bash
set -euo pipefail
REF=ref/GRCh38.fa
QUERY=calls/sample.vcf.gz
SAMPLE=HG002
bcftools view -f PASS,. -s "${SAMPLE}" -Ou "${QUERY}" \
  | bcftools norm -f "${REF}" -m -any -Oz -o query.norm.vcf.gz
tabix -f query.norm.vcf.gz
bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" "${REF}" query.norm.vcf.gz truth/*.vcf.gz truth/*.bed
```

Validate: `bcftools stats query.norm.vcf.gz | grep -E '^SN|^TSTV'` (Ti/Tv sanity before spending compute).

3. Small variants: hap.py with the vcfeval engine (handles representation differences; the default `xcmp` engine does not).

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TRUTH=truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
CONF=truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
STRAT=truth/GRCh38@all/GRCh38-all-stratifications.tsv
QUERY=query.norm.vcf.gz
SEX=male
THREADS=8
OUTDIR=happy_out
mkdir -p "${OUTDIR}"
docker run --rm -u "$(id -u):$(id -g)" -v "${PWD}:/data" -w /data jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py "${TRUTH}" "${QUERY}" -f "${CONF}" -r "${REF}" \
  --engine vcfeval --engine-vcfeval-template ref/GRCh38.sdf \
  --gender "${SEX}" --stratification "${STRAT}" --threads "${THREADS}" \
  --pass-only -o "${OUTDIR}/HG002"
```

Validate: `grep -E '^Type,|,ALL,|,PASS,' happy_out/HG002.summary.csv | cut -d, -f1-3,11-14` (columns: Type, Filter, TRUTH.TOTAL … 11 METRIC.Recall, 12 METRIC.Precision, 13 METRIC.Frac_NA, 14 METRIC.F1_Score). WES/panel: add `-T "${TARGETS}"`. Build the SDF once: `rtg format -o ref/GRCh38.sdf "${REF}"`. Without the SDF hap.py builds one per run (slow). Stratification TSV name differs by release; `ls truth/GRCh38@all/*.tsv`.

4. Small variants, rtg vcfeval directly (faster; gives ROC by QUAL/GQ; use when hap.py's Python 2 image is unavailable).

```bash
set -euo pipefail
SDF=ref/GRCh38.sdf
TRUTH=truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
CONF=truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
QUERY=query.norm.vcf.gz
rtg vcfeval -b "${TRUTH}" -c "${QUERY}" -e "${CONF}" -t "${SDF}" \
  --ref-overlap --vcf-score-field QUAL --output-mode split -o vcfeval_out
```

Validate: `cat vcfeval_out/summary.txt` (precision/recall/F at the best threshold and at `None`); `rtg rocplot vcfeval_out/snp_roc.tsv.gz --png roc.png` for the threshold sweep. `--squash-ploidy` = allele-only matching (report as such).

5. Structural variants: Truvari.

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TRUTH=truth/HG002_GRCh38_v5.0q_stvar.vcf.gz
CONF=truth/HG002_GRCh38_v5.0q_stvar.benchmark.bed
QUERY=sv_calls.vcf.gz
truvari bench -b "${TRUTH}" -c "${QUERY}" -f "${REF}" --includebed "${CONF}" \
  --passonly --pick ac -r 2000 -C 5000 --sizemin 50 --sizemax 50000 -o truvari_out
truvari refine -f "${REF}" --use-original-vcfs --regions truvari_out/candidate.refine.bed truvari_out
```

Validate: `jq '.precision,.recall,.f1,."TP-base",.FP,.FN' truvari_out/summary.json` then `jq '.precision,.recall,.f1' truvari_out/refine.variant_summary.json`. Report the refined numbers for v5.0q and CMRG SV (complex SVs, differing representations); report plain `bench` for v0.6 Tier1 (`-r 2000 -C 3000 --pctseq 0 --pctsize 0.7` was the community default there). Depth-only / symbolic-ALT callers (CNVnator, Delly DEL/DUP without sequence, CNVkit): add `--pctseq 0` and `--dup-to-ins` where appropriate; `--typeignore` when the caller labels DUP vs INS inconsistently. Filter `ALT="*"` from v5.0q truth before un-refined `bench`: `bcftools view -e 'ALT="*"' -Oz -o truth_noast.vcf.gz "${TRUTH}"`.

6. Somatic SNV/indel: SEQC2 with som.py (allele matching, no genotype).

```bash
set -euo pipefail
REF=ref/GRCh38.fa
TRUTH_SNV=truth/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
CONF=truth/High-Confidence_Regions_v1.2.bed
QUERY=somatic.pass.snv.vcf.gz
docker run --rm -u "$(id -u):$(id -g)" -v "${PWD}:/data" -w /data jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/som.py "${TRUTH_SNV}" "${QUERY}" -r "${REF}" -f "${CONF}" \
  --feature-table generic -o sompy_out/snv
```

Validate: `cut -d, -f2-9 sompy_out/snv.stats.csv`. Run SNV and indel separately (two truth VCFs). Stratify by VAF bins (truth carries `VAF`/`TVAF`; query `AF`) with `bcftools view -i 'INFO/VAF<0.1'` subsets of the truth; report ≥0.2, 0.1–0.2, 0.05–0.1 separately. Tumour purity of the SEQC2 cell line is ~100 %: for low-purity performance use the SEQC2 in-silico dilution designs or your own mixtures.

7. Fill the results table (one row per Type × Filter × key stratification) and store `summary.csv`/`summary.json`, the truth URLs, md5s and tool versions next to it.

```text
| Query (caller vX, ref, filter) | Truth (set, version, build, region) | Type | Strat | TP | FP | FN | Precision | Recall | F1 | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| DeepVariant 1.10.0 WGS 35x PASS | GIAB HG002 v4.2.1 GRCh38 noinconsistent | SNP | all | | | | | | | |
| DeepVariant 1.10.0 WGS 35x PASS | GIAB HG002 v4.2.1 GRCh38 noinconsistent | INDEL | all | | | | | | | |
| same | same | SNP | segdup / lowmappability / homopolymer>6 | | | | | | | |
```

## QC thresholds

Defaults for "is this result worth believing"; assay-specific pass marks are in Validation.

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| `ref_sanity.sh` exit | 0 | any mismatch invalidates the comparison |
| Truth records inside confident BED | ≥ 99 % | lower means wrong BED/VCF pairing |
| hap.py `METRIC.Frac_NA` (query calls outside confident region) | WGS ≤ 0.15; WES ≤ 0.6 | high Frac_NA with low recall → wrong BED or build; WES is expected to be mostly outside |
| `FP.gt` share of FP (SNP) | ≤ 30 % | mostly genotype errors → check ploidy/sex (`--gender`), low depth |
| `FP.al` share of FP (INDEL) | ≤ 20 % | allele mismatches → normalisation/representation issue, retry with vcfeval engine |
| Ti/Tv (query, PASS) | WGS 2.0–2.1; WES 2.8–3.3 | cheap proxy before benchmarking; GRCh37 WES ≈ 3.0; non-human differs; canonical values in `.claude/skills/genomics-qc/SKILL.md` |
| Het/hom-alt (WGS) | 1.4–2.0 | ancestry-dependent; far outside → contamination, wrong sample, ploidy |
| Truvari `TP-base`/`TP-comp` gap | ≤ 5 % | large gap → many-to-one matching; inspect `--pick` and `--multimatch` |
| SV calls < `--sizemin` | excluded by design | report size bins (50–100, 100–1k, 1k–10k, >10k) |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| Recall ≈ 0, Frac_NA ≈ 1 | `chr` prefix mismatch or build mismatch | `ref_sanity.sh`; `bcftools annotate --rename-chrs`; never liftover truth if truth exists for the build |
| Recall drops ~5 % in INDEL only | query not left-aligned / multiallelics unsplit; `xcmp` engine | `bcftools norm -f ref -m -any`; `--engine vcfeval` |
| Precision too good to be true | benchmarked against truth with `-R` restriction instead of `-f` confident region, or truth VCF given as query | use `-f`; order is `hap.py TRUTH QUERY` |
| WES F1 looks like WGS | forgot `-T targets.bed`; off-target reads counted | pass both `-f` and `-T` |
| chrX FP.gt inflated for HG002/3/5 | hap.py assumed diploid X | `--gender male`; v4.2.1 truth is chr1–22 only anyway, so restrict query |
| gVCF as query | `<NON_REF>` blocks counted | `bcftools convert --gvcf2vcf` or genotype first |
| Multi-sample query: wrong sample compared | hap.py takes the first sample | `bcftools view -s` before comparing |
| v5.0q numbers lower than v4.2.1 for the same calls | v5.0q includes harder regions and is draft | report both; cite region size |
| Truvari counts a DEL as FP against truth DEL at same locus | symbolic ALT, no sequence | `--pctseq 0`; check `SVLEN` sign/`END` |
| Truvari FN cluster in tandem repeats | representation differs | `truvari refine`, or `truvari collapse` the query first |
| Somatic: F1 fine at high VAF, poor overall | low-VAF truth variants dominate FN | stratify by VAF; report purity/coverage of the query |
| SEQC2 run on GRCh37 data | truth is GRCh38 | re-align or use nf-core/variantbenchmarking `--liftover test` with a chain, and say so |
| Results differ between hap.py and vcfeval | different matching semantics (GT-aware default both, but ROC thresholding and `--ref-overlap`) | report the engine; do not mix engines within one table |
| F1 reported for CNV | no genotype/breakpoint truth of that type | report per-event overlap recall/precision with the overlap rule stated |

## nf-core

`nf-core/variantbenchmarking` 1.6.1 (`https://nf-co.re/variantbenchmarking`): hap.py, rtg vcfeval, Truvari, SVanalyzer, wittyer, som.py, liftover, stratification, multi-caller report. Samplesheet header: `id,test_vcf,caller` (optional columns `subsample`, `liftover`, `missing_genotypes`, Truvari/wittyer per-row params).

```csv
id,test_vcf,caller
dv_hg002,calls/hg002.deepvariant.vcf.gz,deepvariant
gatk_hg002,calls/hg002.gatk.vcf.gz,gatk
```

```bash
set -euo pipefail
nextflow run nf-core/variantbenchmarking -r 1.6.1 -profile docker \
  --input samplesheet.csv --outdir bench_out \
  --analysis germline --variant_type small --method "happy,rtgtools" \
  --truth_id HG002 --truth_vcf truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --regions_bed truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed --high_conf true \
  --fasta ref/GRCh38.fa --stratification_tsv truth/GRCh38@all/GRCh38-all-stratifications.tsv \
  --stratification_bed truth/GRCh38@all --preprocesses "split_multiallelic,normalize,filter_contigs"
```

`--variant_type structural --method truvari` for SV; `--variant_type copynumber` (truvari, wittyer, rtgtools cnveval) for CNV; `--analysis somatic` for som.py; add `--targets_bed` for WES. Execution mechanics (`-profile`, `-resume`, institutional configs) are in `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

Minimum acceptable F1 (PASS calls, confident region, correct build). These are starting defaults from published GIAB/precisionFDA-style results, not pass/fail law; justify any deviation in the results table.

| Assay | SNP F1 | INDEL F1 | Rationale |
|---|---|---|---|
| Short-read WGS ≥ 30x, DeepVariant/DRAGEN/GATK | ≥ 0.995 | ≥ 0.99 | modern callers reach ≈0.999/0.995 on v4.2.1; lower means a pipeline fault, not caller limits |
| Short-read WES (within targets ∩ confident) | ≥ 0.985 | ≥ 0.95 | capture edges and GC bias; indel recall limited by coverage uniformity |
| Targeted panel ≥ 500x | ≥ 0.99 | ≥ 0.97 | few variants: also report counts, not only ratios |
| PacBio HiFi WGS ≥ 20x | ≥ 0.995 | ≥ 0.99 | HiFi indel accuracy ≈ short-read; use v5.0q to show gain in hard regions |
| ONT R10 simplex ≥ 30x | ≥ 0.995 | ≥ 0.90 | homopolymer indels remain the weak spot; stratify `homopolymer_gt6` |
| SV long-read, v0.6 Tier1 / v5.0q refined | — | F1 ≥ 0.90 | sniffles2/pbsv/cuteSV reach ≥0.9 on Tier1 |
| SV short-read, v0.6 Tier1 | — | F1 0.5–0.7 expected | recall-limited; precision should still be ≥ 0.85 |
| Somatic WGS SEQC2, tumour 50x / normal 30x | ≥ 0.90 (VAF ≥ 0.05) | ≥ 0.70 | indel truth is small and hard; report per-VAF |
| Somatic WES SEQC2 | ≥ 0.85 | ≥ 0.60 | fewer reads per locus; low-VAF recall drops |
| CNV ≥ 10 kb vs v5.0q DEL/DUP (`--pctseq 0`) | — | recall ≥ 0.8, precision ≥ 0.7 | depth-only tools have fuzzy breakpoints; state overlap rule |

Cheap proxies when no truth set applies: Ti/Tv and het/hom (table above); dbSNP/gnomAD known-site fraction ≥ 95 % for germline WGS; trio Mendelian violation rate ≤ 1 % of sites; replicate concordance ≥ 99 % SNP. For per-caller validation of the *calls themselves* see the producing skill (`.claude/skills/genomics-wgs/SKILL.md`, `genomics-somatic`, `genomics-longread`, `genomics-cnv`).

## Additional resources

- `.claude/skills/genomics/reference.md` — reference FASTA/SDF, GIAB FASTQ/BAM download and T2T/CHM13 sources; load when the reference or raw benchmark reads are missing.
- `.claude/skills/genomics-workflow/SKILL.md` — nf-core/Nextflow execution (`-profile`, `-resume`, SLURM, containers); load when running nf-core/variantbenchmarking at scale.
- `.claude/skills/genomics-wgs/SKILL.md`, `.claude/skills/genomics-longread/SKILL.md`, `.claude/skills/genomics-somatic/SKILL.md`, `.claude/skills/genomics-cnv/SKILL.md` — producing the query VCFs and caller-specific filters; load before re-running a caller to fix a benchmark failure.
- GA4GH best practices (Krusche et al. 2019, `https://doi.org/10.1038/s41587-019-0054-x`), GIAB stratifications v3.6 (`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/README.md`), v5.0q README (`.../HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md`) — load for definitions of TP/FP.gt/FP.al/UNK and stratification names.
