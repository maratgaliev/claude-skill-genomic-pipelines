---
name: genomics-annotation
description: Germline variant annotation and interpretation — normalisation, VEP/SnpEff/ANNOVAR, population frequency and in-silico scores, split-vep filtering, liftover, ACMG/AMP classification and rare-disease prioritisation. Use when a user has a VCF to annotate, asks about VEP/CSQ fields, gnomAD AF, CADD/REVEL/SpliceAI/AlphaMissense, ClinVar, ACMG criteria, trio/de novo filtering, Exomiser or slivar. Not for somatic tiering (genomics-somatic), CNV/SV calling (genomics-cnv) or variant calling itself (genomics-wgs).
argument-hint: "[annotate|filter|acmg|rare-disease] [vcf]"
paths:
  - "**/*.vep.vcf*"
  - "**/*.annotated.vcf*"
  - "**/*.ped"
  - "**/exomiser*.yml"
---

# Germline variant annotation and interpretation

Expert in germline annotation (VEP 116 / SnpEff / ANNOVAR), consequence-aware filtering with `bcftools +split-vep`, ACMG/AMP + ClinGen SVI classification and HPO-driven rare-disease prioritisation. Outputs are pinned commands, decision tables and validation checks, not narrative.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" bcftools tabix bgzip vep snpEff table_annovar.pl slivar CrossMap gatk java python3 Rscript || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version.

## Inspect first

| Read | Command | Infer |
|---|---|---|
| VCF header | `bcftools view -h in.vcf.gz \| grep -E '^##(reference|contig|source|GATKCommandLine|DeepVariant|INFO=<ID=(CSQ|ANN|BCSQ))' \| head -40` | assembly (contig lengths: chr1 248,956,422 = GRCh38; 249,250,621 = GRCh37), `chr` prefix, caller, prior annotation |
| Samples | `bcftools query -l in.vcf.gz \| wc -l` | single vs multi-sample; trio needs a `.ped` |
| Multi-allelics / norm state | `bcftools view -H in.vcf.gz \| awk '$5~/,/' \| head -3` | whether `bcftools norm` already ran |
| FILTER values | `bcftools query -f '%FILTER\n' in.vcf.gz \| sort \| uniq -c` | is `PASS` populated or all `.` (unfiltered gVCF-derived calls) |
| Existing CSQ keys | `bcftools +split-vep -l in.vcf.gz` | which fields are present; avoids re-running VEP |
| Reference match | `ref_sanity.sh ref.fa in.vcf.gz` | contig names/lengths vs cache assembly |

Ask only if still ambiguous: assembly when header lacks `##reference`/contig lengths; phenotype (HPO terms) and pedigree for rare-disease mode; intended use (research vs. accredited lab); gene panel / inheritance model; whether RefSeq (clinical reports) or Ensembl transcript IDs are required.

## Decision table

| Situation | Choice | Why |
|---|---|---|
| Default germline annotation | VEP 116 cache, `--mane --pick --pick_order mane_select,mane_plus_clinical,canonical,...` | one MANE Select consequence per variant; stable HGVS for reports |
| Need every transcript (splice, UTR, gene-level counts, ClinGen VCEP transcript specs) | VEP **without** `--pick`; filter later with `bcftools +split-vep -s worst` or by `MANE_SELECT` | `--pick` hides non-MANE transcripts a VCEP may mandate |
| RefSeq HGVS required | `--refseq` cache or `--merged` + `--hgvs --hgvsg`; still `--mane` | clinical reports cite NM_ accessions |
| Fast / no Perl | SnpEff 5.4 `GRCh38.mane.1.5.ensembl` (or `.refseq`) | MANE-only db, small, Java-only |
| Need gene-based tabular output, InterVar | ANNOVAR `table_annovar.pl` | InterVar depends on it |
| No cache download possible | `bcftools csq -g gff3 -f ref` | consequence only, haplotype-aware, no scores |
| GRCh37 input, GRCh38 resources | liftover first (step 7) | never mix assemblies |
| Trio / family | rare-disease.md (slivar, Exomiser) | inheritance-aware |
| Somatic VCF | `.claude/skills/genomics-somatic/SKILL.md` | AMP/ASCO/CAP tiering, not ACMG |

MANE Select policy: report on MANE Select (`MANE_SELECT` field non-empty); fall back to MANE Plus Clinical, then Ensembl canonical; record transcript version. Never silently change transcript between runs of the same sample.

## Workflow

### 1. Normalise

```bash
set -euo pipefail
REF=/path/GRCh38.fa
IN=input.vcf.gz
OUT=norm.vcf.gz
bcftools norm -m -both -f "${REF}" --check-ref w "${IN}" -Ou | bcftools norm -d exact -Oz -o "${OUT}" - && tabix -f -p vcf "${OUT}"
```

Validate: `bcftools view -H "${OUT}" | awk '$5~/,/' | wc -l` must print `0`; `--check-ref w` must report `0` REF mismatches on stderr (non-zero → wrong assembly).

Multi-sample warning: splitting multi-allelics rewrites per-sample GT/AD/PL and introduces `*` spanning-deletion alleles; `-d exact` on a merged cohort can drop records that differ only in genotypes. Normalise the joint-called VCF once (after `GenotypeGVCFs`/`glnexus`), never per-sample before merging. Left-alignment shifts POS, so re-index and re-run `ref_sanity.sh` afterwards.

### 2. VEP (Docker, cache 116, GRCh38)

Comments live above the block, never inside the continued command. Cache: `ensemblorg/ensembl-vep:release_116.3` with `homo_sapiens_vep_116_GRCh38` (install via `INSTALL.pl -a cf -s homo_sapiens -y GRCh38` inside the container; download mechanics live in `.claude/skills/genomics/reference.md`). dbNSFP 5.x: VEP-ready tabix file from dbnsfp.org (v5.4a released 2026-08; `<check current: https://www.dbnsfp.org/download>` for the exact `dbNSFP5.4a_grch38.gz` name). CADD v1.7: `whole_genome_SNVs.tsv.gz` + `gnomad.genomes.r4.0.indel.tsv.gz`. LOEUF: Ensembl-processed `loeuf_dataset_grch38.tsv.gz` (gnomAD v4 constraint). ClinVar via `--custom` (no `ClinVar` plugin exists). `gnomADc` is **coverage**, not constraint.

```bash
set -euo pipefail
VEP_IMG=ensemblorg/ensembl-vep:release_116.3
CACHE=/path/vep_cache
PLUG=/path/vep_plugins
DATA=/path/anno_db
REF=/path/GRCh38.fa
IN=norm.vcf.gz
OUT=vep.vcf.gz
THREADS=8
docker run --rm -u "$(id -u)" -v "$PWD":/work -v "${CACHE}":/cache -v "${PLUG}":/plugins -v "${DATA}":/db -v "$(dirname "${REF}")":/ref -w /work "${VEP_IMG}" vep -i "${IN}" -o "${OUT}" --vcf --compress_output bgzip --force_overwrite --fork "${THREADS}" --offline --cache --dir_cache /cache --cache_version 116 --assembly GRCh38 --fasta "/ref/$(basename "${REF}")" --dir_plugins /plugins --mane --pick --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length --symbol --hgvs --hgvsg --numbers --protein --biotype --canonical --domains --variant_class --sift b --polyphen b --af_gnomade --af_gnomadg --max_af --pubmed --no_stats --custom file=/db/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN --plugin CADD,snv=/db/whole_genome_SNVs.tsv.gz,indels=/db/gnomad.genomes.r4.0.indel.tsv.gz --plugin SpliceAI,snv=/db/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/db/spliceai_scores.raw.indel.hg38.vcf.gz --plugin AlphaMissense,file=/db/AlphaMissense_hg38.tsv.gz --plugin REVEL,file=/db/new_tabbed_revel_grch38.tsv.gz --plugin dbNSFP,/db/dbNSFP5.4a_grch38.gz,REVEL_score,gnomAD_LOEUF,gnomAD_lof.oe,MetaRNN_score,phyloP100way_vertebrate --plugin LOEUF,file=/db/loeuf_dataset_grch38.tsv.gz,match_by=gene --plugin LoFtool --plugin MaxEntScan,/plugins/maxentscan/fordownload
tabix -f -p vcf "${OUT}"
```

Validate: `bcftools +split-vep -l "${OUT}" | grep -c -x -E 'MANE_SELECT|ClinVar_CLNSIG|CADD_PHRED|SpliceAI_pred_DS_AG|am_pathogenicity|REVEL|LOEUF'` prints `7` (`-x` anchors whole lines so dbNSFP `REVEL_score`/`LOEUF`-like fields are not double-counted); `bcftools view -H "${OUT}" | wc -l` equals the input count.

Drop `--pick` (and add `--flag_pick` instead) when: ClinGen VCEP specs name a non-MANE transcript, splice/UTR analysis, per-transcript NMD/LoF review, or gene-burden counting. Plugin data files are licensed separately (dbNSFP academic vs commercial branch, CADD, REVEL); record each file's version.

### 3. Filter with split-vep (never raw INFO regex)

```bash
set -euo pipefail
IN=vep.vcf.gz
OUT=rare_impactful.vcf.gz
bcftools +split-vep "${IN}" -c SYMBOL,IMPACT,Consequence,MANE_SELECT,gnomADe_AF:Float,gnomADg_AF:Float,MAX_AF:Float,CADD_PHRED:Float,REVEL:Float,am_pathogenicity:Float,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,ClinVar_CLNSIG -s worst -i 'FILTER="PASS" && (MAX_AF<0.001 || MAX_AF=".") && (IMPACT="HIGH" || IMPACT="MODERATE" || SpliceAI_pred_DS_AG>=0.2 || SpliceAI_pred_DS_AL>=0.2 || SpliceAI_pred_DS_DG>=0.2 || SpliceAI_pred_DS_DL>=0.2) && ClinVar_CLNSIG!~"[Bb]enign"' -Oz -o "${OUT}" && tabix -f -p vcf "${OUT}"
```

Validate: `bcftools view -H "${OUT}" | wc -l` (expect ~100–500 per exome, ~1–3k per genome before inheritance filtering); `bcftools +split-vep -l "${OUT}" >/dev/null` still parses.

Rules: `FILTER="PASS"` (the VCF column), never `INFO/FILTER`; typed columns (`:Float`) or comparisons silently fail; `-s worst` picks the most severe transcript if `--pick` was not used; missing AF is `.` (absent from gnomAD ≠ AF 0) — keep it. Keep `ClinVar_CLNSIG~"Pathogenic"` regardless of AF in a second pass (rescue list).

### 4. Tabular export for R/Python

```bash
set -euo pipefail
IN=rare_impactful.vcf.gz
OUT=variants.tsv
bcftools +split-vep "${IN}" -s worst -d -A tab -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%SYMBOL\t%Consequence\t%IMPACT\t%HGVSc\t%HGVSp\t%MANE_SELECT\t%MAX_AF\t%gnomADe_AF\t%gnomADg_AF\t%CADD_PHRED\t%REVEL\t%am_pathogenicity\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\t%ClinVar_CLNSIG\t%ClinVar_CLNREVSTAT[\t%GT\t%DP\t%AD]\n' -H | sed '1s/^# //; 1s/\[[0-9]*\]//g' > "${OUT}"
```

Validate: `head -1 "${OUT}" | tr '\t' '\n' | wc -l` equals `awk -F'\t' 'NR==2{print NF}' "${OUT}"`.

```r
library(dplyr); library(readr)
v <- read_tsv("variants.tsv", na = c(".", ""), show_col_types = FALSE)
v <- v %>% mutate(SpliceAI_max = pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL, na.rm = TRUE))
cand <- v %>% filter(is.na(MAX_AF) | MAX_AF < 0.001, IMPACT %in% c("HIGH", "MODERATE") | SpliceAI_max >= 0.2, !grepl("benign", ClinVar_CLNSIG, ignore.case = TRUE)) %>% arrange(desc(CADD_PHRED))
write_tsv(cand, "candidates.tsv")
```

```python
import pandas as pd
v = pd.read_csv("variants.tsv", sep="\t", na_values=["."], low_memory=False)
v["SpliceAI_max"] = v[[c for c in v.columns if c.startswith("SpliceAI_pred_DS_")]].max(axis=1)
m = (v["MAX_AF"].isna() | (v["MAX_AF"] < 0.001)) & (v["IMPACT"].isin(["HIGH", "MODERATE"]) | (v["SpliceAI_max"] >= 0.2)) & ~v["ClinVar_CLNSIG"].fillna("").str.contains("benign", case=False)
v[m].sort_values("CADD_PHRED", ascending=False).to_csv("candidates.tsv", sep="\t", index=False)
```

### 5. SnpEff alternative

```bash
set -euo pipefail
DB=GRCh38.mane.1.5.ensembl
IN=norm.vcf.gz
OUT=snpeff.vcf.gz
snpEff download -v "${DB}"
snpEff -Xmx8g ann -v -canon -hgvs -lof -noStats "${DB}" "${IN}" | bgzip -c > "${OUT}" && tabix -f -p vcf "${OUT}"
SnpSift annotate -info CLNSIG,CLNREVSTAT,CLNDN clinvar.vcf.gz "${OUT}" | bgzip -c > snpeff.clinvar.vcf.gz
```

Validate: `bcftools +split-vep -a ANN -l "${OUT}" | head -3` lists `Allele Annotation Annotation_Impact`; db names via `snpEff databases | grep -E '^GRCh38\.(mane|p14)'` (`GRCh38.mane.1.5.{ensembl,refseq}`, `GRCh38.p14` RefSeq are listed in snpEff 5.4 config). Filter with `bcftools +split-vep -a ANN -c Annotation_Impact -i 'Annotation_Impact="HIGH"'` or `SnpSift filter "(ANN[*].IMPACT = 'HIGH')"`.

### 6. ANNOVAR alternative (needed for InterVar)

```bash
set -euo pipefail
HUMANDB=/path/humandb
IN=norm.vcf.gz
for db in refGeneWithVer clinvar_20260829 gnomad41_exome gnomad41_genome dbnsfp54a avsnp151; do annotate_variation.pl -buildver hg38 -downdb -webfrom annovar "${db}" "${HUMANDB}"; done
table_annovar.pl "${IN}" "${HUMANDB}" -buildver hg38 -out annovar -remove -protocol refGeneWithVer,clinvar_20260829,gnomad41_exome,gnomad41_genome,dbnsfp54a,avsnp151 -operation g,f,f,f,f,f -nastring . -vcfinput -polish -thread 8
```

Validate: `head -1 annovar.hg38_multianno.txt | tr '\t' '\n' | grep -c -E 'gnomad41|CLNSIG|REVEL'` ≥ 3. Table names were listed at https://annovar.openbioinformatics.org/en/latest/user-guide/download/ on 2026-11-22; confirm with `annotate_variation.pl -downdb -webfrom annovar avdblist -buildver hg38 .` before relying on a date suffix. ANNOVAR refGene is RefSeq-based and **not** MANE-aware; its left-normalisation differs from bcftools — compare a handful of indels against VEP HGVS.

### 7. Liftover (GRCh37 → GRCh38) before annotation

```bash
set -euo pipefail
CHAIN=/path/GRCh37_to_GRCh38.chain.gz
REF38=/path/GRCh38.fa
IN=hg19.vcf.gz
CrossMap vcf "${CHAIN}" "${IN}" "${REF38}" crossmap38.vcf && bgzip -f crossmap38.vcf && tabix -f -p vcf crossmap38.vcf.gz
gatk LiftoverVcf -I "${IN}" -O gatk38.vcf.gz -CHAIN "${CHAIN}" -REJECT rejected.vcf.gz -R "${REF38}" --RECOVER_SWAPPED_REF_ALT true --WARN_ON_MISSING_CONTIG true
```

Validate: `echo "rejected: $(bcftools view -H rejected.vcf.gz | wc -l) of $(bcftools view -H "${IN}" | wc -l)"` (expect < 1 % for exome calls); re-run step 1 (`--check-ref w`) on the lifted VCF — REF/ALT swaps happen on strand-flipped regions. `CrossMap.py` is the legacy entry point (CrossMap < 0.7); chain-file sources are in `.claude/skills/genomics/reference.md`. Prefer re-calling from GRCh38 alignments over liftover whenever reads exist.

### 8. Interpret

Load `acmg.md` for classification; `rare-disease.md` for trio/HPO/STR/mtDNA/PGx. Record for every report: VEP/SnpEff version, cache version, ClinVar file date, gnomAD release (v4.1 = 730,947 exomes + 76,215 genomes), dbNSFP/CADD/SpliceAI/AlphaMissense/REVEL file versions, MANE release, filter expression.

## QC thresholds

| Metric | Default | Rationale / when to deviate |
|---|---|---|
| `--check-ref w` mismatches | 0 | any >0 → wrong assembly or contig naming; stop |
| Variants annotated / input | 100 % | VEP drops nothing; mismatch → `--check_ref` or contig naming (`chr1` vs `1`: cache tolerates, `--custom` files do not) |
| MAX_AF rare cutoff | <0.001 (AD), <0.01 (AR, per-allele) | loosen for founder populations / recessive carriers; BA1 is 0.05, BS1 disease-specific |
| gnomAD source | v4.1 exomes+genomes, use `MAX_AF` or popmax (`grpmax`) | continental-AF only hides population-specific common alleles |
| IMPACT retained | HIGH, MODERATE, + SpliceAI ≥0.2 | LOW/MODIFIER only with SpliceAI/ClinVar rescue |
| Calibrated scores | REVEL ≥0.644 / ≤0.290; CADD ≥25.3 / ≤22.7; AlphaMissense ≥0.564; SpliceAI ≥0.2 (recall) / ≥0.5 (precision) | ClinGen SVI / Pejaver 2022 (see acmg.md); never sum predictors |
| Genotype quality | GQ ≥20, DP ≥10 (de novo DP ≥12 in all three) | long-read / low-coverage genomes need recalibration |
| ClinVar review status | ≥2 stars (`criteria_provided,_multiple_submitters,_no_conflicts`) for automated acceptance | 0–1 star needs manual review; "Conflicting" is not benign |
| Candidate count per exome after step 3 | 100–500 | thousands → AF filter not applied (field name typo gives `.` for all) |

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| `bcftools view -i 'INFO/gnomAD_AF<0.001'` returns everything | gnomAD AF lives inside CSQ, not a top-level INFO tag | `bcftools +split-vep -c gnomADe_AF:Float -i ...` |
| `INFO/FILTER` expression errors / returns nothing | FILTER is a VCF column | `FILTER="PASS"` |
| VEP: `ERROR: Plugin ClinVar not found` | no such plugin | `--custom file=clinvar.vcf.gz,...,fields=CLNSIG%CLNREVSTAT%CLNDN` |
| VEP: `--plugin gnomADc` gives numbers that look like depth | gnomADc = coverage | constraint via `--plugin LOEUF` or dbNSFP `gnomAD_LOEUF` |
| VEP command fails with `unknown option` after a `#` line | comment inside `\` continuation | single command, comments above |
| All CADD/REVEL/SpliceAI fields `.` | GRCh37 data files with GRCh38 VCF, or `chr` prefix mismatch in tabix file | `tabix -l file.gz \| head -2` vs `bcftools view -h \| grep contig` |
| Multiple CSQ per variant breaks TSV | no `--pick` and no `-s worst` | `-s worst` or `-d` (one row per transcript) |
| Indel absent in ClinVar/gnomAD although known | not normalised / different left-alignment | step 1 on both query and database VCF |
| HGVSp uses Ensembl ENSP, lab wants NP_ | Ensembl cache | `--refseq` or `--merged` cache + `--hgvs` |
| Candidate count thousands per exome | FILTER all `.` (no VQSR/hard filter) or AF field typo | check `bcftools query -f '%FILTER\n' \| sort \| uniq -c`; `+split-vep -l` |
| `bcftools norm` warns `REF_MISMATCH` | wrong FASTA (hg19 vs hg38, or decoys) | match assembly; `ref_sanity.sh` |
| ANNOVAR gnomAD column names differ (`gnomad41_exome_AF`) | db-version-specific headers | read `head -1 *multianno.txt` before scripting |
| SnpEff `ANN` has 16 fields but `split-vep` parsing fails | `-a ANN` not given | `bcftools +split-vep -a ANN` |
| Common allele flagged pathogenic (HFE C282Y, GJB2 V37I) | BA1 exception list | see acmg.md BA1 exceptions |

## nf-core

- `nf-core/raredisease` — alignment → SNV/SV/STR/mtDNA calling → VEP → Genmod ranking; samplesheet header: `sample,lane,fastq_1,fastq_2,sex,phenotype,paternal_id,maternal_id,case_id` (see rare-disease.md).
- `nf-core/sarek` — germline calling with `--tools haplotypecaller,deepvariant --vep_cache ... ` emits VEP-annotated VCFs; samplesheet header: `patient,sample,lane,fastq_1,fastq_2`.
- `nf-core/longraredisease` — long-read rare disease; see `.claude/skills/genomics-longread/SKILL.md`.

Run mechanics (`-profile`, `-resume`, configs): `.claude/skills/genomics-workflow/SKILL.md`.

## Validation

- Round-trip a known variant: `bcftools view -r chr17:43092919 vep.vcf.gz | bcftools +split-vep -f '%SYMBOL %HGVSc %HGVSp %ClinVar_CLNSIG\n'` must print BRCA1 HGVS and ClinVar status consistent with ClinVar web.
- Consequence sanity: `bcftools +split-vep -s worst -f '%Consequence\n' vep.vcf.gz | sort | uniq -c | sort -k1,1nr | head` — missense ≈ synonymous ×1.0–1.3 in exomes; stop_gained ≈ 0.5–1 % of coding.
- Ti/Tv unchanged by annotation (`bcftools stats`), record counts identical pre/post.
- Positive controls: re-annotate ClinVar 2-star pathogenic subset (`bcftools view -i 'INFO/CLNREVSTAT~"multiple"' clinvar.vcf.gz`) and confirm >95 % classify P/LP by your rules.
- Truth-set benchmarking of the calls themselves: `.claude/skills/genomics-benchmark/SKILL.md`.

## Additional resources

- `acmg.md` — ACMG/AMP 2015 criteria without PP5/BP6, ClinGen SVI modifications (PVS1 tree, PM2_Supporting, BA1 exceptions), Pejaver 2022 thresholds, Tavtigian 2020 points, InterVar, VCEP pointer, conflict handling. Load when classifying variants.
- `rare-disease.md` — Exomiser/LIRICAL HPO prioritisation, slivar trio expressions, PED format, STR/mtDNA/SMN/PGx, `nf-core/raredisease`, reporting checklist. Load for family/trio or undiagnosed-disease work.
- `.claude/skills/genomics/reference.md` — downloads: VEP cache, ClinVar, gnomAD, CADD, dbNSFP, chain files.
- `.claude/skills/genomics-qc/SKILL.md` — sex check, relatedness, contamination before trusting inheritance filters.
- `.claude/skills/genomics-somatic/SKILL.md` — somatic annotation and AMP/ASCO/CAP tiering.
- `.claude/skills/genomics-cnv/SKILL.md` — CNV/SV annotation (AnnotSV) and ACMG CNV scoring.
- `.claude/skills/genomics-benchmark/SKILL.md` — GIAB truth sets for the underlying calls.
