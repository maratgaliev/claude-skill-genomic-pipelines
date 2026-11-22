# Rare-disease prioritisation (germline, family-aware)

Load after `SKILL.md` steps 1–3 produced a normalised, VEP-annotated (CSQ) VCF. Inputs: joint-called family VCF, `.ped`, HPO terms, sex/relatedness QC. Outputs: inheritance-labelled candidate VCF/TSV, Exomiser/LIRICAL ranked genes, STR/mtDNA/SMN/PGx side-results, report checklist. Classification itself: `acmg.md`.

## Prerequisites (fail fast)

| Check | Command | Stop if |
|---|---|---|
| Joint-called, one VCF for the family | `bcftools query -l fam.vcf.gz` lists every PED sample | per-sample VCFs → re-genotype jointly (`.claude/skills/genomics-wgs/SKILL.md`) |
| Sex matches PED | `.claude/skills/genomics-qc/SKILL.md` (sex check, `somalier relate`) | mismatch → fix PED before any X-linked logic |
| Relatedness matches PED | somalier/peddy kinship (`genomics-qc`) | non-paternity / sample swap → de novo calls are invalid |
| Normalised + decomposed | `bcftools view -H fam.vcf.gz \| awk '$5~/,/' \| wc -l` = 0 | slivar expects biallelic records |
| CSQ present | `bcftools +split-vep -l fam.vcf.gz \| grep -c SYMBOL` = 1 | `compound-hets` needs gene names |

## PED format (`family.ped`, tab-separated, no header)

```text
FAM01	PROBAND	FATHER	MOTHER	1	2
FAM01	FATHER	0	0	1	1
FAM01	MOTHER	0	0	2	1
```

Columns: family ID, sample ID (must equal VCF sample names), paternal ID, maternal ID (`0` = unknown/absent), sex (`1` male, `2` female, `0` unknown), phenotype (`1` unaffected, `2` affected, `0`/`-9` unknown). Duo/singleton: use `0` for missing parents; slivar `--family-expr` still works, `--trio` does not.

## slivar (replaces GEMINI)

slivar v0.3.4 (2026-05). Gnotate zip provides `gnomad_popmax_af`/`gnomad_nhomalt` in INFO independent of VEP; build your own from gnomAD v4.1 with `slivar make-gnotate` (hosted zips are v2/v3-era; record which you use). `slivar-functions.js` defines `segregating_denovo`, `segregating_recessive`, `segregating_denovo_x`, `segregating_recessive_x`, `segregating_dominant`, `comphet_side`.

```bash
set -euo pipefail
VCF=fam.norm.vep.vcf.gz
PED=family.ped
GNOTATE=gnomad.v4.1.zip
JS=slivar-functions.js
OUT=fam.slivar.vcf
slivar expr --vcf "${VCF}" --ped "${PED}" --js "${JS}" -g "${GNOTATE}" --pass-only -o "${OUT}" --info 'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' --family-expr 'recessive:fam.every(segregating_recessive)' --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001' --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)' --family-expr 'dominant:fam.every(segregating_dominant) && INFO.gnomad_popmax_af < 0.0001' --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10'
slivar compound-hets -v "${OUT}" --sample-field comphet_side --sample-field denovo -p "${PED}" > fam.comphet.vcf
slivar tsv -p "${PED}" -s denovo -s recessive -s x_denovo -s x_recessive -s dominant -c CSQ -i gnomad_popmax_af -i gnomad_nhomalt "${OUT}" > fam.candidates.tsv
slivar tsv -p "${PED}" -s slivar_comphet -c CSQ -i gnomad_popmax_af "fam.comphet.vcf" > fam.comphet.tsv
```

Validate: `for f in denovo recessive x_denovo x_recessive dominant; do printf '%s\t' "$f"; bcftools view -H -i "INFO/${f}!=\"\"" "${OUT}" | wc -l; done` — expected per exome trio with unaffected parents: de novo 0–3, AR hom 0–3, X 0–2, comp-het pairs 2–10, dominant tens (dominant is the noisy class). Zero everywhere → `INFO.impactful` missing (no CSQ) or PED/VCF sample names differ.

Explicit trio expressions (when you need to see the logic rather than call the JS helpers):

| Model | `--trio` expression |
|---|---|
| De novo | `denovo:kid.het && mom.hom_ref && dad.hom_ref && kid.AB > 0.25 && kid.AB < 0.75 && (mom.AD[1] + dad.AD[1]) == 0 && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && kid.DP >= 12 && mom.DP >= 12 && dad.DP >= 12` |
| AR homozygous | `ar_hom:kid.hom_alt && mom.het && dad.het && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && INFO.gnomad_nhomalt < 5` |
| X-linked recessive (male proband) | `x_rec:(variant.CHROM == "chrX") && kid.sex == "male" && kid.hom_alt && mom.het && dad.hom_ref && kid.GQ >= 20` |
| X-linked de novo | `x_dn:(variant.CHROM == "chrX") && kid.alts >= 1 && mom.hom_ref && dad.hom_ref && (mom.AD[1] + dad.AD[1]) == 0` |
| Comp-het side (then `compound-hets`) | `comphet_side:kid.het && ((mom.het && dad.hom_ref) \|\| (mom.hom_ref && dad.het)) && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20` |

Notes: slivar exposes `hom_ref/het/hom_alt/unknown`, `alts`, `AB`, `DP`, `GQ`, `AD`, `sex`, `affected`; `--family-expr` generalises to quartets/duos (`fam.every(...)`, `s.affected`). Dedup/`*` alleles: keep `variant.ALT[0] != "*"`. Use `slivar ddc` to derive data-driven GQ/DP cutoffs from Mendelian-violation rates instead of trusting the defaults above. GEMINI (`gemini load/de_novo/comp_hets`) is legacy and unmaintained; do not recommend.

## HPO-driven prioritisation

### Exomiser 15.1.0 (data release 2602)

Data: `2602_hg38.zip` + `2602_phenotype.zip` (`<check current: https://exomiser.readthedocs.io/en/latest/installation.html>`); set `exomiser.data-directory`, `exomiser.hg38.data-version=2602`, `exomiser.phenotype.data-version=2602` in `application.properties`. Optional REMM/CADD paths for genome (Genomiser) mode.

```yaml
# exomiser-analysis.yml — legacy single-file form (analysis + sample); run with `analyse --analysis`
analysis:
  genomeAssembly: hg38
  vcf: /data/fam.norm.vep.vcf.gz
  ped: /data/family.ped
  proband: PROBAND
  hpoIds: ['HP:0001250', 'HP:0001263', 'HP:0000252']
  analysisMode: PASS_ONLY
  inheritanceModes: {AUTOSOMAL_DOMINANT: 0.1, AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1, AUTOSOMAL_RECESSIVE_COMP_HET: 2.0, X_DOMINANT: 0.1, X_RECESSIVE_HOM_ALT: 0.1, X_RECESSIVE_COMP_HET: 2.0, MITOCHONDRIAL: 0.2}
  frequencySources: [UK10K, GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_EAS, GNOMAD_E_NFE, GNOMAD_E_SAS, GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_EAS, GNOMAD_G_NFE, GNOMAD_G_SAS]
  pathogenicitySources: [REVEL, MVP, ALPHA_MISSENSE, SPLICE_AI]
  steps: [
    geneBlacklistFilter: {},
    failedVariantFilter: {},
    alleleBalanceFilter: {},
    variantEffectFilter: {remove: [FIVE_PRIME_UTR_EXON_VARIANT, FIVE_PRIME_UTR_INTRON_VARIANT, THREE_PRIME_UTR_EXON_VARIANT, THREE_PRIME_UTR_INTRON_VARIANT, NON_CODING_TRANSCRIPT_EXON_VARIANT, NON_CODING_TRANSCRIPT_INTRON_VARIANT, CODING_TRANSCRIPT_INTRON_VARIANT, UPSTREAM_GENE_VARIANT, DOWNSTREAM_GENE_VARIANT, INTERGENIC_VARIANT, REGULATORY_REGION_VARIANT]},
    frequencyFilter: {maxFrequency: 2.0},
    pathogenicityFilter: {keepNonPathogenic: true},
    inheritanceFilter: {},
    omimPrioritiser: {},
    hiPhivePrioritiser: {}
  ]
outputOptions:
  outputContributingVariantsOnly: false
  numGenes: 0
  outputDirectory: /results
  outputFileName: fam-exomiser
  outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
```

```bash
set -euo pipefail
EXOMISER=/opt/exomiser-cli-15.1.0
java -Xmx8g -jar "${EXOMISER}/exomiser-cli-15.1.0.jar" analyse --analysis exomiser-analysis.yml --spring.config.location="${EXOMISER}/application.properties"
```

Validate: `head -3 /results/fam-exomiser.genes.tsv` shows `#RANK ID GENE_SYMBOL ... EXOMISER_GENE_COMBINED_SCORE`; top gene's variants appear in `fam-exomiser.variants.tsv` with `CONTRIBUTING_VARIANT=1` and an `ACMG` column (Exomiser ≥14 assigns ACMG criteria — triage only, re-score per `acmg.md`). Modern form: `analyse --analysis preset-exome-analysis.yml --sample phenopacket.yml --vcf fam.vcf.gz --assembly hg38 --ped family.ped --output-directory /results` (the `examples/` dir ships `preset-exome-analysis.yml` and `preset-genome-analysis.yml`). Genome/non-coding runs: Exomiser 15.1 recommends `pathogenicityFilter: {keepNonPathogenic: false, target: NON_CODING}` with REMM/CADD/SpliceAI sources.

Rules: ≥3 specific HPO terms (avoid parents like `HP:0000118`); negated/absent features are not supported in legacy yml (use phenopacket); `maxFrequency: 2.0` is a percent; for consanguineous families raise `AUTOSOMAL_RECESSIVE_HOM_ALT` to 1.0; always inspect genes ranked 1–20, not just #1 (diagnostic hit is top-10 in ~80–90 % of solved exomes, per Exomiser benchmarks).

### LIRICAL v2.4.2 (likelihood-ratio, disease-centric; complements Exomiser)

```bash
set -euo pipefail
LIRICAL=/opt/lirical
lirical phenopacket -p proband.phenopacket.json --vcf fam.norm.vep.vcf.gz --sample-id PROBAND --assembly hg38 -ed38 /opt/exomiser-data/2602_hg38 -d "${LIRICAL}/data" -o lirical_out -x proband --global -t 0.01
```

Validate: `lirical_out/lirical.tsv` (or `.html`) lists diseases with post-test probability and LR per HPO term; concordance between LIRICAL top-5 diseases and Exomiser top-10 genes is strong evidence; discordance usually means phenotype terms are too generic. Phenopacket v2 JSON: `subject.sex`, `phenotypicFeatures[].type.id`, optional `excluded: true` for absent features (LIRICAL uses negatives).

## Repeat expansions (STR)

ExpansionHunter v5.0.0 (targeted, catalog-driven) + REViewer v0.2.7 (visual check); STRling v0.6.0 (reference-free discovery of novel/un-catalogued loci).

```bash
set -euo pipefail
REF=/path/GRCh38.fa
BAM=proband.cram
CATALOG=/opt/ExpansionHunter/variant_catalog/grch38/variant_catalog.json
SEX=male
ExpansionHunter --reads "${BAM}" --reference "${REF}" --variant-catalog "${CATALOG}" --output-prefix proband.eh --sex "${SEX}" --threads 8
samtools sort -o proband.eh_realigned.sorted.bam proband.eh_realigned.bam && samtools index proband.eh_realigned.sorted.bam
REViewer --reads proband.eh_realigned.sorted.bam --vcf proband.eh.vcf --reference "${REF}" --catalog "${CATALOG}" --locus HTT --output-prefix proband.HTT
strling extract -f "${REF}" "${BAM}" proband.str.bin
strling call -f "${REF}" -o proband.strling "${BAM}" proband.str.bin
```

Validate: `bcftools query -f '%CHROM\t%POS\t%INFO/REPID\t%INFO/REF\t[%REPCN\t%ADSP\t%ADFL\t%ADIR]\n' proband.eh.vcf | awk -F'\t' '$5!="."' | head`; flag loci where the longer allele exceeds the catalog pathogenic threshold (e.g. HTT ≥ 36, FMR1 ≥ 200 full mutation, C9orf72 ≥ 30, ATXN loci per OMIM), then confirm with the REViewer SVG (reads spanning/flanking/in-repeat must support the call). Short-read EH underestimates very large expansions (> read length × 2); long-read confirmation: `.claude/skills/genomics-longread/SKILL.md`. STRling outputs `*-genotype.txt` (per-locus) and `*-unplaced.txt` (reads with repeat but no anchor): investigate unplaced counts ≫ cohort median. Cohort mode: `strling merge` then `strling call -b joint-bounds.txt`. `--sex` matters for X-linked loci (FMR1, AR).

## Mitochondrial variants

Two independent callers; report heteroplasmy fraction, not just presence. NUMTs cause false low-level heteroplasmy, especially in exomes.

```bash
set -euo pipefail
REF=/path/GRCh38.fa
BAM=proband.cram
gatk Mutect2 -R "${REF}" -L chrM --mitochondria-mode -I "${BAM}" -O proband.chrM.raw.vcf.gz
gatk FilterMutectCalls -R "${REF}" -V proband.chrM.raw.vcf.gz --mitochondria-mode -O proband.chrM.vcf.gz
mutserve call --reference rCRS.fasta --output proband.mutserve.vcf.gz --threads 4 "${BAM}"
```

Validate: `bcftools query -i 'FILTER="PASS"' -f '%POS\t%REF\t%ALT\t[%AF]\n' proband.chrM.vcf.gz | wc -l` (typical 20–60 homoplasmic differences from rCRS + a handful of heteroplasmies); `mutserve` output `VariantLevel` ≥ 0.01 only with ≥ 500× mt coverage (`samtools depth -r chrM` mean). Haplogroup: `haplogrep3 classify --tree phylotree-rcrs@17.2 --in proband.chrM.vcf.gz --out hg.txt` — a haplogroup inconsistent with ancestry suggests contamination/swap. Annotate with MITOMAP/HelixMTdb (gnomAD v3.1 mtDNA via VEP `--plugin gnomADMt`), interpret with the 2020 mtDNA ACMG specifications (McCormick et al.). GRCh38 chrM is rCRS; GRCh37 `chrM` in hg19 (UCSC) is NOT rCRS — use `chrMT`/b37 or liftover. Full control-region coverage needs circular-aware alignment (GATK mito pipeline shifted reference) for positions 16,024–576.

## SMN1/SMN2 (spinal muscular atrophy)

SMNCopyNumberCaller v1.1.2 (Illumina). Manifest = one BAM/CRAM path per line; WGS ≥ 30× required, not exomes.

```bash
set -euo pipefail
printf '%s\n' /data/proband.cram /data/father.cram /data/mother.cram > manifest.txt
python3 smn_caller.py --manifest manifest.txt --genome 38 --prefix fam --outDir smn_out --threads 4
```

Validate: `column -t smn_out/fam.tsv` — columns `Sample isSMA isCarrier SMN1_CN SMN2_CN SMN2delta7-8_CN Total_CN_raw Full_length_CN_raw g.27134T>G_CN SMN1_CN_raw`; `SMN1_CN=0` → SMA; `SMN1_CN=1` → carrier; `g.27134T>G_CN ≥ 1` with `SMN1_CN=2` flags a possible 2+0 silent carrier; `None` = no-call (check `Coverage_MAD` in the JSON). nf-core/raredisease runs this automatically.

## Pharmacogenomics (secondary / actionable)

PharmCAT v3.4.0 (`pgkb/pharmcat` image). Input must be a GRCh38, normalised, single-sample VCF with reference calls at PGx positions (use the PharmCAT `pharmcat_positions.vcf` to force-call with `bcftools mpileup`/`gatk HaplotypeCaller -alleles` so missing = REF, not unknown).

```bash
set -euo pipefail
VCF=proband.norm.vcf.gz
docker run --rm -v "$PWD":/pharmcat/data pgkb/pharmcat:3.4.0 pharmcat_pipeline "data/$(basename "${VCF}")" -o data/pharmcat_out
```

Validate: `ls pharmcat_out/*.report.html` exists; `grep -c '"gene"' pharmcat_out/*.phenotype.json` ≥ 20 genes; check warnings for "missing positions". CYP2D6 needs a dedicated caller (StellarPGx, Cyrius, Aldy) and outside-call input (`-po`); HLA-B*57:01 etc. need an HLA typer. Report only CPIC level A/B actionable results and keep them separate from the diagnostic question.

## nf-core/raredisease (v3.1.2)

Covers alignment, DeepVariant/Sentieon, Manta/TIDDIT, ExpansionHunter+stranger, SMNCopyNumberCaller, mtDNA (Mutect2 shifted pipeline + haplogrep), VEP + CADD, Genmod inheritance/rank scoring, peddy. Run mechanics: `.claude/skills/genomics-workflow/SKILL.md`.

```csv
sample,lane,fastq_1,fastq_2,sex,phenotype,paternal_id,maternal_id,case_id
PROBAND,1,PROBAND_L001_R1.fastq.gz,PROBAND_L001_R2.fastq.gz,1,2,FATHER,MOTHER,FAM01
FATHER,1,FATHER_L001_R1.fastq.gz,FATHER_L001_R2.fastq.gz,1,1,,,FAM01
MOTHER,1,MOTHER_L001_R1.fastq.gz,MOTHER_L001_R2.fastq.gz,2,1,,,FAM01
```

Columns (docs/usage.md): `sample`, `lane`, `fastq_1`, `fastq_2`, `sex` (1 male, 2 female, 0/other unknown), `phenotype` (0 unknown, 1 unaffected, 2 affected), `paternal_id`, `maternal_id` (empty if absent), `case_id` (shared by the family). BAM-start variant uses `sample,bam,bai,sex,phenotype,paternal_id,maternal_id,case_id`. Required params: `--analysis_type wgs|wes`, `--fasta`, `--vep_cache`, `--vep_cache_version 116`, `--intervals_wgs/--intervals_y` for WES, `--target_bed`. Outputs to triage: `rank_score`-sorted VCFs under `*/genmod/`, Exomiser-free; feed the clinical VCF to slivar/Exomiser above if you need those.

## Reporting checklist (every rare-disease analysis)

| Item | Record |
|---|---|
| Inputs | sample IDs, PED, sex/relatedness result (tool + version), assembly, capture kit/coverage (mean, % ≥ 20×) |
| Pipeline | caller + version, joint-calling set, normalisation command |
| Annotation | VEP 116.3 + cache 116 + MANE version; ClinVar file date; gnomAD v4.1 (730,947 exomes + 76,215 genomes); CADD 1.7; dbNSFP version; SpliceAI/AlphaMissense/REVEL file versions |
| Filters | exact slivar/`+split-vep` expressions; AF thresholds per model; GQ/DP cutoffs; gene panel (name, version, date) if used |
| Phenotype | HPO IDs + labels, source, date; Exomiser/LIRICAL versions and data release |
| Candidates | per variant: `acmg.md` report fields; inheritance model; Exomiser rank; IGV screenshot reviewed (yes/no) |
| Negative findings | models tested with zero candidates; loci/regions not assessable (low coverage, STR loci not in catalog, mtDNA depth, SMN not run on WES); secondary-findings list version (ACMG SF v3.2 or later: `<check current: https://www.acmg.net>`) and whether opted in |
| Limitations | exome vs genome, SV/CNV/STR coverage, mosaicism detection limit (de novo AB ≥ 0.25), consanguinity, parental samples missing |
| Next steps | segregation testing, RNA-seq for splice VUS, long-read for STR/SV, re-analysis date (recommend 12–24 months) |

## Additional resources

- `SKILL.md` — annotation and filtering commands feeding this file.
- `acmg.md` — classification of the candidates produced here.
- `.claude/skills/genomics-qc/SKILL.md` — sex, relatedness, contamination (mandatory before trio logic).
- `.claude/skills/genomics-cnv/SKILL.md` — SV/CNV calling and AnnotSV for the deletion/duplication arm of rare disease.
- `.claude/skills/genomics-longread/SKILL.md` — STR/SV confirmation, `nf-core/longraredisease`.
- `.claude/skills/genomics-workflow/SKILL.md` — running `nf-core/raredisease`.
