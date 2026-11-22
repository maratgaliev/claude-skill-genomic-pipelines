# ACMG/AMP germline classification (2015 framework + ClinGen SVI updates)

Load after annotation (`SKILL.md` steps 1–4). Input: one variant with MANE Select consequence, gnomAD v4.1 AF/popmax/homozygote count, ClinVar status, calibrated scores, phenotype and inheritance. Output: criteria list with strengths, points total, class, and the evidence lines used.

> **Research use only.** Classifications produced here are not clinical results. Clinical reporting requires review by a qualified professional under an accredited workflow (CLIA/CAP, ISO 15189, or local equivalent), with the exact ClinVar/gnomAD/transcript-set versions recorded.

## Order of operations

1. Gene–disease validity first (ClinGen gene validity ≥ Moderate; GenCC). No validated gene–disease pair → stop at VUS, do not apply PVS1/PS1/PM1.
2. Check for a ClinGen VCEP specification for the gene (https://cspec.genome.network/cspec/ui/svi/). If one exists it overrides everything below.
3. Apply criteria in this order: frequency (BA1/BS1/PM2) → variant type (PVS1 tree) → prior reports (PS1/PM5) → segregation/de novo (PS2/PM6/PP1/BS4) → functional (PS3/BS3) → computational (PP3/BP4, one predictor only) → phenotype (PP4) → in-trans (PM3/BP2).
4. Sum points (Tavtigian 2020), map to class, record every criterion with its evidence line.

## Criteria (Richards 2015, with SVI amendments; PP5 and BP6 removed)

| Code | Strength | Definition | Default evidence source | SVI note |
|---|---|---|---|---|
| PVS1 | Very strong (tunable) | Null variant (nonsense, frameshift, canonical ±1/2 splice, initiation codon, single/multi-exon deletion) in a gene where LoF is an established mechanism | IMPACT=HIGH + LOEUF < 0.35 or ClinGen haploinsufficiency 3 | use the Abou Tayoun 2018 tree (below); strength VS/S/M/Supporting/none |
| PS1 | Strong | Same amino-acid change as an established pathogenic variant, different nucleotide | ClinVar P/LP ≥2★ at same codon with same HGVSp | confirm splicing is not the mechanism of the known variant |
| PS2 | Strong (tunable) | De novo, both parents confirmed (relatedness check: `genomics-qc`) | slivar `denovo` + maternity/paternity verified | SVI points: 0.5 (phenotype non-specific) → 2 per proband; multiple de novos sum to PS2_VeryStrong |
| PS3 | Strong (tunable) | Well-established functional assay shows damaging effect | literature, MaveDB | Brnich 2019: strength scales with assay validation (OddsPath); default Supporting unless ≥11 validated controls |
| PS4 | Strong | Prevalence in affected ≫ controls (OR > 5, CI excludes 1) | case–control; rare-variant counts | may be used as PS4_Moderate/Supporting for case counts in rare dominant disease |
| PM1 | Moderate | Mutational hotspot or critical functional domain without benign variation | VEP `DOMAINS`, UniProt, VCEP domain lists | do not combine with PP3 for the same domain argument in some VCEPs |
| PM2 | **Supporting** (SVI 2020) | Absent (or extremely rare) in gnomAD v4.1 | `MAX_AF` = `.`; AR: AF < 1e-4 and no homozygotes | never Moderate by default; check coverage (`gnomADc` coverage plugin) before calling "absent" |
| PM3 | Moderate (tunable) | AR: detected in trans with a pathogenic variant | phased comp-het or parental testing | SVI points table (0.5–2 per observation); phase unknown → halve |
| PM4 | Moderate | Protein length change (in-frame indel in non-repeat, stop-loss) | Consequence `inframe_*`, `stop_lost` | not for repeat regions (`rmsk`) |
| PM5 | Moderate | Novel missense at a residue with a different established pathogenic missense | ClinVar P/LP ≥2★ same codon, different HGVSp | downgrade to Supporting if only one prior P/LP or Grantham distance small |
| PM6 | Moderate | Assumed de novo, parents not both confirmed | trio without relatedness check | |
| PP1 | Supporting (tunable) | Co-segregation in affected family members | LOD or meioses count | Jarvik & Browning 2016: 3 meioses Supporting, 5 Moderate, 7 Strong |
| PP2 | Supporting | Missense in gene with low benign missense rate and missense as mechanism | gnomAD `mis.oe`/MOEUF (dbNSFP `gnomAD_MOEUF`) < ~0.6 | weak; many VCEPs drop it |
| PP3 | Supporting→Strong (calibrated) | Computational evidence supports deleterious effect | one calibrated predictor (table below) | Pejaver 2022; do not stack predictors; do not apply PP3 when PVS1 is applied for the same splice effect |
| PP4 | Supporting | Phenotype highly specific for single-gene disease | HPO match (Exomiser/LIRICAL) | strengthen only per VCEP |
| ~~PP5~~ | — | Reputable source reports pathogenic without evidence | — | **Removed** (Biesecker & Harrison 2018). Use the underlying evidence instead. |
| BA1 | Stand-alone | AF > 5 % in any gnomAD continental population (≥ 2,000 alleles) | `gnomADe_AF_*`, `gnomADg_AF_*` popmax | exception list below; use filtering AF (`faf95`) when available |
| BS1 | Strong | AF greater than expected for the disorder | disease-specific threshold (prevalence × penetrance × allelic heterogeneity; calculator: cardiodb.org/allelefrequencyapp) | default 0.1 % AD / 1 % AR only if no disease value |
| BS2 | Strong | Observed in healthy adult (hom for AR, het for AD full-penetrance, hemizygous X-linked) | gnomAD `nhomalt` | penetrance and age of onset must support it |
| BS3 | Strong | Functional study shows no damaging effect | literature | same Brnich scaling as PS3 |
| BS4 | Strong | Lack of segregation in affected members | pedigree | beware phenocopies, reduced penetrance |
| BP1 | Supporting | Missense in gene where truncating variants are the only mechanism | gene mechanism | |
| BP2 | Supporting | In trans with pathogenic variant in full-penetrance AD gene, or in cis with pathogenic variant | phasing | |
| BP3 | Supporting | In-frame indel in repetitive region without known function | `rmsk`, low-complexity | |
| BP4 | Supporting→Very strong (calibrated) | Computational evidence suggests no impact | one calibrated predictor | Pejaver 2022 allows BP4_Moderate/Strong; splice: SpliceAI ≤ 0.1 |
| BP5 | Supporting | Found in a case with an alternate molecular basis | solved-case data | |
| ~~BP6~~ | — | Reputable source reports benign | — | **Removed** (SVI 2018). |
| BP7 | Supporting | Synonymous/intronic (outside ±1/2) with no predicted splice impact, not conserved | SpliceAI ≤ 0.1 + phyloP < 2 (or MaxEntScan delta small) | SVI splicing 2023: BP7 may apply to deep-intronic beyond +7/−21 with SpliceAI ≤ 0.1 |

## PVS1 decision tree (Abou Tayoun 2018, summary)

| Variant type | Condition | PVS1 strength |
|---|---|---|
| Nonsense / frameshift | predicted NMD (PTC > 50–55 nt upstream of last exon–exon junction) and exon present in biologically relevant transcript(s) | Very strong |
| Nonsense / frameshift | NMD-escape, truncated region critical to function or removes > 10 % of protein | Strong |
| Nonsense / frameshift | NMD-escape, region not critical, removes < 10 % | Moderate |
| Nonsense / frameshift | exon absent from relevant transcripts, or alternative exon | none |
| Canonical ±1/2 splice | exon skipping/cryptic site → frameshift + NMD | Very strong |
| Canonical ±1/2 splice | in-frame skip of critical region, or > 10 % removed | Strong |
| Canonical ±1/2 splice | in-frame skip of non-critical < 10 % | Moderate |
| Single/multi-exon deletion | out-of-frame + NMD, or whole-gene deletion with LoF mechanism | Very strong (whole gene: Very strong only if LoF mechanism; else none) |
| Single/multi-exon deletion | in-frame, critical | Strong; non-critical → Moderate |
| Duplication (tandem, in situ) | out-of-frame | Strong (Moderate if tandemness unproven) |
| Initiation codon | alternative start in-frame downstream with ≥1 pathogenic variant upstream of it | Moderate (Supporting if no such variant) |

Prerequisite in every branch: LoF is an established disease mechanism for the gene (ClinGen haploinsufficiency score 3, or ≥3 unrelated P/LP LoF variants, LOEUF < 0.35 is supportive, not sufficient). Never add PP3 (splice) when PVS1 already captures the splice effect. VEP `LoF`/NMD plugins (`--plugin NMD`) help but do not replace exon-level review.

## PM2_Supporting and gnomAD v4.1 practice

- gnomAD v4.1: 730,947 exomes + 76,215 genomes (GRCh38). Use `grpmax` (popmax) AF and `faf95` (filtering AF) where available; `MAX_AF` from VEP covers gnomAD e/g + 1000G.
- "Absent" requires coverage: mean depth ≥ 20 at the site in gnomAD (coverage via VEP `--plugin gnomADc` or gnomAD coverage TSV). Low coverage → PM2 not applicable.
- AR disease: PM2_Supporting if AF < 1e-4 (or disease-specific) **and** no homozygotes; AD: AF below BS1 threshold / ~1e-5.
- Founder variants in bottlenecked populations (Finnish, Ashkenazi, Amish) can exceed PM2 thresholds; use the disease-specific BS1 calculator rather than a blanket rule.

## BA1 exceptions (ClinGen SVI, Ghosh 2018)

BA1 (AF > 5 %) is stand-alone benign **except** for the curated list of common pathogenic/low-penetrance alleles (9 variants at publication). Confirmed members: HFE NM_000410 c.845G>A p.Cys282Tyr; HFE c.187C>G p.His63Asp; GJB2 NM_004004 c.109G>A p.Val37Ile. Full current list: `<check current: https://clinicalgenome.org/working-groups/sequence-variant-interpretation/>` (SVI "BA1 exception list"). Also do not apply BA1 when the population with AF > 5 % has < 2,000 alleles, to non-continental sub-populations, or to bottlenecked founder groups.

## Calibrated computational evidence (Pejaver 2022 missense; SVI Splicing 2023)

Use exactly one missense predictor per variant (choose the lab's default, record it). Thresholds are posterior-calibrated; intervals outside both directions give no evidence.

| Predictor | BP4 Very strong | BP4 Strong | BP4 Moderate | BP4 Supporting | PP3 Supporting | PP3 Moderate | PP3 Strong |
|---|---|---|---|---|---|---|---|
| REVEL | ≤ 0.003 | ≤ 0.016 | ≤ 0.183 | ≤ 0.290 | ≥ 0.644 | ≥ 0.773 | ≥ 0.932 |
| CADD v1.6/1.7 PHRED | — | ≤ 0.15 | ≤ 17.3 | ≤ 22.7 | ≥ 25.3 | ≥ 28.1 | — |
| BayesDel, MetaRNN, others | — | see Pejaver 2022 Table 2 for the predictor you pre-declare; do not transcribe from memory | | | | | |
| AlphaMissense (developer bins, not SVI-calibrated) | — | — | — | likely benign ≤ 0.34 | likely pathogenic ≥ 0.564 | — | — |

- Default: REVEL. CADD only when REVEL unavailable (CADD scores all variant types but is weaker for missense). VEP `am_class` gives the AlphaMissense bin; treat as Supporting at most unless a VCEP calibrates it.
- SpliceAI (max of DS_AG/AL/DG/DL, masked or raw — record which): ≥ 0.2 → PP3 (Supporting) for predicted splice impact (high recall); ≥ 0.5 recommended operating point for prioritisation; ≤ 0.1 → BP4/BP7 (no splice impact). Canonical ±1/2 variants go through PVS1, not PP3.
- Conservation (phyloP100way) is not a stand-alone criterion; use only for BP7 eligibility.
- Never count SIFT/PolyPhen/MutationTaster votes; the 2015 "multiple lines" language is superseded.

## Points system (Tavtigian 2020, Bayesian)

| Strength | Pathogenic points | Benign points |
|---|---|---|
| Very strong | +8 (PVS1) | −8 (BA1 is stand-alone; BP4_VeryStrong −8) |
| Strong | +4 | −4 |
| Moderate | +2 | −2 |
| Supporting | +1 | −1 |

| Total | Class |
|---|---|
| ≥ 10 | Pathogenic |
| 6 – 9 | Likely pathogenic |
| 0 – 5 | VUS (report sub-tier: 0–1 "cold", 2–3 "tepid", 4–5 "warm/hot" — Fowler/ClinGen VUS tiers, optional) |
| −1 – −6 | Likely benign |
| ≤ −7 | Benign |

BA1 alone → Benign regardless of points (unless on the exception list). Points replace the 2015 combination rules; when a VCEP publishes its own rules, follow the VCEP.

## Conflict handling

| Conflict | Resolution |
|---|---|
| Pathogenic and benign evidence both present (e.g. PVS1 + BS1) | Re-check mechanism and transcript; LoF in gene with high LoF tolerance (LOEUF > 1) voids PVS1; if both survive, sum points — a mixed sum near 0 is VUS, never "likely". |
| ClinVar P/LP vs your rules → B/LB | Compare review status and date; ≥3★ (expert panel) wins unless new population data; 0–1★ is not evidence (PP5 removed). Record the discordance. |
| gnomAD AF disagrees between exomes and genomes (> 10×) | Suspect mapping/segdup artefact; check gnomAD flags (`lcr`, `segdup`, `non_par`), use the dataset with better coverage at the site. |
| Predictor disagreement (REVEL benign, AlphaMissense pathogenic) | Apply the pre-declared default predictor only; note the other in comments; never pick post hoc. |
| Phenotype mismatch with gene | PP4 not applicable; consider dual diagnosis, phenotypic expansion requires literature; stays VUS. |
| Multiple transcripts give different consequences | Classify on MANE Select unless VCEP names another; report the alternative in comments. |
| Variant classified P in one disease, B in another for the same gene | Classification is per gene–disease pair; report pair explicitly. |

## InterVar (automated first pass only)

Repo: https://github.com/WGLab/InterVar (latest release v2.2.1, 2021-08; master receives occasional commits). Requires ANNOVAR and its `intervardb/`.

```bash
set -euo pipefail
INTERVAR=/path/InterVar
HUMANDB=/path/humandb
IN=norm.vcf.gz
bcftools view -Ov -o input.vcf "${IN}"
python3 "${INTERVAR}/Intervar.py" -i input.vcf --input_type=VCF -o intervar_out -b hg38 -d "${HUMANDB}" -t "${INTERVAR}/intervardb" --table_annovar="$(command -v table_annovar.pl)" --convert2annovar="$(command -v convert2annovar.pl)" --annotate_variation="$(command -v annotate_variation.pl)"
```

Validate: `grep -c . intervar_out.hg38_multianno.txt.intervar` ≈ variant count; column ` InterVar: InterVar and Evidence ` is parseable.

Caveats (why it is a triage tool, not a classifier):
- Implements 2015 rules literally: applies PP5/BP6, PM2 at Moderate, PP3/BP4 by predictor vote counts (dbNSFP `dbnsfp42a`/`dbnsfp47a`), no points system, no VCEP specs, no PVS1 tree.
- Shipped `config.ini` pins old ANNOVAR tables (`dbnsfp42a clinvar_20210501 gnomad_genome ...`, `current_version = Intervar_20210727`); update `database_names` to `gnomad41_exome gnomad41_genome clinvar_<current> dbnsfp54a` only after checking the code's expected column headers, or AF/ClinVar rules silently fail. ANNOVAR's `intervar_20250721` (hg38) is the newer evidence table.
- Indels in `intervar_*` evidence tables unsupported; splice/BP7 logic absent.
- Use its evidence string as a checklist; re-score with the tables above. Alternatives with maintained rule sets: Exomiser ≥14 ACMG assignment (rare-disease.md), Franklin/VarSome (web, licence/consent constraints for patient data).

## ClinGen VCEP specifications

- Registry: https://cspec.genome.network/cspec/ui/svi/ — gene/disease-specific strengths, thresholds (e.g. BS1/PM2 AF cutoffs, PP3 predictor + cutoff, PVS1 exon maps, PS3 assay lists).
- When a VCEP exists, document its version in the report; it supersedes this file.
- Expert-panel ClinVar submissions (3★) carry the VCEP classification; a 3★ record is the strongest prior available.

## Report fields per variant

Gene (HGNC), MANE Select transcript with version, HGVSg/c/p, zygosity, inheritance, gnomAD v4.1 grpmax AF + nhomalt (+ coverage), ClinVar ID + review status + date, predictor name + score + strength, every criterion with strength and one-line evidence, points total, class, gene–disease pair, VCEP used, software and database versions.
