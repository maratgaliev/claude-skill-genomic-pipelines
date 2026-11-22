#!/usr/bin/env bash
# Reference sanity gate: verify that every input agrees with the reference on contig naming and length.
# Usage: ref_sanity.sh <reference.fa|.fa.gz> <file...>
#   files: BAM/CRAM/SAM, VCF/BCF, GTF/GFF, BED-like, .fai, .dict
# Exit 0 = consistent (WARN rows allowed); exit 1 = MISMATCH found; exit 2 = cannot evaluate (missing tool/index).
# Writes only inside its own temp directory; never modifies or indexes files in place.
set -uo pipefail
export LC_ALL=C   # deterministic sort/join regardless of the caller's locale

ref=${1:?usage: ref_sanity.sh <reference.fa> <file...>}
shift
[ $# -ge 1 ] || { echo "ref_sanity: give at least one input file" >&2; exit 2; }
[ -e "$ref" ] || { echo "ref_sanity: reference not found: $ref" >&2; exit 2; }

tmp=$(mktemp -d) || { echo "ref_sanity: cannot create a temp directory" >&2; exit 2; }
trap 'rm -rf "$tmp"' EXIT
have() { command -v "$1" >/dev/null 2>&1; }

# Locate or build (in tmp) the .fai
fai=""
for cand in "$ref.fai" "${ref%.gz}.fai" "${ref%.fa}.fai" "${ref%.fasta}.fai"; do
  [ -e "$cand" ] && { fai=$cand; break; }
done
if [ -z "$fai" ]; then
  have samtools || { echo "ref_sanity: no .fai next to $ref and samtools not found; run 'samtools faidx $ref' first" >&2; exit 2; }
  echo "ref_sanity: no .fai next to $ref; indexing into a temp file (run 'samtools faidx' yourself to keep it)" >&2
  fai="$tmp/ref.fai"
  case "$ref" in
    *.gz) samtools faidx --fai-idx "$fai" --gzi-idx "$tmp/ref.gzi" "$ref" ;;
    *)    samtools faidx --fai-idx "$fai" "$ref" ;;
  esac || { echo "ref_sanity: faidx failed (bgzip-compressed? plain gzip is not indexable)" >&2; exit 2; }
fi

cut -f1,2 "$fai" | sort -k1,1 > "$tmp/ref.tsv"
ref_n=$(wc -l < "$tmp/ref.tsv" | tr -d ' ')

# Naming style = majority of primary-looking contigs (ignores HLA-*, *_alt, *_decoy, *_random, chrUn, GL/KI/NC accessions).
style_of() {
  awk -F'\t' '
    $1 ~ /^(HLA-|GL|KI|KN|KQ|KV|KZ|ML|NC_|NT_|NW_|JH|MU)/ || $1 ~ /(_alt|_decoy|_random|_fix|_hap)/ || $1 ~ /^chrUn/ || $1 ~ /^Un/ {next}
    {t++; if ($1 ~ /^chr/) c++}
    END {if (t==0) print "unknown"; else if (c*2>t) print "chr"; else print "nochr"}' "$1"
}
guess_build() {   # arg: TSV name<TAB>len ; assembly guess from chr1/1 length
  local l
  l=$(awk -F'\t' '($1=="chr1"||$1=="1") && $2!="" {print $2; exit}' "$1")
  case "${l:-}" in
    248956422) echo GRCh38 ;;
    249250621) echo GRCh37/hg19 ;;
    248387328) echo T2T-CHM13v2.0 ;;
    195154279) echo GRCm39 ;;
    195471971) echo GRCm38/mm10 ;;
    "")        echo "unknown" ;;
    *)         echo "unknown(chr1=$l)" ;;
  esac
}
ref_style=$(style_of "$tmp/ref.tsv")
ref_build=$(guess_build "$tmp/ref.tsv")
echo "reference: $ref  contigs=$ref_n  style=$ref_style  build=$ref_build"

reader() { case "$1" in *.gz|*.bgz) gzip -cd "$1" ;; *) cat "$1" ;; esac; }
sq_awk='$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)} print n"\t"l}'
contig_awk='/^##contig/{n="";l="";for(i=1;i<=NF;i++){if($i=="ID")n=$(i+1);if($i=="length")l=$(i+1)} print n"\t"l} /^#CHROM/{exit}'
status=0

# contigs_of <file>  -> name<TAB>len rows on stdout; returns 2 when a required tool is missing
contigs_of() {
  local f=$1
  case "$f" in
    *.bam|*.cram|*.sam)
      have samtools || return 2
      samtools view -H "$f" | awk -F'\t' "$sq_awk" ;;
    *.bcf)
      have bcftools || return 2
      bcftools view -h "$f" | awk -F'[<>,=]' "$contig_awk" ;;
    *.vcf|*.vcf.gz|*.vcf.bgz)
      if have bcftools; then bcftools view -h "$f" | awk -F'[<>,=]' "$contig_awk"
      else reader "$f" | awk -F'[<>,=]' "$contig_awk"; fi ;;
    *.fai)  cut -f1,2 "$f" ;;
    *.dict) awk -F'\t' "$sq_awk" "$f" ;;
    *)      # BED/GTF/GFF/narrowPeak/TSV: first column, skipping comment/track/browser lines
      reader "$f" | grep -Ev '^(#|track|browser)' | head -n 2000000 | cut -f1 | sort -u | awk '{print $1"\t"}' ;;
  esac
}

for f in "$@"; do
  [ -e "$f" ] || { echo "MISSING  $f"; status=1; continue; }
  if ! contigs_of "$f" > "$tmp/raw.tsv"; then
    echo "SKIP     $f  (need samtools/bcftools to read this format; install or convert, then rerun)"
    status=2; continue
  fi
  awk -F'\t' '$1!=""' "$tmp/raw.tsv" | sort -k1,1 -u > "$tmp/in.tsv"
  n=$(wc -l < "$tmp/in.tsv" | tr -d ' ')
  if [ "$n" -eq 0 ]; then
    case "$f" in
      *.vcf|*.vcf.gz|*.vcf.bgz|*.bcf)   # header lacks ##contig: use body contigs
        { if have bcftools; then bcftools view -H "$f"; else reader "$f" | grep -v '^#'; fi; } \
          | head -n 2000000 | cut -f1 | sort -u | awk '{print $1"\t"}' > "$tmp/in.tsv"
        n=$(wc -l < "$tmp/in.tsv" | tr -d ' ') ;;
    esac
  fi
  if [ "$n" -eq 0 ]; then echo "EMPTY    $f  (no contigs found)"; status=1; continue; fi

  style=$(style_of "$tmp/in.tsv")
  build=$(guess_build "$tmp/in.tsv")
  join -t $'\t' -v1 "$tmp/in.tsv" "$tmp/ref.tsv" | cut -f1 > "$tmp/missing.txt"
  # Alt/HLA/decoy/unplaced contigs absent from a no-alt reference: harmless for BAM/CRAM headers, a WARN only.
  grep -Ev '^(HLA-|chrUn|Un)|(_alt|_decoy|_random|_fix|_hap)' "$tmp/missing.txt" > "$tmp/missing_primary.txt" || true
  missing_n=$(grep -c . "$tmp/missing.txt" || true)
  missing_primary_n=$(grep -c . "$tmp/missing_primary.txt" || true)
  lenmis=$(join -t $'\t' "$tmp/in.tsv" "$tmp/ref.tsv" | awk -F'\t' '$2!="" && $2!=$3 {print $1"("$2" vs "$3")"}')
  lenmis_n=$(printf '%s' "$lenmis" | grep -c . || true)

  verdict=OK; notes=""
  if [ "$style" != unknown ] && [ "$ref_style" != unknown ] && [ "$style" != "$ref_style" ]; then
    verdict=MISMATCH; notes+=" naming:$style-vs-$ref_style"
  fi
  if [ "$missing_primary_n" -gt 0 ]; then
    verdict=MISMATCH; notes+=" ${missing_primary_n}_contigs_not_in_ref($(head -n 3 "$tmp/missing_primary.txt" | tr '\n' ',')…)"
  elif [ "$missing_n" -gt 0 ]; then
    [ "$verdict" = OK ] && verdict=WARN
    notes+=" ${missing_n}_alt/decoy_contigs_not_in_ref"
  fi
  [ "$lenmis_n" -gt 0 ] && { verdict=MISMATCH; notes+=" ${lenmis_n}_length_mismatches($(printf '%s' "$lenmis" | head -n 2 | tr '\n' ',')…)"; }
  case "$build:$ref_build" in
    unknown*:*|*:unknown*) ;;
    *) [ "$build" != "$ref_build" ] && { verdict=MISMATCH; notes+=" build:$build"; } ;;
  esac
  [ "$verdict" = MISMATCH ] && status=1
  printf '%-8s %s  contigs=%s style=%s build=%s%s\n' "$verdict" "$f" "$n" "$style" "$build" "$notes"
done

if [ "$status" -eq 1 ]; then
  cat >&2 <<'EOF'
ref_sanity: at least one input disagrees with the reference.
  naming  -> bcftools annotate --rename-chrs map.txt (VCF) | samtools reheader (BAM) | sed on BED/GTF; never mix styles
  build   -> stop: realign/lift over (CrossMap, gatk LiftoverVcf) rather than renaming
  missing -> primary contigs absent from the reference: wrong reference or wrong file; alt/decoy-only gaps are reported as WARN
EOF
elif [ "$status" -eq 2 ]; then
  echo "ref_sanity: some inputs could not be evaluated (missing tools); do not treat them as checked" >&2
fi
exit "$status"
