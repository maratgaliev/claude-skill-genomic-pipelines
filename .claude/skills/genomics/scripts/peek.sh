#!/usr/bin/env bash
# Safe, bounded preview of genomic files so an agent never loads a whole file into context.
# Usage: peek.sh <file> [n_records]   (default n=5)
# Prints size, type-specific header summary, and the first n records. Never reads more than needed.
set -uo pipefail

f=${1:?usage: peek.sh <file> [n]}
n=${2:-5}
[ -e "$f" ] || { echo "peek: no such file: $f" >&2; exit 1; }
real=$(readlink -f "$f" 2>/dev/null || printf '%s' "$f")

size=$(du -hL "$f" 2>/dev/null | cut -f1)
echo "== $f  ($size)"

has() { command -v "$1" >/dev/null 2>&1; }
reader() {            # emit decompressed text stream
  case "$1" in
    *.gz|*.bgz) gzip -cd "$1" ;;
    *.bz2)      bzip2 -cd "$1" ;;
    *.zst)      zstd -cdq "$1" ;;
    *)          cat "$1" ;;
  esac
}
fq_quick() {          # read count estimate for FASTQ without full decompression
  local bytes lines total
  if [[ $1 == *.gz ]]; then
    lines=$(gzip -cd "$1" | head -n 40000 | wc -l | tr -d ' ')
    bytes=$(gzip -cd "$1" | head -n 40000 | wc -c | tr -d ' ')
    total=$(gzip -l "$real" 2>/dev/null | awk 'NR==2{print $2}')
    if [ "${lines:-0}" -ge 4 ] && [ "${bytes:-0}" -gt 0 ] && [ -n "${total:-}" ]; then
      echo "  ~$(( total / (bytes / (lines / 4)) )) reads (gzip-size estimate; unreliable >4 GB)"
    fi
  fi
}

case "$f" in
  *.bam|*.cram)
    has samtools || { echo "samtools not found; falling back to file(1)"; file "$f"; exit 0; }
    samtools quickcheck "$f" && echo "  quickcheck: OK" || echo "  quickcheck: FAILED"
    echo "-- header (@HD/@RG/@PG, first 3 @SQ, total @SQ)"
    samtools view -H "$f" | grep -E '^@(HD|RG)' | head -n 6
    samtools view -H "$f" | grep '^@SQ' | head -n 3
    echo "  @SQ count: $(samtools view -H "$f" | grep -c '^@SQ')"
    samtools view -H "$f" | grep '^@PG' | cut -f1-3 | head -n 5
    echo "-- first $n records (cut to 9 columns)"
    samtools view "$f" | head -n "$n" | cut -f1-9
    [ -e "$f.bai" ] || [ -e "$f.crai" ] || [ -e "${f%.*}.bai" ] || echo "  WARNING: no index"
    ;;
  *.vcf|*.vcf.gz|*.vcf.bgz|*.bcf)
    if has bcftools; then
      echo "-- header summary"
      bcftools view -h "$f" | grep -E '^##(fileformat|reference|source|GATKCommandLine|DeepVariant|contig=<ID=(chr)?(1|X|M|MT)[,>])' | cut -c1-200 | head -n 12
      echo "  contigs: $(bcftools view -h "$f" | grep -c '^##contig')  INFO: $(bcftools view -h "$f" | grep -c '^##INFO')  FORMAT: $(bcftools view -h "$f" | grep -c '^##FORMAT')"
      echo "  samples: $(bcftools query -l "$f" | wc -l | tr -d ' ') -> $(bcftools query -l "$f" | head -n 5 | tr '\n' ' ')"
      echo "-- first $n records (cut to 150 chars)"
      bcftools view -H "$f" | head -n "$n" | cut -c1-150
      [ -e "$f.tbi" ] || [ -e "$f.csi" ] || echo "  WARNING: no index"
    elif [[ $f == *.bcf ]]; then
      echo "  binary BCF; install bcftools to inspect (no text fallback: never dump BGZF bytes into context)"
    else
      reader "$f" | grep -v '^##' | head -n "$((n + 1))" | cut -c1-150
    fi
    ;;
  *.fastq|*.fq|*.fastq.gz|*.fq.gz)
    echo "-- first $n reads"
    reader "$f" | head -n "$((n * 4))" | cut -c1-120
    first=$(reader "$f" | head -n 1)
    len=$(reader "$f" | sed -n 2p | tr -d '\n' | wc -c | tr -d ' ')
    echo "  read length (first read): $len"
    case "$first" in
      *" 1:N:"*|*" 2:N:"*) echo "  Illumina CASAVA 1.8+ header; mate=$(echo "$first" | awk '{print substr($2,1,1)}')" ;;
      @*runid=*|@*flow_cell_id=*) echo "  ONT header" ;;
      @m*/ccs*) echo "  PacBio HiFi header" ;;
    esac
    fq_quick "$f"
    ;;
  *.fa|*.fasta|*.fna|*.fa.gz|*.fasta.gz|*.fna.gz)
    echo "-- sequence names (first 10)"
    reader "$f" | grep '^>' | head -n 10 | cut -c1-100
    if [ -e "$f.fai" ]; then
      echo "  .fai present: $(wc -l < "$f.fai" | tr -d ' ') contigs; first: $(head -n1 "$f.fai" | cut -f1,2)"
    else
      echo "  WARNING: no .fai index (samtools faidx)"
    fi
    ;;
  *.gtf|*.gff|*.gff3|*.gtf.gz|*.gff.gz|*.gff3.gz)
    echo "-- header + first $n features"
    reader "$f" | grep '^#' | head -n 5
    reader "$f" | grep -v '^#' | head -n "$n" | cut -c1-150
    echo "  contigs (first 10): $(reader "$f" | grep -v '^#' | head -n 200000 | cut -f1 | sort -u | head -n 10 | tr '\n' ' ')"
    ;;
  *.bed|*.bed.gz|*.narrowPeak|*.broadPeak|*.bedGraph|*.bedgraph)
    echo "-- first $n lines"
    reader "$f" | head -n "$n" | cut -c1-150
    echo "  lines: $(reader "$f" | wc -l | tr -d ' ')  contigs (first 10): $(reader "$f" | cut -f1 | sort -u | head -n 10 | tr '\n' ' ')"
    ;;
  *.bw|*.bigWig|*.bigwig)
    has bigWigInfo && bigWigInfo "$f" | head -n 12 || echo "  binary bigWig; install UCSC bigWigInfo for summary"
    ;;
  *.h5ad|*.h5|*.hdf5)
    has h5ls && h5ls "$f" | head -n 20 || echo "  HDF5 container; inspect with h5ls or in Python (anndata.read_h5ad(backed='r'))"
    ;;
  *.pod5|*.fast5)
    echo "  nanopore raw signal; inspect with 'pod5 inspect summary' / 'pod5 view --limit 5'"
    ;;
  *.csv|*.tsv|*.txt|*.csv.gz|*.tsv.gz|*.txt.gz|*.sf|*.tab)
    echo "-- first $n lines"
    reader "$f" | head -n "$n" | cut -c1-200
    echo "  lines: $(reader "$f" | wc -l | tr -d ' ')"
    ;;
  *)
    file "$f"
    echo "-- first $n lines"
    reader "$f" | head -n "$n" | cut -c1-200
    ;;
esac
exit 0
