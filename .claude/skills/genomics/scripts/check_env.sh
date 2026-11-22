#!/usr/bin/env bash
# Print "tool<TAB>version|MISSING" for each tool named on the command line.
# Always exits 0 so it is safe inside skill dynamic-context injection.
# Each probe is time-limited and its output truncated to one short line, so a hung or
# hostile binary cannot stall skill loading or inject multi-line text into the agent context.
# Usage: check_env.sh gatk bwa-mem2 samtools bcftools nextflow

# Refuse to probe through relative PATH entries ("." or empty): they let a file in the
# current directory impersonate a tool at skill-load time.
case ":${PATH}:" in
  *:.:*|*::*|*:./*) echo "check_env: PATH contains a relative entry; refusing to probe tools" >&2; exit 0 ;;
esac

# Run a probe with a wall-clock limit when timeout(1)/gtimeout exist; otherwise run it directly.
if command -v timeout >/dev/null 2>&1; then
  probe() { timeout 10 "$@"; }
elif command -v gtimeout >/dev/null 2>&1; then
  probe() { gtimeout 10 "$@"; }
else
  probe() { "$@"; }
fi

bound() { head -c 80 | tr -d '\n\r\t' ; }

version_of() {
  local tool=$1 out
  case "$tool" in
    gatk)            out=$(gatk --version 2>/dev/null | grep -m1 -oE 'v?[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    picard)          out=$(picard MarkDuplicates --version 2>&1 | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+') ;;
    bwa)             out=$(bwa 2>&1 | grep -m1 -oE 'Version: [0-9.a-z-]+' | cut -d' ' -f2) ;;
    bwa-mem2)        out=$(bwa-mem2 version 2>/dev/null | tail -n1) ;;
    STAR)            out=$(STAR --version 2>/dev/null) ;;
    salmon)          out=$(salmon --version 2>&1 | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+') ;;
    kallisto)        out=$(kallisto version 2>/dev/null | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+') ;;
    bowtie2)         out=$(bowtie2 --version 2>/dev/null | grep -m1 -oE 'version [0-9.]+' | cut -d' ' -f2) ;;
    minimap2)        out=$(minimap2 --version 2>/dev/null) ;;
    samtools|bcftools|bedtools)
                     out=$("$tool" --version 2>/dev/null | head -n1 | awk '{print $2}') ;;
    tabix|bgzip)     out=$("$tool" --version 2>/dev/null | head -n1 | awk '{print $NF}') ;;
    mosdepth)        out=$(mosdepth --version 2>/dev/null | awk '{print $2}') ;;
    fastqc)          out=$(fastqc --version 2>/dev/null | awk '{print $2}') ;;
    fastp)           out=$(fastp --version 2>&1 | awk '{print $2}') ;;
    multiqc)         out=$(multiqc --version 2>/dev/null | grep -m1 -oE '[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    macs3|macs2)     out=$("$tool" --version 2>&1 | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+[a-z0-9]*') ;;
    vep)             out=$(vep --help 2>/dev/null | grep -m1 -oE 'ensembl-vep *: *[0-9.]+' | grep -oE '[0-9.]+$') ;;
    snpEff)          out=$(snpEff -version 2>/dev/null | awk '{print $2}') ;;
    nextflow)        out=$(nextflow -v 2>/dev/null | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+') ;;
    snakemake)       out=$(snakemake --version 2>/dev/null) ;;
    docker|singularity|apptainer)
                     out=$("$tool" --version 2>/dev/null | grep -m1 -oE '[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    conda|mamba|micromamba)
                     out=$("$tool" --version 2>/dev/null | grep -m1 -oE '[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    Rscript)         out=$(Rscript --version 2>&1 | grep -m1 -oE '[0-9]+\.[0-9]+\.[0-9]+') ;;
    python|python3)  out=$("$tool" --version 2>&1 | awk '{print $2}') ;;
    java)            out=$(java -version 2>&1 | grep -m1 -oE '"[0-9._]+"' | tr -d '"') ;;
    *)               out=$("$tool" --version 2>&1 | grep -m1 -oE '[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
  esac
  printf '%s' "${out:-present}"
}

if [ $# -eq 0 ]; then
  echo "usage: check_env.sh <tool> [tool...]" >&2
  exit 0
fi

export -f version_of
for t in "$@"; do
  case "$t" in
    *[!A-Za-z0-9._+-]*) printf '%s\tINVALID\n' "$t"; continue ;;   # tool names only; no paths or shell metacharacters
  esac
  if command -v "$t" >/dev/null 2>&1; then
    v=$(probe bash -c 'version_of "$1"' _ "$t" 2>/dev/null | bound)
    printf '%s\t%s\n' "$t" "${v:-present}"
  else
    printf '%s\tMISSING\n' "$t"
  fi
done

# Scheduler / container hints that change how commands should be issued.
for hint in sbatch qsub bsub; do
  command -v "$hint" >/dev/null 2>&1 && printf 'scheduler\t%s\n' "$hint"
done
printf 'cpus\t%s\n' "$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo '?')"
if [ -r /proc/meminfo ]; then
  awk '/MemTotal/ {printf "mem_gb\t%d\n", $2/1024/1024}' /proc/meminfo
elif command -v sysctl >/dev/null 2>&1; then
  printf 'mem_gb\t%d\n' "$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1024 / 1024 / 1024 ))"
fi
exit 0
