# Workflow templates

Copy-ready scaffolds referenced from `SKILL.md`. Every template pins software, declares resources, stages through scratch and records versions. Replace `${UPPER_CASE}` placeholders; keep the structure.

## 1. Nextflow DSL2 module + workflow

Layout:

```text
mypipe/
├── main.nf
├── nextflow.config
├── conf/{base,slurm,test}.config
├── modules/local/samtools_stats/main.nf
└── assets/samplesheet.csv
```

`modules/local/samtools_stats/main.nf` — nf-core conventions: `tuple val(meta), path(...)` in/out, `task.cpus`/`task.memory`, `ext.args`/`ext.prefix`, `versions.yml`, `stub:` block.

```groovy
process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::samtools=1.25"
    container "quay.io/biocontainers/samtools:1.25--h9dcdb79_0"

    cpus   4
    memory { 4.GB * task.attempt }
    time   { 2.h * task.attempt }

    publishDir "${params.outdir}/samtools_stats", mode: 'copy', pattern: '*.stats'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools stats --threads ${task.cpus} ${args} --reference ${fasta} ${bam} > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    touch versions.yml
    """
}
```

Pin check before committing: `skopeo inspect docker://quay.io/biocontainers/samtools:1.25--h9dcdb79_0 | jq -r .Digest` and replace the tag with `@sha256:<digest>` for production. Verify the tag exists first (`https://quay.io/repository/biocontainers/samtools?tab=tags`); Nextflow rewrites `quay.io/...` to `docker://` for Singularity/Apptainer automatically.

`main.nf`:

```groovy
include { SAMTOOLS_STATS } from './modules/local/samtools_stats/main'

workflow {
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).first()

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, sex: row.sex ?: 'NA']
            tuple(meta, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true))
        }
        .set { ch_bam }

    SAMTOOLS_STATS(ch_bam, ch_fasta)

    SAMTOOLS_STATS.out.versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}
```

`assets/samplesheet.csv`:

```csv
sample,sex,bam,bai
HG002,male,/data/HG002.cram,/data/HG002.cram.crai
```

Run: `nextflow run . -profile docker --input assets/samplesheet.csv --fasta ref/GRCh38.fa --outdir results -resume`. Lint: `nextflow lint .` and `nf-core pipelines lint` when the repo follows the nf-core template.

## 2. `nextflow.config` with profiles

```groovy
params {
    input          = null
    fasta          = null
    outdir         = 'results'
    max_cpus       = 32
    max_memory     = '192.GB'
    max_time       = '72.h'
}

manifest {
    name            = 'myorg/mypipe'
    description     = 'Example genomics pipeline'
    version         = '0.1.0'
    nextflowVersion = '!>=25.04.0'
}

process {
    errorStrategy = { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'finish' }
    maxRetries    = 2
    cache         = 'lenient'

    withLabel: process_low    { cpus = 2;  memory = { 8.GB  * task.attempt }; time = { 4.h  * task.attempt } }
    withLabel: process_medium { cpus = 8;  memory = { 32.GB * task.attempt }; time = { 12.h * task.attempt } }
    withLabel: process_high   { cpus = 16; memory = { 96.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'SAMTOOLS_STATS' { ext.args = '--remove-overlaps' }
}

profiles {
    docker {
        docker.enabled      = true
        docker.runOptions   = '-u $(id -u):$(id -g)'
        singularity.enabled = false
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = System.getenv('NXF_SINGULARITY_CACHEDIR') ?: "${System.getenv('HOME')}/.nxf_singularity"
        singularity.pullTimeout = '2h'
        docker.enabled         = false
    }
    apptainer {
        apptainer.enabled    = true
        apptainer.autoMounts = true
        apptainer.cacheDir   = System.getenv('NXF_APPTAINER_CACHEDIR') ?: "${System.getenv('HOME')}/.nxf_apptainer"
    }
    conda {
        conda.enabled  = true
        conda.useMamba = true
        conda.cacheDir = System.getenv('NXF_CONDA_CACHEDIR') ?: "${System.getenv('HOME')}/.nxf_conda"
    }
    slurm {
        process.executor       = 'slurm'
        process.queue          = { task.time > 24.h ? 'long' : 'core' }
        process.clusterOptions = "--account=${System.getenv('SLURM_ACCOUNT') ?: 'myproj'}"
        process.scratch        = '$SLURM_TMPDIR'
        executor.queueSize     = 200
        executor.submitRateLimit = '20/1min'
        executor.pollInterval  = '30 sec'
    }
    wave {
        wave.enabled   = true
        wave.strategy  = ['conda', 'container']
        wave.freeze    = true
    }
    fusion {
        fusion.enabled = true
        wave.enabled   = true
    }
    test {
        params.input  = "${projectDir}/assets/samplesheet.test.csv"
        params.fasta  = "${projectDir}/assets/chr22.fa"
        params.max_cpus = 2
        params.max_memory = '6.GB'
    }
}

timeline { enabled = true; file = "${params.outdir}/pipeline_info/timeline.html"; overwrite = true }
report   { enabled = true; file = "${params.outdir}/pipeline_info/report.html";   overwrite = true }
trace    { enabled = true; file = "${params.outdir}/pipeline_info/trace.txt";     overwrite = true
           fields = 'task_id,hash,name,status,exit,attempt,realtime,cpus,%cpu,memory,peak_rss,rchar,wchar,workdir' }
dag      { enabled = true; file = "${params.outdir}/pipeline_info/dag.html";      overwrite = true }
```

Use: `-profile singularity,slurm`; add `-c conf/site.config` for anything institution-specific instead of editing this file. For nf-core pipelines prefer the institutional profile from `nf-core/configs` over a local copy.

## 3. Snakemake 9 `Snakefile` + rule + SLURM profile

`Snakefile`:

```python
import os
import pandas as pd

configfile: "config/config.yaml"

SAMPLES = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)

wildcard_constraints:
    sample="|".join(SAMPLES.index),

rule all:
    input:
        expand("results/stats/{sample}.stats", sample=SAMPLES.index),
        "results/multiqc/multiqc_report.html",

rule samtools_stats:
    input:
        bam=lambda wc: SAMPLES.loc[wc.sample, "bam"],
        bai=lambda wc: SAMPLES.loc[wc.sample, "bai"],
        fasta=config["fasta"],
    output:
        "results/stats/{sample}.stats",
    log:
        "logs/samtools_stats/{sample}.log",
    conda:
        "envs/samtools.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.25--h9dcdb79_0"
    threads: 4
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: 120 * attempt,
        tmpdir=lambda wc: os.environ.get("SLURM_TMPDIR", "/tmp"),
        slurm_partition="core",
    retries: 2
    shell:
        "samtools stats --threads {threads} --reference {input.fasta} {input.bam} > {output} 2> {log}"

rule multiqc:
    input:
        expand("results/stats/{sample}.stats", sample=SAMPLES.index),
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "envs/multiqc.yaml"
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "multiqc --force --outdir results/multiqc results/stats > {log} 2>&1"
```

`envs/samtools.yaml`:

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - samtools=1.25
```

`profiles/slurm/config.yaml` (Snakemake ≥ 8 executor-plugin style; `--profile profiles/slurm`):

```yaml
executor: slurm
jobs: 200
latency-wait: 120
rerun-incomplete: true
keep-going: true
printshellcmds: true
software-deployment-method:
  - conda
  - apptainer
conda-frontend: mamba
apptainer-args: "--bind /scratch,/data"
default-resources:
  slurm_account: myproj
  slurm_partition: core
  mem_mb: 4000
  runtime: 60
  tmpdir: "/scratch/$USER/tmp"
set-resources:
  samtools_stats:
    mem_mb: 8000
set-threads:
  samtools_stats: 8
```

Resource names the SLURM plugin reads: `mem_mb` (or `mem_mb_per_cpu`), `runtime` (minutes), `slurm_account`, `slurm_partition`, `slurm_extra` (e.g. `"'--qos=short --gres=gpu:1'"`), `cpus_per_task` (defaults to `threads`), `nodes`, `tasks`, `constraint`. GPU rules: `slurm_extra` plus `apptainer-args: "--nv"`. Dry-run: `snakemake --profile profiles/slurm -n`. Workflow-local defaults go in `workflow/profiles/default/config.yaml` (picked up automatically).

## 4. `sbatch` array script

`scripts/align_array.sbatch` — one sample per task; reads `samples.txt` (`sample<TAB>r1.fq.gz<TAB>r2.fq.gz`).

```bash
#!/usr/bin/env bash
#SBATCH --job-name=align
#SBATCH --account=myproj
#SBATCH --partition=core
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --tmp=300G
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err
set -euo pipefail

SAMPLES=${1:?usage: sbatch --array=1-N align_array.sbatch samples.txt}
REF=/data/ref/GRCh38.fa
IMG=/shared/containers/bwa-mem2_2.3_samtools_1.25.sif
OUTDIR=/data/results/bam
THREADS=${SLURM_CPUS_PER_TASK}
# Per-job scratch: always a fresh directory so the EXIT trap can never remove a shared SLURM_TMPDIR.
SCRATCH=$(mktemp -d -p "${SLURM_TMPDIR:-/scratch/${USER}}" "align.${SLURM_JOB_ID}.XXXXXX")

mkdir -p "${OUTDIR}"
trap 'rm -rf "${SCRATCH}"' EXIT

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
SAMPLE=$(cut -f1 <<< "${LINE}")
R1=$(cut -f2 <<< "${LINE}")
R2=$(cut -f3 <<< "${LINE}")
# Samplesheet fields are data, never shell: restrict the sample ID and pass paths as arguments.
[[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ && -s "${R1}" && -s "${R2}" ]] || { echo "bad line ${SLURM_ARRAY_TASK_ID}: ${LINE}" >&2; exit 2; }

echo "$(date -Is) start ${SAMPLE} on ${HOSTNAME} cpus=${THREADS} mem=${SLURM_MEM_PER_NODE:-NA}"
apptainer exec --bind /data,"${SCRATCH}" "${IMG}" bash -c '
  set -euo pipefail
  threads=$1; sample=$2; ref=$3; r1=$4; r2=$5; scratch=$6
  bwa-mem2 mem -t "${threads}" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "${ref}" "${r1}" "${r2}" \
    | samtools sort -@ 4 -m 2G -T "${scratch}/${sample}.sort" -o "${scratch}/${sample}.bam" -
  samtools index "${scratch}/${sample}.bam"
  samtools quickcheck "${scratch}/${sample}.bam"
' _ "${THREADS}" "${SAMPLE}" "${REF}" "${R1}" "${R2}" "${SCRATCH}"
mv "${SCRATCH}/${SAMPLE}.bam" "${SCRATCH}/${SAMPLE}.bam.bai" "${OUTDIR}/"
apptainer exec "${IMG}" bash -c 'bwa-mem2 version; samtools --version | head -1' > "${OUTDIR}/${SAMPLE}.versions.txt"
echo "$(date -Is) done ${SAMPLE}"
```

Submit: `mkdir -p logs; sbatch --array=1-$(wc -l < samples.txt)%50 scripts/align_array.sbatch samples.txt`. Post-check: `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS -P | awk -F'|' '$2!="COMPLETED"'`. Inspect `MaxRSS` to shrink `--mem` next time. Array limits: `scontrol show config | grep MaxArraySize`.

## 5. Pinned `environment.yml` (conda/mamba) and lock files

```yaml
name: wgs-germline
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - python=3.12.*
  - samtools=1.25
  - bcftools=1.25
  - bwa-mem2=2.3
  - gatk4=4.6.2.0
  - multiqc=1.35
  - snakemake-minimal=9.27.2
  - snakemake-executor-plugin-slurm
  - pip
  - pip:
      - truvari==5.5.0
```

Check each pin exists before shipping: `mamba search -c bioconda samtools=1.25` (or `mamba repoquery search`). bioconda lags GitHub: gatk4 4.6.2.0 on bioconda vs 4.7.0.0 upstream on 2026-11-22; pin whichever channel you install from and record it. Other bioconda versions seen the same day: samtools 1.25, bcftools 1.25, bwa-mem2 2.3, multiqc 1.35, snakemake-minimal 9.27.2, truvari 5.5.0. Re-check before copying (`https://anaconda.org/bioconda/<pkg>`).

Lock and reproduce:

```bash
set -euo pipefail
ENV=environment.yml
conda-lock lock -f "${ENV}" -p linux-64 --lockfile conda-lock.yml
conda-lock install -n wgs-germline conda-lock.yml
conda env export -n wgs-germline --no-builds > environment.lock.export.yml
```

Pixi equivalent: `pixi init && pixi add samtools=1.25 bcftools=1.25` → commit `pixi.toml` + `pixi.lock`; run with `pixi run samtools …`.

Container equivalent of an env: `wave` with `conda.enabled` builds and freezes an image from the same spec (`-profile wave,conda`), or build once with `micromamba` in a Dockerfile and pin the resulting digest.

## 6. `versions.txt` capture snippet

Append to any ad-hoc bash script after `set -euo pipefail`; produces a sorted, de-duplicated `versions.txt` in the output directory plus the reference checksum.

```bash
set -euo pipefail
OUTDIR=results
REF=/data/ref/GRCh38.fa
mkdir -p "${OUTDIR}"
{
  echo "date	$(date -Is)"
  echo "host	$(hostname)"
  echo "user	${USER}"
  echo "git_commit	$(git rev-parse HEAD 2>/dev/null || echo NA)"
  echo "reference	${REF}	$(md5sum "${REF}" | cut -d' ' -f1)"
  for t in samtools bcftools bwa-mem2 gatk nextflow snakemake truvari rtg; do
    if command -v "${t}" >/dev/null 2>&1; then
      v=$("${t}" --version 2>&1 | head -n1 | tr -s ' ' | tr '\t' ' ')
      echo "${t}	${v}"
    else
      echo "${t}	MISSING"
    fi
  done
  if [[ -n "${SINGULARITY_CONTAINER:-}${APPTAINER_CONTAINER:-}" ]]; then
    echo "container	${SINGULARITY_CONTAINER:-${APPTAINER_CONTAINER}}"
  fi
  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    echo "conda_env	${CONDA_PREFIX}"
    conda list --export 2>/dev/null | grep -v '^#' | awk '{print "conda_pkg\t"$0}'
  fi
} | sort -u > "${OUTDIR}/versions.txt"
```

Validate: `column -t -s$'\t' "${OUTDIR}/versions.txt" | head -20`. Nextflow users: rely on `versions.yml` per module (template 1) gathered into `pipeline_info/software_versions.yml`; Snakemake users: `conda env export` per `conda:` env (`snakemake --conda-list-envs`, then export each prefix) or `--report` which embeds them.
