---
name: genomics-workflow
description: Runs and authors genomics pipelines reproducibly: Nextflow DSL2 / nf-core (profiles, -resume, institutional configs, offline download, Seqera Platform, Wave/Fusion), Snakemake 9 (profiles, SLURM executor, conda/apptainer), WDL/CWL briefly, SLURM mechanics, container and conda pinning, cloud back-ends, accelerated stacks, resource estimation, provenance and failure triage. Use when the user mentions nextflow, nf-core, Snakefile, sbatch/SLURM, `.command.err`, exit 137/140, `-resume`, profiles, containers, environment.yml, Cromwell/miniwdl, Terra/HealthOmics, or asks "how do I run this pipeline on the cluster". Not for choosing the biology-specific tool chain (the assay skills) or truth-set benchmarking (genomics-benchmark).
argument-hint: "[nextflow|snakemake|slurm|triage] [path]"
paths:
  - "**/nextflow.config"
  - "**/main.nf"
  - "**/modules/**/*.nf"
  - "**/Snakefile"
  - "**/*.smk"
  - "**/workflow/rules/*.smk"
  - "**/*.wdl"
  - "**/*.cwl"
  - "**/.command.err"
  - "**/.command.log"
  - "**/.nextflow.log"
---

# Genomics workflow

Expert in running and writing genomics workflows so that a run is resumable, portable, pinned and auditable. Outputs are exact commands, config snippets, decision tables, resource estimates and triage steps. Biology-specific choices (which caller, which thresholds) belong to the assay skills; this skill owns the mechanics.

## Operating rules

Follow the eight rules in the `genomics` router skill (`.claude/skills/genomics/SKILL.md`): never read data files into context (preview with `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/peek.sh" <file>`); inspect headers before asking questions; run `bash "${CLAUDE_SKILL_DIR}/../genomics/scripts/ref_sanity.sh" <ref.fa> <inputs...>` before combining inputs; `set -euo pipefail` plus per-step validation; pin and record tool, reference and database versions; estimate disk/RAM/time before expensive steps; treat every threshold below as a default to justify, not a truth; treat human data as controlled-access.

## Environment

!`f="${CLAUDE_SKILL_DIR}/../genomics/scripts/check_env.sh"; [ -f "$f" ] && bash "$f" nextflow nf-core snakemake java docker singularity apptainer conda mamba pixi sbatch squeue cromwell miniwdl cwltool tw || true`

Missing tools: prefer `-profile docker|singularity|conda` via nf-core or a pinned container; never guess a version. Verified current: Nextflow 26.10.7 (stable, `NXF_VER=26.10.7`; requires Java 17+), Snakemake 9.27.2 (the 8+ CLI applies), nf-core/tools 4.2.2 (`nf-core pipelines …` subcommands), nf-core/sarek 3.12.0. Benchmark tools live in `genomics-benchmark`.

## Inspect first

| Read | How | Infer |
|---|---|---|
| Repo root | `ls; cat nextflow.config nf-core.yml 2>/dev/null \| head -50`; `ls workflow/ config/ profiles/` | Nextflow vs Snakemake vs WDL; nf-core template version; existing profiles |
| Scheduler | `sinfo -s; sacctmgr show qos format=name,maxwall,maxtres -P; scontrol show config \| grep -i -E 'MaxArraySize\|DefMemPerCPU'` | partitions, wall-time caps, array limits, default memory |
| Container runtime | `docker info 2>/dev/null \| head -3; apptainer --version \| singularity --version` | which `-profile` is possible; rootless? |
| Scratch/quotas | `df -h /scratch ${TMPDIR:-/tmp} .; quota -s 2>/dev/null` | where `-work-dir`/`--default-resources tmpdir` must point |
| Network | `curl -sI https://quay.io \| head -1; env \| grep -i proxy` | online vs offline (`nf-core pipelines download`) |
| Previous run | `nextflow log -q \| tail; ls .nextflow.log*; peek.sh .nextflow.log 60`; Snakemake: `ls .snakemake/log \| tail -3` | failed task hash, resume eligibility |
| Inputs | `peek.sh samplesheet.csv`; `ref_sanity.sh ref.fa <bam/vcf/bed>` | columns, build consistency before launching |

Ask only if still ambiguous: SLURM account/partition and QOS to charge; whether internet is reachable from compute nodes; institutional nf-core config name (`https://github.com/nf-core/configs`); storage for the work dir (needs 3–5× final output); GPU availability for accelerated stacks.

## Decision table

| Situation | Choice | Why |
|---|---|---|
| An nf-core pipeline exists for the assay | `nextflow run nf-core/<pipe> -r <tag> -profile <container>,<institution>` | tested, containerised, provenance (MultiQC, versions) built in; never fork for a parameter you can set in a config |
| nf-core exists but one step must change | `-c custom.config` with `process { withName: 'X' { ext.args = … } }`; else nf-core `modules` swap in a fork | keeps upgrade path |
| No nf-core pipeline; team writes Groovy/Java-tolerant, cloud or Seqera target | Nextflow DSL2 with nf-core modules (`nf-core modules install <tool>`) | module registry, Wave/Fusion, Platform |
| No nf-core; Python team, single HPC, exploratory | Snakemake 9 workflow with `conda:`/`container:` per rule and a SLURM profile | lower ceremony; wrappers catalog |
| Broad/Terra ecosystem or GATK WDL already exists | WDL on Cromwell (Terra, server) or `miniwdl run` locally | reuse published WDL; WDL↔Nextflow ports are a source of bugs |
| Regulatory/CWL requirement (ELIXIR, Arvados, Seven Bridges) | CWL via `cwltool` or `toil-cwl-runner` | portability mandated; otherwise avoid |
| One-off per-sample command over N samples on SLURM | `sbatch --array` script (templates.md) | no framework overhead; still pin container |
| Nodes have no internet | `nf-core pipelines download … --container-system singularity`, `NXF_OFFLINE=true` | pre-fetch code + images |
| GPU nodes, licence available | Parabricks (`fq2bam`, `deepvariant`, `haplotypecaller`), DRAGEN, Sentieon | drop-in replacements for BWA/GATK/DeepVariant; benchmark once with `genomics-benchmark` before trusting |
| Cloud | see table below | managed executors; keep code identical, swap `-profile` |

Cloud back-ends (one row each; executor is set in config, not code):

| Platform | Nextflow setting | Notes |
|---|---|---|
| AWS Batch / HealthOmics | `process.executor = 'awsbatch'`, `aws.batch.cliPath`, S3 work dir; HealthOmics runs Nextflow/WDL/CWL natively (`aws omics start-run`) | HealthOmics needs the pipeline packaged as a zip + ECR images |
| Google Batch | `process.executor = 'google-batch'`, `google.project`, `google.location`, `gs://` work dir | preemptible via `google.batch.spot = true` |
| Azure Batch | `process.executor = 'azurebatch'`, `azure.storage.accountName`, `az://` work dir | pool autoscale in `azure.batch.pools` |
| Terra / AnVIL | WDL on Cromwell via Terra UI/API; Nextflow via Seqera Platform compute env | controlled-access data stays in the workspace; moving data in or out needs the DUA/consent check (rule 8) |
| Seqera Platform (any of the above + HPC) | `nextflow launch` or `tw launch <pipeline> --compute-env <ce> --params-file params.yaml`; `-with-tower` for monitoring | Wave builds/augments containers; Fusion mounts object storage. `-with-tower` ships task names, `tag` values (often sample IDs), paths and logs to Seqera's service: confirm the DUA allows it or use an on-prem Platform/Enterprise instance |

## Workflow

1. Nextflow / nf-core run (canonical form).

```bash
set -euo pipefail
export NXF_VER=26.10.7
export NXF_SINGULARITY_CACHEDIR=/shared/containers/nxf
export NXF_WORK=/scratch/${USER}/nxf_work
PIPE=nf-core/sarek
TAG=3.12.0
INST=uppmax
nextflow run "${PIPE}" -r "${TAG}" -profile singularity,"${INST}" \
  -params-file params.yaml -work-dir "${NXF_WORK}" -resume \
  -with-report report.html -with-trace trace.txt -with-timeline timeline.html -with-dag dag.html
```

Validate: `nextflow log -q | tail -1` then `nextflow log <run_name> -f name,status,exit,realtime,peak_rss -F 'status != "COMPLETED"'` (empty = success). Check the pipeline tag exists first: `nextflow pull "${PIPE}"; nextflow info "${PIPE}" | grep revisions -A20`. Institutional configs: `https://github.com/nf-core/configs/tree/master/conf` (`-profile <name>` works for any nf-core pipeline). Build `params.yaml` interactively with `nf-core pipelines launch "${PIPE}" -r "${TAG}"`. Re-run with `-resume` from the same directory; `-resume <session_id>` to pick an older session. Replace `singularity` by `docker`, `apptainer`, `conda`, `podman`, `charliecloud` as available; Wave + Fusion for object-store work dirs: `-with-wave` plus `fusion.enabled = true` in config (needs Wave).

2. Offline / air-gapped nf-core.

```bash
set -euo pipefail
PIPE=sarek
TAG=3.12.0
OUT=/shared/pipelines
nf-core pipelines download "${PIPE}" -r "${TAG}" --outdir "${OUT}/${PIPE}-${TAG}" \
  --container-system singularity --container-cache-utilisation amend --compress none
export NXF_OFFLINE=true
export NXF_SINGULARITY_CACHEDIR=/shared/containers/nxf
nextflow run "${OUT}/${PIPE}-${TAG}/${TAG//./_}" -profile singularity -params-file params.yaml -resume
```

Validate: `ls "${NXF_SINGULARITY_CACHEDIR}" | wc -l` ≥ number of distinct containers (`nextflow inspect <pipe> -profile singularity | jq '.processes[].container' | sort -u | wc -l`).

3. Snakemake 9 run with SLURM executor and pinned software.

```bash
set -euo pipefail
PROFILE=profiles/slurm
snakemake --profile "${PROFILE}" --software-deployment-method conda apptainer \
  --rerun-incomplete --keep-going --printshellcmds --jobs 200 --rerun-triggers mtime params input software-env
```

Validate: `snakemake --profile "${PROFILE}" -n --quiet rules` (dry run, nothing missing) then `snakemake --summary | awk '$5=="missing"' | wc -l` = 0 after the run. Install the executor once: `pip install snakemake-executor-plugin-slurm`. `--sdm` is the short form. Report: `snakemake --report report.html`. Profile layout is in `templates.md`.

4. WDL / CWL (brief).

```bash
set -euo pipefail
WDL=workflows/germline.wdl
INPUTS=inputs.json
miniwdl check "${WDL}"
miniwdl run "${WDL}" -i "${INPUTS}" --dir runs/ --cfg miniwdl.cfg
```

Validate: `jq '.outputs' runs/*/outputs.json`. Server mode: `java -jar cromwell.jar run "${WDL}" -i "${INPUTS}" -o options.json` (`options.json` sets `final_workflow_outputs_dir`, `use_relative_output_paths`). CWL: `cwltool --validate wf.cwl` then `cwltool --singularity --outdir out wf.cwl job.yml`.

5. SLURM direct (array over a sample list).

```bash
set -euo pipefail
SAMPLES=samples.txt
N=$(wc -l < "${SAMPLES}")
JOB=$(sbatch --parsable --array=1-"${N}"%50 scripts/align_array.sbatch "${SAMPLES}")
echo "submitted ${JOB}"
```

Validate: `sacct -j "${JOB}" --format=JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem -P | awk -F'|' '$2!="COMPLETED"'` (empty = all good; `MaxRSS` vs `ReqMem` tells you how to right-size). Script in `templates.md`: `--cpus-per-task`, `--mem`, `--time`, `--tmp`, `$SLURM_TMPDIR`/`/scratch` staging, `trap` cleanup, container exec.

6. Pin software and record provenance in every run directory.

```bash
set -euo pipefail
OUTDIR=results
mkdir -p "${OUTDIR}/provenance"
nextflow -version > "${OUTDIR}/provenance/nextflow_version.txt"
cp "${OUTDIR}/pipeline_info/"*software_versions.yml "${OUTDIR}/provenance/" 2>/dev/null || true
git -C . rev-parse HEAD > "${OUTDIR}/provenance/repo_commit.txt" 2>/dev/null || true
conda env export --no-builds > "${OUTDIR}/provenance/environment.lock.yml" 2>/dev/null || true
```

Validate: `grep -c ':' "${OUTDIR}/provenance/"*software_versions.yml` > 0; nf-core pipelines also emit `pipeline_info/execution_report_*.html`, `execution_trace_*.txt`, `params_*.json`. For containers pin by digest (`image@sha256:…`, resolve with `docker buildx imagetools inspect <image:tag>` or `skopeo inspect docker://<image:tag> | jq -r .Digest`); for conda use `conda-lock -f environment.yml -p linux-64` and ship `conda-lock.yml`; for pixi commit `pixi.lock`. The `versions.txt` capture snippet for ad-hoc scripts is in `templates.md`.

7. Clean up — only after the user confirms no `-resume` is needed for the session being cleaned: `nextflow clean -f -before <run_name>` (run name from `nextflow log -q`; `-k` also drops the run's log entry).

Validate: `du -sh "${NXF_WORK:-work}"`. Snakemake: `snakemake --delete-temp-output`, `rm -rf .snakemake/conda` only after `conda env export`. Never delete the latest successful session's work dir if you still want `-resume`.

## QC thresholds

Approximate resource defaults for 30× human WGS unless stated; measure with `trace.txt`/`sacct` and right-size on the second run.

| Step (tool) | CPUs | RAM | Wall | Scratch | Rationale / when to deviate |
|---|---|---|---|---|---|
| fastp / FastQC per sample | 4 | 4 GB | 0.5 h | 1× FASTQ | I/O bound; more threads rarely help |
| bwa-mem2 + sort (30× WGS) | 16 | 64 GB | 5 h | 3× FASTQ | bwa-mem2 index needs ~30 GB RAM alone; bwa 0.7.x fits 32 GB |
| minimap2 ONT/HiFi 30× | 16 | 32 GB | 6 h | 3× FASTQ | `-x map-ont`/`map-hifi`; sort separately |
| GATK MarkDuplicates / `samtools markdup` | 8 | 32 GB | 3 h | 2× BAM | MarkDuplicatesSpark halves wall at 16 CPUs |
| HaplotypeCaller per scatter shard (50 shards) | 2 | 8 GB | 1 h | small | scatter by `--intervals`; gather with GatherVcfs |
| DeepVariant 1.10.0 CPU | 32 | 64 GB | 4 h | 1× BAM | GPU ≈ 1 h; `make_examples` scales linearly |
| Mutect2 per shard | 2 | 8 GB | 2 h | small | PoN/gnomAD lookups add I/O |
| STAR (2-pass, per RNA sample) | 8 | 40 GB | 0.5 h | 2× FASTQ | index ≈ 32 GB RAM for human; `--genomeLoad` shares it |
| salmon / kallisto quant | 8 | 16 GB | 0.3 h | small | — |
| MACS3 callpeak | 1 | 8 GB | 0.5 h | small | single-threaded |
| VEP 116 (cache, `--fork 8`) | 8 | 16 GB | 1 h / 5 M variants | cache 20+ GB | plugins (CADD, SpliceAI) add files, not RAM |
| MultiQC | 1 | 4 GB | 0.2 h | small | thousands of samples: 16 GB |
| Work dir total, nf-core/sarek WGS | — | — | — | 4–5× final outputs | plan disk before launching; `-resume` needs it kept |

Scheduler rules: request RAM at 1.3× observed `peak_rss`; request time at 1.5× observed; nf-core `-profile` retries (`errorStrategy retry`, `maxRetries 2`) escalate automatically for exit 137/140/143 — do not disable them. Hundreds of tiny tasks → `executor.queueSize = 100`, `executor.submitRateLimit = '10 sec'`.

## Pitfalls

| Symptom | Cause | Fix |
|---|---|---|
| Exit 137, `.command.log` ends abruptly, `sacct` shows `OUT_OF_MEMORY` | cgroup OOM kill | raise `memory` (×2) or lower threads; check `peak_rss` in trace; Java tools: `-Xmx` must be < requested mem |
| Exit 140 (SLURM) / `TIMEOUT` | wall time exceeded | raise `time`; split intervals; on nf-core use `--max_time` or a `withName` override |
| Exit 143 | SIGTERM from scheduler (preemption, node drain, QOS) | retry; spot/preemptible → `errorStrategy 'retry'` |
| Exit 139 | segfault, often wrong CPU arch for container or corrupt index | rebuild index on that node type; pin image digest |
| Exit 127 / `command not found` inside task | container lacks tool or `-profile conda` not activated | `cat .command.run` to see the container used; `nextflow inspect` |
| `-resume` reruns everything | work dir moved/deleted, inputs touched (mtime), different `-params`, process with non-deterministic input order | keep work dir; use `-work-dir` explicitly; `cache 'lenient'` for shared FS with flaky mtimes |
| `Process requirement exceeds available memory` | local executor, `memory` > host | use the cluster profile or `--max_memory` |
| Singularity `FATAL: while pulling … unauthorized` or disk full | no cache dir, pulls into `work/singularity` per run | set `NXF_SINGULARITY_CACHEDIR`; pre-pull with `nf-core pipelines download` |
| `Failed to pull singularity image … exit 255 / timeout` | registry rate limit, no internet on compute node | pull from login node once; `singularity.pullTimeout = '2h'` |
| Snakemake `IncompleteFilesException` | previous run killed mid-write | `--rerun-incomplete` |
| Snakemake `MissingOutputException` after job finished | NFS latency | `--latency-wait 120` |
| Snakemake SLURM jobs stay `PENDING` | missing `slurm_account`/`slurm_partition` resources | set them in the profile `default-resources` |
| Conda env solve takes hours / fails | unpinned channels, libmamba disabled | pin versions, `channel_priority: strict`, `conda-lock` |
| Same pipeline, different results on two clusters | floating tags (`:latest`), different reference, thread-dependent tools | pin by digest; record `software_versions.yml`; set `--seed` where tools support it |
| Pipeline fills `$HOME` | `NXF_HOME`, `~/.singularity`, conda pkgs in home | export `NXF_HOME`, `APPTAINER_CACHEDIR`, `CONDA_PKGS_DIRS` to scratch |
| nf-core run ignores my param | typo not in schema, or set in config `params {}` after CLI | `nf-core pipelines schema validate`; CLI > `-params-file` > config |
| Cloud run: tasks fail copying inputs | missing object-store permissions / `aws.region` | test `nextflow fs ls s3://bucket/`; use Fusion to avoid staging |

Failure triage, in order: `nextflow log <run> -f name,status,exit,workdir -F 'status == "FAILED"'` → `cd <workdir>` → `cat .command.sh` (exact command) → `tail -50 .command.err` (tool stderr) → `.command.log` (combined) → `.exitcode` → `cat .command.run | grep -E 'singularity|docker|conda'` (environment) → re-run `bash .command.run` interactively on a compute node with the same resources. Snakemake: `.snakemake/log/*.log`, then the rule's `log:` file, then `--debug-dag`.

## nf-core

This skill runs any nf-core pipeline; the assay skill picks which one and owns its samplesheet header (verify existence at `https://nf-co.re/pipelines.json`; list releases with `nextflow info nf-core/<pipe>`). Pipelines verified to exist today and their owning skill:

| Pipeline | Owning skill (samplesheet header there) |
|---|---|
| nf-core/sarek (3.12.0) | `.claude/skills/genomics-wgs/SKILL.md`, `genomics-somatic`, `genomics-cnv` |
| nf-core/rnaseq, rnafusion, smrnaseq, scrnaseq, spatialvi | `.claude/skills/genomics-rnaseq/SKILL.md` |
| nf-core/chipseq, atacseq, cutandrun | `.claude/skills/genomics-chipseq/SKILL.md` |
| nf-core/methylseq | `.claude/skills/genomics-methylation/SKILL.md` |
| nf-core/nanoseq, longraredisease | `.claude/skills/genomics-longread/SKILL.md` |
| nf-core/raredisease | `.claude/skills/genomics-annotation/SKILL.md` |
| nf-core/oncoanalyser | `.claude/skills/genomics-somatic/SKILL.md` |
| nf-core/variantbenchmarking (1.6.1) | `.claude/skills/genomics-benchmark/SKILL.md` (header `id,test_vcf,caller`) |

Generic launch for any of them:

```bash
set -euo pipefail
PIPE=nf-core/rnaseq
TAG=$(nextflow info "${PIPE}" | sed -n '/revisions/,$p' | tr -d ' *' | grep -E '^[0-9]+\.[0-9]+' | cut -d'[' -f1 | sort -V | tail -1)
nextflow run "${PIPE}" -r "${TAG}" -profile test,docker --outdir test_out
```

`nf-core/copynumber` and `nf-core/fastqc` do not exist. Always `-r <tag>`; never `-r dev` for production. Test a new setup with `-profile test,<container>` first.

## Validation

- Run completed: `nextflow log <run> -f status | sort | uniq -c` shows only `COMPLETED`/`CACHED`; Snakemake exits 0 and `--summary` has no `missing`.
- Reproducibility proof: re-launch with `-resume` → every task `CACHED`; re-launch on a second machine with the same digests → identical `md5sum` of deterministic outputs (BAM headers may differ in `@PG`; compare `samtools view` bodies).
- Provenance present: `pipeline_info/*software_versions.yml`, `params_*.json`, `execution_trace_*.txt`, MultiQC report, repo commit, container digests, reference checksums.
- Resource fit: no task retried for 137/140 in the final run; `peak_rss`/`requested` between 0.5 and 0.8.
- Scientific correctness is validated by the assay skill and, for variant calls, by `.claude/skills/genomics-benchmark/SKILL.md` (GIAB/SEQC2 F1 tables): a green pipeline run is not a correct result.

## Additional resources

- `.claude/skills/genomics-workflow/templates.md` — copy-ready Nextflow DSL2 module + workflow, `nextflow.config` profiles (docker/singularity/slurm), Snakemake `Snakefile` + SLURM profile `config.yaml`, `sbatch` array script, pinned `environment.yml`, `versions.txt` capture; load when writing or scaffolding a pipeline.
- `.claude/skills/genomics/SKILL.md` — router and the eight operating rules; `.claude/skills/genomics/reference.md` for reference/annotation acquisition (iGenomes vs refgenie vs custom `--fasta`).
- `.claude/skills/genomics-benchmark/SKILL.md` — proving the pipeline's calls are right after it runs green.
- `.claude/skills/genomics-qc/SKILL.md` — read-level QC/trimming decisions that feed `--skip_trimming`-style params.
- Upstream references: Nextflow docs `https://docs.seqera.io/nextflow`, nf-core configs `https://github.com/nf-core/configs`, nf-core tools CLI `https://nf-co.re/docs/nf-core-tools`, Snakemake CLI `https://snakemake.readthedocs.io/en/stable/executing/cli.html`, SLURM executor plugin `https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html`, Parabricks 4.8 `https://docs.nvidia.com/clara/parabricks/`.
