# preprocess_hic.wdl without Docker (SLURM) and with Arima libraries

`preprocess_hic.wdl` was written for Docker images on Google Cloud. This branch lets it run on a SLURM cluster without Docker, with HiC-Pro installed natively. It also documents how to run Arima (DpnII + HinfI) libraries. Leaving the new inputs unset gives the original Docker behaviour.

## Native toolchain

| Component | Source |
|---|---|
| HiC-Pro 3.1.0 | built from `nservant/HiC-Pro` tag `v3.1.0` (`make configure && make install`, `CLUSTER_SYS = SLURM`) |
| HiC-Pro dependencies | `envs/hicpro-3.1.0.yml`: upstream `environment.yml` pins, without the `defaults` channel; includes bowtie2 2.4.4, samtools 1.12, cooler 0.8.11 |
| Helpers | `envs/hic-utils.yml`: Java 17, pigz, pairix, htslib |
| Cromwell, juicer_tools | release jars |

The exact solved envs are committed as `envs/*.linux-64.explicit.txt`. The system also needs `zip`, `unzip` and `bc`.

## New workflow inputs

| Input | Default (unset) | Native value |
|---|---|---|
| `preprocess_hic.hicpro_dir` | `/HiC-Pro` | `<prefix>/HiC-Pro_3.1.0` |
| `preprocess_hic.annotation_dir` | `/HiC-Pro/annotation` (cooler: `/annotation`) | not needed if `genome_size` is an absolute path |
| `preprocess_hic.juicer_jar` | `/usr/local/juicer/juicer_tools.1.7.6_jcuda.0.8.jar` | path to a juicer_tools jar |
| `preprocess_hic.run_sparsehic` | `true` | `false` (the sparseHiC R script exists only in the Docker image) |
| `preprocess_hic.hicpro_align.genome_index_dir` | unset (untar `genome_index_tgz` per chunk) | directory with `<genome_name>.*.bt2`; `genome_index_tgz` can then be omitted |

`preprocess_hic.hicpro_align.min_mapq` was already an input (default `"20"`).

## Other changes

- Every HiC-Pro task edits a task-local copy of `config-hicpro.txt` instead of the copy in the install. In Docker each task had its own image copy. With a shared native install, parallel chunks would otherwise overwrite each other's settings.
- `split_fastq_files` uses `pigz` when it is on PATH and splits R1 and R2 concurrently. The chunk contents are unchanged.
- `cooler` chooses the smallest bin with a numeric sort. Before, it used a lexicographic sort, which picks 10000 over 5000.

## Arima (DpnII + HinfI)

- `hicpro_align.genome_fragment` is an absolute path to a fragment BED from HiC-Pro's `bin/utils/digest_genome.py -r ^GATC G^ANTC -o <out>.bed <genome>.fa`.
- `hicpro_align.ligation_site` is `"GATCGATC,GANTGATC,GANTANTC,GATCANTC"`.
- `genome_size` is an absolute path to a chrom.sizes file.

## Running on SLURM

```
export PATH=<envs>/hicpro-3.1.0/bin:<envs>/hic-utils/bin:$PATH
java -Dconfig.file=cromwell/slurm-native.conf -jar cromwell.jar run preprocess_hic.wdl -i inputs.json -o options.json
```

Run Cromwell itself as a long SLURM job; it submits each task with `sbatch`. Task `cpu` and `memory` runtime attributes become `-c` and `--mem`. `disks`, `preemptible` and `docker` are ignored.

## Notes from the first native run (2026-10-24)

- **HiC-Pro 3.x output names.** Valid pairs are written as `<sample>.allValidPairs` (2.x used `<sample>_allValidPairs`), and the merge stats go to `hic_results/stats/<sample>/`. The WDL accepts both layouts. The original Docker image probably shipped HiC-Pro 2.x.
- **pandas pin.** HiC-Pro's upstream `environment.yml` leaves pandas unpinned, and a fresh solve picks pandas 2.x. cooler 0.8.11 then fails to import (`is_categorical`). `envs/hicpro-3.1.0.yml` pins `pandas=1.3.5`.
- **cooler 0.8.x.** It has no `zoomify --no-balance` flag; balancing is off by default there. The WDL passes the flag only if the installed cooler supports it.
- **Where to compile.** HiC-Pro's C++ helpers (`scripts/cutsite_trimming`, `scripts/build_matrix`) must be compiled on the OS the jobs run on. Built on a newer login or interactive host, they need a newer libstdc++ (`GLIBCXX_3.4.20+ not found`). Recompile with `g++ -O2 -std=c++0x` on a compute node.
- **Cromwell config.**
  - Do not declare `docker` in the backend's `runtime-attributes`. If it is declared, Cromwell writes container paths (`/cromwell-executions/...`) into the job scripts.
  - Give `memory_mb` a default (8 GB here). An optional attribute with no default makes Cromwell reject tasks that set no memory.
- **Test results.** With mm10 Arima settings and MAPQ 30, a 2 M read-pair test on mouse lung tumour Arima Hi-C (2×51 bp) gave 56% deduplicated valid pairs, of which 20% were trans. The run took 14 min on SLURM with 1 M-pair chunks.
