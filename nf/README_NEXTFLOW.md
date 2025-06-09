# Micro-C Pipeline - Nextflow DSL2 Implementation

A comprehensive Nextflow pipeline for processing Micro-C sequencing data, converted from the original Python implementation.

**Note**: This Nextflow implementation is located in the `nf/` subdirectory to maintain clear separation from the original Python pipeline in the project root.

## Overview

This pipeline processes Micro-C sequencing data through the following steps:

1. **BWA Alignment** - Align paired-end reads with BWA and process with pairtools
2. **Merge & Deduplication** - Merge multiple samples, remove duplicates, and split into pairs/BAM
3. **Cooler Processing** - Generate contact matrices in cooler format (raw and balanced mcool files)
4. **Juicer HiC** - Generate HiC files using Juicer tools (with error handling)
5. **QC Metrics** - Extract quality control statistics

## Quick Start

### Prerequisites

- Nextflow (≥21.10.3)
- Docker, Singularity, or Conda/Mamba
- Java 11+ (for Juicer tools)

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd microc-pipeline

# Test the pipeline
nextflow run nf/main.nf -profile test,docker
```

### Basic Usage

```bash
nextflow run nf/main.nf \
  --sample_id my_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  --output_dir results \
  -profile docker
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--sample_id` | Sample identifier |
| `--fastq_r1` | R1 FASTQ file(s) (comma-separated for multiple) |
| `--fastq_r2` | R2 FASTQ file(s) (comma-separated for multiple) |
| `--reference_bwa_idx` | BWA index tar.gz file |
| `--chrom_sizes` | Chromosome sizes file |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | `./results` | Output directory |
| `--resolution` | `1000` | Contact matrix resolution (bp) |
| `--bwa_cores` | `2` | Number of BWA alignment cores |
| `--mapq` | `20` | Minimum mapping quality |

## Profiles

### Execution Profiles

- `local` - Run locally (default)
- `cluster` - Run on SLURM cluster
- `cloud` - Run on cloud (Google Batch)

### Container Profiles

- `docker` - Use Docker containers
- `singularity` - Use Singularity containers
- `conda` - Use Conda environments
- `mamba` - Use Mamba environments

### Test Profiles

- `test` - Run with small test dataset
- `test_full` - Run with larger test dataset

## Usage Examples

### Single Sample

```bash
nextflow run nf/main.nf \
  --sample_id sample1 \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  -profile docker
```

### Multiple FASTQ Pairs

```bash
nextflow run nf/main.nf \
  --sample_id sample1 \
  --fastq_r1 "reads1_R1.fq.gz,reads2_R1.fq.gz" \
  --fastq_r2 "reads1_R2.fq.gz,reads2_R2.fq.gz" \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  -profile docker
```

### High-Performance Run

```bash
nextflow run nf/main.nf \
  --sample_id large_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  --bwa_cores 16 \
  --resolution 5000 \
  -profile cluster
```

### With Custom Resources

```bash
nextflow run nf/main.nf \
  --sample_id sample1 \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  --max_cpus 32 \
  --max_memory 256.GB \
  --max_time 48.h \
  -profile cluster
```

## Output Structure

```
results/
├── align/
│   └── chunk_*/
│       ├── genome_index/     # Extracted BWA index
│       └── *.pairsam.gz      # Aligned and parsed reads
├── merge/
│   ├── *.mapped.pairs        # Merged and deduplicated pairs
│   ├── *.bam                 # Sorted BAM file
│   ├── *.bam.bai             # BAM index
│   └── *.stats.txt           # Pairtools statistics
├── cooler/
│   ├── *.cool                # Contact matrix
│   ├── *.raw.mcool           # Multi-resolution raw matrix
│   └── *.balanced.mcool      # Multi-resolution balanced matrix
├── hic/
│   ├── *.hic                 # Juicer HiC file (may be empty if failed)
│   └── *.juicer.pairs        # Juicer-formatted pairs
├── qc/
│   ├── *_qc.txt              # QC metrics (text format)
│   └── *_qc_summary.json     # QC metrics (JSON format)
└── pipeline_info/
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.svg
```

## Resource Requirements

### Minimum Requirements

- **CPU**: 4 cores
- **Memory**: 8 GB RAM
- **Storage**: 50 GB free space

### Recommended for Production

- **CPU**: 16+ cores
- **Memory**: 64+ GB RAM
- **Storage**: 500+ GB free space (NVMe SSD preferred)

### Process-Specific Requirements

| Process | CPU | Memory | Time |
|---------|-----|--------|------|
| Extract BWA Index | 2 | 4 GB | 1h |
| BWA Align | 2-16 | 16 GB | 8h |
| Merge/Dedup | 12 | 32 GB | 12h |
| Cooler | 8 | 24 GB | 6h |
| Juicer HiC | 2-16 | 128 GB | 4h |
| QC Metrics | 2 | 4 GB | 1h |

## Configuration

### Custom Configuration

Create a custom configuration file:

```groovy
// custom.config
params {
    max_cpus = 32
    max_memory = '256.GB'
    max_time = '48.h'
}

process {
    withName: 'BWA_ALIGN_PAIRTOOLS' {
        cpus = 16
        memory = '32.GB'
    }
}
```

Use with: `nextflow run nf/main.nf -c custom.config`

### Cluster Configuration

For SLURM clusters:

```groovy
// cluster.config
process {
    executor = 'slurm'
    queue = 'normal'
    clusterOptions = '--account=your_account --partition=compute'
}
```

## Troubleshooting

### Common Issues

1. **Juicer HiC fails**: This is expected with some data formats. The pipeline continues and cooler files provide equivalent functionality.

2. **Out of memory**: Increase memory allocation or reduce the number of parallel processes.

3. **Slow performance**: Use faster storage (SSD) and increase CPU allocation.

### Error Handling

The pipeline includes robust error handling:

- **Automatic retries** for transient failures
- **Graceful degradation** for Juicer HiC failures
- **Resource scaling** on retry attempts
- **Detailed error logs** in work directories

### Debugging

Enable debug mode:

```bash
nextflow run nf/main.nf -profile debug,docker --sample_id test ...
```

Check work directories for detailed logs:

```bash
ls -la work/
```

## Performance Optimization

### For Large Datasets

1. **Increase parallelization**:
   ```bash
   --bwa_cores 16 --max_cpus 32
   ```

2. **Use faster storage**:
   - NVMe SSD for work directory
   - High-bandwidth network storage

3. **Optimize memory**:
   ```bash
   --max_memory 256.GB
   ```

4. **Use appropriate resolution**:
   ```bash
   --resolution 5000  # Lower resolution = faster processing
   ```

## Comparison with Python Pipeline

| Feature | Python Pipeline | Nextflow Pipeline |
|---------|----------------|-------------------|
| **Parallelization** | Limited | Full workflow parallelization |
| **Resource Management** | Manual | Automatic with retry logic |
| **Scalability** | Single machine | Local to cloud |
| **Reproducibility** | Environment dependent | Container/conda based |
| **Monitoring** | Basic logging | Rich reports and visualization |
| **Error Handling** | Basic | Robust with retries |
| **Resume Capability** | None | Full workflow resume |

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes and test
4. Submit a pull request

## License

This project is licensed under the MIT License.

## Citation

If you use this pipeline, please cite:

```
Micro-C Pipeline Nextflow Implementation
https://github.com/your-org/microc-pipeline-nextflow
```
