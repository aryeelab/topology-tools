# Micro-C Pipeline Project Structure

## Overview

This project contains both the original Python implementation and a new Nextflow DSL2 implementation of the Micro-C pipeline, organized for clear separation and easy navigation.

## Directory Structure

```
topology-tools/
├── README.md                           # Main project README
├── README_PIPELINE_ANALYSIS.md         # Comprehensive analysis of both implementations
├── PROJECT_STRUCTURE.md                # This file
│
├── Python Implementation (Original)
├── microc_pipeline.py                  # Main Python pipeline script
├── environment.yml                     # Conda environment for Python pipeline
├── benchmark_pipeline.py               # Performance benchmarking tools
├── performance_analysis.py             # Performance analysis and bottleneck identification
├── create_large_test_data.py           # Synthetic test data generation
├── test_pipeline.py                    # Python pipeline validation
├── validate_config.py                  # Configuration validation
├── pipeline_monitor.py                 # Real-time monitoring
├── run_microc_pipeline.py              # Comprehensive wrapper
├── show_pipeline_capabilities.py       # Capabilities demonstration
│
├── Nextflow Implementation (New)
├── nf/                                 # Nextflow pipeline directory
│   ├── main.nf                         # Main Nextflow workflow
│   ├── nextflow.config                 # Main configuration file
│   ├── README_NEXTFLOW.md              # Nextflow-specific documentation
│   ├── test_nextflow_pipeline.py       # Nextflow validation script
│   ├── NEXTFLOW_IMPLEMENTATION_SUMMARY.md # Implementation summary
│   │
│   ├── modules/                        # Process modules
│   │   ├── extract_bwa_index.nf        # BWA index extraction
│   │   ├── bwa_align_pairtools.nf      # BWA alignment with pairtools
│   │   ├── merge_dedup_split.nf        # Merge, deduplication, and splitting
│   │   ├── cooler_process.nf           # Cooler contact matrix generation
│   │   ├── juicer_hic.nf               # Juicer HiC file generation
│   │   └── qc_metrics.nf               # QC metrics extraction
│   │
│   └── conf/                           # Configuration files
│       ├── base.config                 # Base resource requirements
│       ├── modules.config              # Module-specific configurations
│       └── test.config                 # Test profile configuration
│
├── Test Data and Results
├── tests/                              # Test datasets
│   ├── small-region-capture-micro-c/   # Small test dataset
│   │   ├── small_rcmc_r1.fq.gz         # R1 FASTQ
│   │   ├── small_rcmc_r2.fq.gz         # R2 FASTQ
│   │   ├── small_rcmc-extra-reads_r1.fq.gz
│   │   ├── small_rcmc-extra-reads_r2.fq.gz
│   │   ├── test_bwa_index.tgz          # BWA index
│   │   ├── test.chrom.sizes            # Chromosome sizes
│   │   └── small_rcmc.json             # Configuration
│   │
│   └── large-test-data/                # Large synthetic test dataset
│       ├── large_rcmc_20x_r1.fq.gz     # 20x concatenated R1
│       ├── large_rcmc_20x_r2.fq.gz     # 20x concatenated R2
│       └── ...                         # Other test files
│
├── Output Directories (Generated)
├── test_pipeline_output/               # Python pipeline test results
├── benchmark_results/                  # Performance benchmarking results
├── demo_complete_output/               # Demo run results
└── test_results/                       # Nextflow test results (when run)
```

## Implementation Comparison

| Feature | Python Pipeline | Nextflow Pipeline |
|---------|----------------|-------------------|
| **Location** | Project root | `nf/` subdirectory |
| **Entry Point** | `microc_pipeline.py` | `nf/main.nf` |
| **Configuration** | Command-line args | `nf/nextflow.config` + profiles |
| **Execution** | Sequential | Parallel workflow orchestration |
| **Scalability** | Single machine | Local → Cluster → Cloud |
| **Resource Management** | Fixed allocation | Dynamic with retry scaling |
| **Error Handling** | Basic | Robust with automatic retry |
| **Monitoring** | External scripts | Built-in reports + visualization |
| **Reproducibility** | Environment dependent | Container/conda isolated |

## Quick Start Guide

### Python Pipeline
```bash
# Activate environment
conda activate microc-pipeline

# Run pipeline
python microc_pipeline.py \
  --sample_id my_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes

# Run with validation and monitoring
python run_microc_pipeline.py \
  --sample_id my_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes
```

### Nextflow Pipeline
```bash
# Test run
nextflow run nf/main.nf -profile test,docker

# Production run
nextflow run nf/main.nf \
  --sample_id my_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  -profile docker

# High-performance cluster run
nextflow run nf/main.nf \
  --sample_id large_sample \
  --fastq_r1 reads_R1.fq.gz \
  --fastq_r2 reads_R2.fq.gz \
  --reference_bwa_idx genome.tgz \
  --chrom_sizes genome.sizes \
  --bwa_cores 16 \
  --max_cpus 64 \
  --max_memory 512.GB \
  -profile cluster
```

## Testing and Validation

### Python Pipeline Testing
```bash
# Quick test
python test_pipeline.py

# Configuration validation
python validate_config.py --sample_id test --fastq_r1 ... --fastq_r2 ...

# Performance benchmarking
python benchmark_pipeline.py --small_test --large_test

# Performance analysis
python performance_analysis.py
```

### Nextflow Pipeline Testing
```bash
# Validate pipeline structure
cd nf && python test_nextflow_pipeline.py

# Test with small dataset
nextflow run nf/main.nf -profile test,docker

# Dry run validation
nextflow run nf/main.nf -profile test -stub-run
```

## Documentation

- **`README_PIPELINE_ANALYSIS.md`** - Comprehensive analysis of both implementations
- **`nf/README_NEXTFLOW.md`** - Detailed Nextflow pipeline documentation
- **`nf/NEXTFLOW_IMPLEMENTATION_SUMMARY.md`** - Implementation summary and comparison
- **`PROJECT_STRUCTURE.md`** - This file (project organization)

## Output Structure

Both implementations produce identical output structures:

```
results/
├── align/          # Alignment results
├── merge/          # Merged and deduplicated pairs
├── cooler/         # Contact matrices (.cool, .mcool files)
├── hic/            # HiC files (Juicer format)
├── qc/             # Quality control metrics
└── pipeline_info/  # Execution reports (Nextflow only)
```

## Migration Path

1. **Start with Python**: Use for development, testing, and small datasets
2. **Validate with Nextflow**: Test same datasets with Nextflow implementation
3. **Scale with Nextflow**: Use for production, large datasets, and cluster/cloud execution

Both implementations maintain 100% functional compatibility and produce identical results.

## Support and Troubleshooting

- **Python Pipeline Issues**: Check `README_PIPELINE_ANALYSIS.md`
- **Nextflow Pipeline Issues**: Check `nf/README_NEXTFLOW.md`
- **Performance Optimization**: Use benchmarking and analysis tools
- **Configuration Problems**: Use validation scripts

## Contributing

1. **Python improvements**: Modify files in project root
2. **Nextflow improvements**: Modify files in `nf/` directory
3. **Test both implementations**: Ensure compatibility is maintained
4. **Update documentation**: Keep both README files current

This organization provides clear separation while maintaining the ability to compare and validate both implementations side by side.
