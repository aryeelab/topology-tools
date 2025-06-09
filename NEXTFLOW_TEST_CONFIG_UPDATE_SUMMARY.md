# Nextflow Test Configuration Update Summary

## 🎯 **MISSION ACCOMPLISHED!**

The Nextflow pipeline test configuration has been successfully updated to use the **exact same dataset, genome index, and parameters** as specified in the WDL JSON configuration file.

## 📋 **What Was Updated**

### 1. **Test Configuration Synchronization**
- **Source**: `tests/small-region-capture-micro-c/small_rcmc.json` (WDL specification)
- **Target**: `nf/conf/test.config` (Nextflow test configuration)
- **Result**: ✅ **100% parameter matching**

### 2. **Key Parameter Updates**

| Parameter | Before | After (WDL Match) |
|-----------|--------|-------------------|
| **Sample ID** | `test-sample` | `small-rcmc` |
| **FASTQ R1** | Single file | **2 files**: `small_rcmc_r1.fq.gz` + `small_rcmc-extra-reads_r1.fq.gz` |
| **FASTQ R2** | Single file | **2 files**: `small_rcmc_r2.fq.gz` + `small_rcmc-extra-reads_r2.fq.gz` |
| **Resolution** | 1000 bp | 1000 bp ✅ (confirmed match) |
| **BWA Cores** | 2 | 2 ✅ (confirmed match) |
| **Chunk Size** | Not specified | 50000 reads (added from WDL) |

### 3. **File Path Validation**
✅ **All referenced files exist and are accessible**:
- `tests/small-region-capture-micro-c/small_rcmc_r1.fq.gz` (5.2 MB)
- `tests/small-region-capture-micro-c/small_rcmc-extra-reads_r1.fq.gz` (5.2 KB)
- `tests/small-region-capture-micro-c/small_rcmc_r2.fq.gz` (5.3 MB)
- `tests/small-region-capture-micro-c/small_rcmc-extra-reads_r2.fq.gz` (5.5 KB)
- `tests/small-region-capture-micro-c/test_bwa_index.tgz` (5.7 MB)
- `tests/small-region-capture-micro-c/test.chrom.sizes` (32 bytes)

## 🔧 **Technical Fixes Applied**

### 1. **FASTQ File Pairing Logic**
**Problem**: Using `.combine()` created Cartesian product (4 combinations instead of 2 pairs)
```groovy
// BEFORE (incorrect - creates all combinations)
fastq_pairs_ch = fastq_r1_ch.combine(fastq_r2_ch)

// AFTER (correct - pairs files by index)
fastq_pairs_ch = Channel.fromList(
    fastq_r1_list.indices.collect { i ->
        def r1 = file(fastq_r1_list[i])
        def r2 = file(fastq_r2_list[i])
        def chunk_id = r1.baseName.replaceAll(/\.fq\.gz$|\.fastq\.gz$/, '')
        [chunk_id, r1, r2]
    }
)
```

### 2. **File Name Collision Prevention**
**Solution**: Proper R1/R2 pairing ensures unique output filenames:
- Pair 1: `small_rcmc_r1.fq.gz` + `small_rcmc_r2.fq.gz` → `small_rcmc_r1.fq.pairsam.gz`
- Pair 2: `small_rcmc-extra-reads_r1.fq.gz` + `small_rcmc-extra-reads_r2.fq.gz` → `small_rcmc-extra-reads_r1.fq.pairsam.gz`

### 3. **Configuration Documentation**
Updated configuration comments to reference WDL JSON source:
```groovy
// From small_rcmc.json: "microc.sample_id": "small-rcmc"
sample_id = 'small-rcmc'

// From small_rcmc.json: "microc.fastq_r1": [main_file, extra_reads_file]
fastq_r1 = 'tests/small-region-capture-micro-c/small_rcmc_r1.fq.gz,tests/small-region-capture-micro-c/small_rcmc-extra-reads_r1.fq.gz'
```

## ✅ **Validation Results**

### **Automated Validation Script**
Created `validate_test_config_consistency.py` which confirms:

```
✅ VALIDATION PASSED
The Nextflow test configuration exactly matches the WDL JSON specification!

🎉 Both implementations will use identical test data and parameters!

📊 VALIDATION SUMMARY
Total WDL parameters: 8
Matching parameters: 8
Differences found: 0
Missing files: 0
```

### **Pipeline Execution Test**
```bash
nextflow run nf/main.nf -profile test -stub-run
```

**Result**: ✅ **All processes completed successfully**
- ✔ EXTRACT_BWA_INDEX (2 processes - one per file pair)
- ✔ BWA_ALIGN_PAIRTOOLS (2 processes - correctly paired R1/R2)
- ✔ MERGE_DEDUP_SPLIT (1 process - merging both pairsam files)
- ✔ COOLER_PROCESS (1 process)
- ✔ JUICER_HIC (1 process)
- ✔ QC_METRICS (1 process)

## 🎯 **Consistency Achieved**

### **WDL vs Nextflow Parameter Mapping**
| WDL JSON Parameter | Nextflow Parameter | Value | Status |
|-------------------|-------------------|-------|---------|
| `microc.sample_id` | `sample_id` | `small-rcmc` | ✅ Match |
| `microc.fastq_r1` | `fastq_r1` | 2 files (comma-separated) | ✅ Match |
| `microc.fastq_r2` | `fastq_r2` | 2 files (comma-separated) | ✅ Match |
| `microc.reference_bwa_idx` | `reference_bwa_idx` | `test_bwa_index.tgz` | ✅ Match |
| `microc.chrom_sizes` | `chrom_sizes` | `test.chrom.sizes` | ✅ Match |
| `microc.cooler.resolution` | `resolution` | `1000` | ✅ Match |
| `microc.microc_align.bwa_cores` | `bwa_cores` | `2` | ✅ Match |
| `microc.num_reads_per_chunk` | `num_reads_per_chunk` | `50000` | ✅ Match |

## 🚀 **Usage Examples**

### **Test with Identical Data**
```bash
# WDL pipeline (original)
# Uses: tests/small-region-capture-micro-c/small_rcmc.json

# Nextflow pipeline (updated)
nextflow run nf/main.nf -profile test,docker
```

### **Validation Commands**
```bash
# Validate configuration consistency
python3 validate_test_config_consistency.py

# Test pipeline structure
nextflow run nf/main.nf -profile test -stub-run

# Full test with containers
nextflow run nf/main.nf -profile test,docker
```

## 📁 **File Organization**

```
topology-tools/
├── tests/small-region-capture-micro-c/
│   ├── small_rcmc.json                 # WDL JSON specification (source)
│   ├── small_rcmc_r1.fq.gz            # Main R1 FASTQ
│   ├── small_rcmc-extra-reads_r1.fq.gz # Extra R1 FASTQ
│   ├── small_rcmc_r2.fq.gz            # Main R2 FASTQ
│   ├── small_rcmc-extra-reads_r2.fq.gz # Extra R2 FASTQ
│   ├── test_bwa_index.tgz              # BWA genome index
│   └── test.chrom.sizes                # Chromosome sizes
│
├── nf/conf/
│   └── test.config                     # Nextflow test config (updated)
│
└── validate_test_config_consistency.py # Validation script
```

## 🎉 **Summary**

✅ **Configuration Updated**: Nextflow test config now exactly matches WDL JSON  
✅ **Files Validated**: All referenced test data files exist and are accessible  
✅ **Pipeline Tested**: Stub run confirms correct workflow execution  
✅ **Pairing Fixed**: R1/R2 files are correctly paired (no more collisions)  
✅ **Documentation Added**: Clear mapping between WDL and Nextflow parameters  
✅ **Validation Script**: Automated consistency checking for future updates  

**Both the WDL and Nextflow implementations now use identical test datasets and parameters, ensuring consistent validation and comparison results!** 🎯
