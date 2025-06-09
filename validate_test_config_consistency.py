#!/usr/bin/env python3
"""
Validation script to ensure Nextflow test configuration matches WDL JSON specification.
This script compares the test parameters in nf/conf/test.config with tests/small-region-capture-micro-c/small_rcmc.json
"""

import json
import os
import sys
from pathlib import Path

def parse_nextflow_config(config_path):
    """Parse Nextflow test configuration to extract parameters."""
    params = {}
    
    with open(config_path, 'r') as f:
        content = f.read()
    
    # Extract parameters from the config file
    in_params_block = False
    for line in content.split('\n'):
        line = line.strip()
        
        if line.startswith('params {'):
            in_params_block = True
            continue
        elif line == '}' and in_params_block:
            in_params_block = False
            continue
        
        if in_params_block and '=' in line and not line.startswith('//'):
            # Parse parameter lines like: sample_id = 'small-rcmc'
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip().rstrip(',').strip("'\"")
            params[key] = value
    
    return params

def parse_wdl_json(json_path):
    """Parse WDL JSON configuration to extract parameters."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Extract relevant parameters from WDL JSON
    wdl_params = {}
    
    # Map WDL parameter names to our parameter names
    param_mapping = {
        'microc.sample_id': 'sample_id',
        'microc.fastq_r1': 'fastq_r1',
        'microc.fastq_r2': 'fastq_r2',
        'microc.reference_bwa_idx': 'reference_bwa_idx',
        'microc.chrom_sizes': 'chrom_sizes',
        'microc.cooler.resolution': 'resolution',
        'microc.microc_align.bwa_cores': 'bwa_cores',
        'microc.num_reads_per_chunk': 'num_reads_per_chunk'
    }
    
    for wdl_key, nf_key in param_mapping.items():
        if wdl_key in data:
            value = data[wdl_key]
            if isinstance(value, list):
                # Convert list to comma-separated string for comparison
                value = ','.join(value)
            wdl_params[nf_key] = str(value)
    
    return wdl_params

def validate_file_paths(params, base_dir):
    """Validate that all file paths in parameters exist."""
    file_params = ['fastq_r1', 'fastq_r2', 'reference_bwa_idx', 'chrom_sizes']
    missing_files = []
    
    for param in file_params:
        if param in params:
            if ',' in params[param]:
                # Multiple files
                files = params[param].split(',')
            else:
                files = [params[param]]
            
            for file_path in files:
                full_path = os.path.join(base_dir, file_path)
                if not os.path.exists(full_path):
                    missing_files.append(full_path)
    
    return missing_files

def compare_parameters(nf_params, wdl_params):
    """Compare Nextflow and WDL parameters."""
    differences = []
    matches = []
    
    # Check each WDL parameter
    for key, wdl_value in wdl_params.items():
        if key in nf_params:
            nf_value = nf_params[key]
            if nf_value == wdl_value:
                matches.append(f"✅ {key}: {nf_value}")
            else:
                differences.append(f"❌ {key}: NF='{nf_value}' vs WDL='{wdl_value}'")
        else:
            differences.append(f"❌ {key}: Missing in Nextflow config (WDL value: '{wdl_value}')")
    
    # Check for extra Nextflow parameters
    for key, nf_value in nf_params.items():
        if key not in wdl_params and key not in ['config_profile_name', 'config_profile_description', 'max_cpus', 'max_memory', 'max_time', 'output_dir', 'mapq']:
            differences.append(f"⚠️  {key}: Extra in Nextflow config (value: '{nf_value}')")
    
    return matches, differences

def main():
    print("🔍 Validating Nextflow test configuration against WDL JSON specification...")
    print("="*80)
    
    # Define paths
    base_dir = Path.cwd()
    nf_config_path = base_dir / "nf" / "conf" / "test.config"
    wdl_json_path = base_dir / "tests" / "small-region-capture-micro-c" / "small_rcmc.json"
    
    # Check if files exist
    if not nf_config_path.exists():
        print(f"❌ ERROR: Nextflow test config not found: {nf_config_path}")
        sys.exit(1)
    
    if not wdl_json_path.exists():
        print(f"❌ ERROR: WDL JSON config not found: {wdl_json_path}")
        sys.exit(1)
    
    print(f"📁 Nextflow config: {nf_config_path}")
    print(f"📁 WDL JSON config: {wdl_json_path}")
    print()
    
    # Parse configurations
    try:
        nf_params = parse_nextflow_config(nf_config_path)
        print(f"✅ Parsed Nextflow config: {len(nf_params)} parameters")
    except Exception as e:
        print(f"❌ ERROR parsing Nextflow config: {e}")
        sys.exit(1)
    
    try:
        wdl_params = parse_wdl_json(wdl_json_path)
        print(f"✅ Parsed WDL JSON config: {len(wdl_params)} parameters")
    except Exception as e:
        print(f"❌ ERROR parsing WDL JSON config: {e}")
        sys.exit(1)
    
    print()
    
    # Validate file paths
    print("📂 Validating file paths...")
    missing_files = validate_file_paths(nf_params, base_dir)
    if missing_files:
        print("❌ Missing files:")
        for file_path in missing_files:
            print(f"   - {file_path}")
        print()
    else:
        print("✅ All referenced files exist")
        print()
    
    # Compare parameters
    print("🔄 Comparing parameters...")
    matches, differences = compare_parameters(nf_params, wdl_params)
    
    if matches:
        print("✅ Matching parameters:")
        for match in matches:
            print(f"   {match}")
        print()
    
    if differences:
        print("❌ Parameter differences:")
        for diff in differences:
            print(f"   {diff}")
        print()
    
    # Summary
    print("="*80)
    print("📊 VALIDATION SUMMARY")
    print("="*80)
    
    total_wdl_params = len(wdl_params)
    matching_params = len(matches)
    
    print(f"Total WDL parameters: {total_wdl_params}")
    print(f"Matching parameters: {matching_params}")
    print(f"Differences found: {len(differences)}")
    print(f"Missing files: {len(missing_files)}")
    
    if differences or missing_files:
        print("\n❌ VALIDATION FAILED")
        print("The Nextflow test configuration does not exactly match the WDL JSON specification.")
        sys.exit(1)
    else:
        print("\n✅ VALIDATION PASSED")
        print("The Nextflow test configuration exactly matches the WDL JSON specification!")
        print("\n🎉 Both implementations will use identical test data and parameters!")
    
    # Show key configuration details
    print("\n📋 Key Configuration Details:")
    print(f"   Sample ID: {nf_params.get('sample_id', 'N/A')}")
    print(f"   FASTQ R1 files: {nf_params.get('fastq_r1', 'N/A').count(',') + 1}")
    print(f"   FASTQ R2 files: {nf_params.get('fastq_r2', 'N/A').count(',') + 1}")
    print(f"   Resolution: {nf_params.get('resolution', 'N/A')} bp")
    print(f"   BWA cores: {nf_params.get('bwa_cores', 'N/A')}")

if __name__ == "__main__":
    main()
