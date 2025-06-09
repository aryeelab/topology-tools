#!/usr/bin/env python3
"""
Test script for the Nextflow Micro-C pipeline.
Validates the pipeline structure, configuration, and basic functionality.
"""

import os
import sys
import subprocess
import json
from pathlib import Path

class NextflowPipelineValidator:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.info = []
        
    def add_error(self, message):
        self.errors.append(f"❌ ERROR: {message}")
        
    def add_warning(self, message):
        self.warnings.append(f"⚠️  WARNING: {message}")
        
    def add_info(self, message):
        self.info.append(f"ℹ️  INFO: {message}")
    
    def check_file_exists(self, filepath, description):
        """Check if a file exists"""
        if os.path.exists(filepath):
            self.add_info(f"{description}: {filepath} ✓")
            return True
        else:
            self.add_error(f"{description} missing: {filepath}")
            return False
    
    def validate_pipeline_structure(self):
        """Validate the Nextflow pipeline file structure"""
        self.add_info("Validating pipeline structure...")
        
        required_files = [
            ("main.nf", "Main workflow file"),
            ("nextflow.config", "Main configuration file"),
            ("conf/base.config", "Base configuration"),
            ("conf/modules.config", "Module configuration"),
            ("conf/test.config", "Test configuration"),
            ("modules/extract_bwa_index.nf", "BWA index extraction module"),
            ("modules/bwa_align_pairtools.nf", "BWA alignment module"),
            ("modules/merge_dedup_split.nf", "Merge/dedup module"),
            ("modules/cooler_process.nf", "Cooler processing module"),
            ("modules/juicer_hic.nf", "Juicer HiC module"),
            ("modules/qc_metrics.nf", "QC metrics module"),
            ("README_NEXTFLOW.md", "Nextflow documentation")
        ]
        
        all_exist = True
        for filepath, description in required_files:
            if not self.check_file_exists(filepath, description):
                all_exist = False
        
        return all_exist
    
    def validate_nextflow_syntax(self):
        """Validate Nextflow syntax"""
        self.add_info("Validating Nextflow syntax...")
        
        try:
            # Check if nextflow is available
            result = subprocess.run(['nextflow', '-version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                version_info = result.stdout.strip()
                self.add_info(f"Nextflow version: {version_info}")
            else:
                self.add_warning("Nextflow not found - syntax validation skipped")
                return False
        except Exception as e:
            self.add_warning(f"Could not check Nextflow version: {e}")
            return False
        
        # Validate main workflow syntax
        try:
            result = subprocess.run(['nextflow', 'config', 'main.nf'], 
                                  capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                self.add_info("Main workflow syntax is valid ✓")
                return True
            else:
                self.add_error(f"Syntax error in main.nf: {result.stderr}")
                return False
        except Exception as e:
            self.add_error(f"Error validating syntax: {e}")
            return False
    
    def validate_module_syntax(self):
        """Validate individual module syntax"""
        self.add_info("Validating module syntax...")
        
        modules = [
            "modules/extract_bwa_index.nf",
            "modules/bwa_align_pairtools.nf", 
            "modules/merge_dedup_split.nf",
            "modules/cooler_process.nf",
            "modules/juicer_hic.nf",
            "modules/qc_metrics.nf"
        ]
        
        all_valid = True
        for module in modules:
            if os.path.exists(module):
                # Basic syntax check - look for required elements
                with open(module, 'r') as f:
                    content = f.read()
                
                required_elements = ['process', 'input:', 'output:', 'script:']
                missing_elements = []
                
                for element in required_elements:
                    if element not in content:
                        missing_elements.append(element)
                
                if missing_elements:
                    self.add_error(f"Module {module} missing: {missing_elements}")
                    all_valid = False
                else:
                    self.add_info(f"Module {module} structure is valid ✓")
            else:
                self.add_error(f"Module {module} not found")
                all_valid = False
        
        return all_valid
    
    def validate_configuration(self):
        """Validate configuration files"""
        self.add_info("Validating configuration...")
        
        # Check main config
        if os.path.exists("nextflow.config"):
            with open("nextflow.config", 'r') as f:
                config_content = f.read()
            
            required_sections = ['params', 'profiles', 'process']
            missing_sections = []
            
            for section in required_sections:
                if section not in config_content:
                    missing_sections.append(section)
            
            if missing_sections:
                self.add_error(f"nextflow.config missing sections: {missing_sections}")
                return False
            else:
                self.add_info("nextflow.config structure is valid ✓")
        
        # Check if test data exists for test profile (relative to project root)
        test_files = [
            "../tests/small-region-capture-micro-c/small_rcmc_r1.fq.gz",
            "../tests/small-region-capture-micro-c/small_rcmc_r2.fq.gz",
            "../tests/small-region-capture-micro-c/test_bwa_index.tgz",
            "../tests/small-region-capture-micro-c/test.chrom.sizes"
        ]
        
        test_data_available = True
        for test_file in test_files:
            if not os.path.exists(test_file):
                self.add_warning(f"Test data missing: {test_file}")
                test_data_available = False
        
        if test_data_available:
            self.add_info("Test data is available ✓")
        
        return True
    
    def validate_containers_conda(self):
        """Validate container and conda specifications"""
        self.add_info("Validating container/conda specifications...")
        
        modules = [
            "modules/extract_bwa_index.nf",
            "modules/bwa_align_pairtools.nf", 
            "modules/merge_dedup_split.nf",
            "modules/cooler_process.nf",
            "modules/juicer_hic.nf",
            "modules/qc_metrics.nf"
        ]
        
        all_valid = True
        for module in modules:
            if os.path.exists(module):
                with open(module, 'r') as f:
                    content = f.read()
                
                has_conda = 'conda' in content
                has_container = 'container' in content
                
                if not (has_conda or has_container):
                    self.add_warning(f"Module {module} has no conda or container specification")
                    all_valid = False
                else:
                    self.add_info(f"Module {module} has environment specification ✓")
        
        return all_valid
    
    def run_dry_run_test(self):
        """Run a dry-run test of the pipeline"""
        self.add_info("Running dry-run test...")
        
        # Check if test data exists (relative to project root)
        test_files = [
            "../tests/small-region-capture-micro-c/small_rcmc_r1.fq.gz",
            "../tests/small-region-capture-micro-c/small_rcmc_r2.fq.gz",
            "../tests/small-region-capture-micro-c/test_bwa_index.tgz",
            "../tests/small-region-capture-micro-c/test.chrom.sizes"
        ]
        
        if not all(os.path.exists(f) for f in test_files):
            self.add_warning("Test data not available - skipping dry-run test")
            return False
        
        try:
            # Run nextflow with dry-run option
            cmd = [
                'nextflow', 'run', 'main.nf',
                '-profile', 'test',
                '-stub-run',  # Use stub mode for quick validation
                '--sample_id', 'test-validation',
                '--fastq_r1', '../tests/small-region-capture-micro-c/small_rcmc_r1.fq.gz',
                '--fastq_r2', '../tests/small-region-capture-micro-c/small_rcmc_r2.fq.gz',
                '--reference_bwa_idx', '../tests/small-region-capture-micro-c/test_bwa_index.tgz',
                '--chrom_sizes', '../tests/small-region-capture-micro-c/test.chrom.sizes',
                '--output_dir', 'test_validation_output'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0:
                self.add_info("Dry-run test completed successfully ✓")
                return True
            else:
                self.add_error(f"Dry-run test failed: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.add_error("Dry-run test timed out")
            return False
        except Exception as e:
            self.add_error(f"Error running dry-run test: {e}")
            return False
    
    def print_summary(self):
        """Print validation summary"""
        print("\n" + "="*80)
        print("🔍 NEXTFLOW PIPELINE VALIDATION SUMMARY")
        print("="*80)
        
        if self.info:
            print("\nℹ️  INFORMATION:")
            for msg in self.info:
                print(f"   {msg}")
        
        if self.warnings:
            print("\n⚠️  WARNINGS:")
            for msg in self.warnings:
                print(f"   {msg}")
        
        if self.errors:
            print("\n❌ ERRORS:")
            for msg in self.errors:
                print(f"   {msg}")
        
        print("\n" + "="*80)
        
        if self.errors:
            print("❌ VALIDATION FAILED - Please fix the errors above")
            return False
        elif self.warnings:
            print("⚠️  VALIDATION PASSED WITH WARNINGS - Pipeline should work but check warnings")
            return True
        else:
            print("✅ VALIDATION PASSED - Nextflow pipeline looks good!")
            return True

def main():
    print("🔍 Validating Nextflow Micro-C Pipeline...")
    
    validator = NextflowPipelineValidator()
    
    # Run all validation checks
    structure_ok = validator.validate_pipeline_structure()
    syntax_ok = validator.validate_nextflow_syntax()
    modules_ok = validator.validate_module_syntax()
    config_ok = validator.validate_configuration()
    containers_ok = validator.validate_containers_conda()
    
    # Only run dry-run if basic validation passes
    if structure_ok and syntax_ok:
        dry_run_ok = validator.run_dry_run_test()
    else:
        validator.add_warning("Skipping dry-run test due to validation errors")
        dry_run_ok = False
    
    # Print summary and exit
    success = validator.print_summary()
    
    if success:
        print("\n🎉 Nextflow pipeline validation completed successfully!")
        print("\nNext steps:")
        print("1. Run test: nextflow run nf/main.nf -profile test,docker")
        print("2. Run with your data: nextflow run nf/main.nf --sample_id ... --fastq_r1 ... --fastq_r2 ...")
        print("3. Check nf/README_NEXTFLOW.md for detailed usage instructions")
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
