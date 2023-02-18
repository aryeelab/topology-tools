import os
import pytest
import hashlib

# This test checks the md5 of mapped.pairs but ignores the header lines
# which contain run-specific FASTQ paths (which would change the md5)

@pytest.mark.workflow('small-region-capture-micro-c')
def test_mapped_pairs(workflow_dir):
    mapped_pairs_file = os.path.join(workflow_dir, "test-output/small-rcmc.mapped.pairs")
    f = open(mapped_pairs_file, "r")
    pairs = f.readlines()[19:]
    f.close()
    pairs.sort()
    pairs_string = "".join(pairs)
    md5 = hashlib.md5(pairs_string.encode()).hexdigest()
    assert md5 == '98340fc4018d9fc2dee42c14be973137'


