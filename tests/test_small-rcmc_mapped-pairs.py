import os
import pytest
import hashlib

# This test checks the md5 of mapped.pairs but ignores the header lines
# which contain run-specific FASTQ paths (which would change the md5)

@pytest.mark.workflow('small-region-capture-micro-c')
def test_mapped_pairs(workflow_dir):
    mapped_pairs_file = os.path.join(workflow_dir, "test-output/mapped.pairs")
    f = open(mapped_pairs_file, "r")
    pairs = f.readlines()[12:]
    f.close()
    pairs_string = "".join(pairs)
    md5 = hashlib.md5(pairs_string.encode()).hexdigest()
    assert md5 == 'df7e1de6de6b10cf139d388d4d297526'
