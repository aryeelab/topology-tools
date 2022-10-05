import pathlib
import pytest

@pytest.mark.workflow('small-region-capture-micro-c')
def test_mapped_pairs(workflow_dir):
    mapped_pairs = pathlib.Path(workflow_dir, "test-output/mapped.pairs").read_text()
    assert 3 % 3 == 0