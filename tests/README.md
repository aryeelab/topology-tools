# Topology tools WDL pipeline tests

The tests involve running the pipelines on small test datasets. 
The `small-region-capture-micro-c` dataset, for example, is created by `small-region-capture-micro-c/make_small_rcmc.sh`.
 
To use locally install make sure `pytest-workflow` is installed. 

	conda create --name topology-tools
	conda activate topology-tools
	conda install -c conda-forge pytest-workflow

To run tests:

	conda activate topology-tools # Once per shell session
	pytest --git-aware --keep-workflow-wd-on-fail
