# Topology tools WDL pipeline tests

The tests involve running the pipelines on small test datasets. 
The `small-region-capture-micro-c` dataset, for example, is created by `small-region-capture-micro-c/make_small_rcmc.sh`.
 
To use locally install make sure `pytest-workflow` is installed. 

	mamba create --name topology-tools
	mamba activate topology-tools
	mamba install -c conda-forge pytest-workflow

To run tests:

	# From the topology-tools directory
	mamba activate topology-tools # Once per shell session
	pytest --git-aware --keep-workflow-wd-on-fail
