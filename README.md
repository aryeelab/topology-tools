# topology_tools

### Setup

1. Install cromwell

On Macs cromwell can be installed with 
```
brew install cromwell
```

2. Clone this repository

3. Build the docker image(s)

```
cd Docker/hicpro
docker build -t aryeelab/hicpro .
```

### Running the WDL workflow in Cromwell

Preprocess HiC:
```
$ cromwell run preprocess_hic.wdl imr90-rep1_hic_small.json 
```