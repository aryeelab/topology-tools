# topology-tools

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

### Running the WDL Hi-C workflow in Cromwell

Preprocess HiC:
```
./cromwell run preprocess_hic.wdl imr90-rep1_hic_small.json 
./cromwell run preprocess_hic.wdl imr90-rep2_hic_small.json 
```
    
### Using a database to store workflow state and history

Using a MySQL database allows Cromwell to persist job state and history between restarts. This enables call caching which can vastly speed up multiple runs (which is especially useful during development/debugging).

First start a mysql docker container to host the cromwell database:

    docker run -p 3306:3306 --name cromwell_mysql -e MYSQL_ROOT_PASSWORD=root_pass \
    -e MYSQL_DATABASE=cromwell -e MYSQL_USER=cromwell -e MYSQL_PASSWORD=cromwell_pass \
    -d mysql/mysql-server:5.5

A config file specifies cromwell options, including the database connection and usernames/passwords.

	java -Dconfig.file=cromwell.local.conf \
	-jar /usr/local/Cellar/cromwell/79/libexec/cromwell.jar \
	run -i imr90-rep1_hic_small.json preprocess_hic.wdl


# Merging reads across multiple samples
```
./cromwell run merge_hic_replicates.wdl merge_hic_replicates_imr90.json 
```
