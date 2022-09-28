# topology-tools

## Setup

#### Install cromwell

On Macs Cromwell can be installed with 
```
brew install cromwell
```

Otherwise you can download the Cromwell Jar from the [releases](https://github.com/broadinstitute/cromwell/releases) page.

#### Set up Google Cloud Tools for access to docker images:

First download and install the Google Cloud CLI: https://cloud.google.com/sdk/docs/install-sdk

Then authenticate with Google for access to the docker image registry

	gcloud auth login
	gcloud auth configure-docker us-central1-docker.pkg.dev

## Running the WDL workflows in Cromwell

#### Micro-C

```
cromwell run -i tests/small-region-capture-micro-c/small_rcmc.json microc.wdl
```

#### Hi-C

```
./cromwell run preprocess_hic.wdl imr90-rep1_hic_small.json 
./cromwell run preprocess_hic.wdl imr90-rep2_hic_small.json 
```
    
## Using a database to store workflow state and history

Using a MySQL database allows Cromwell to persist job state and history between restarts. This enables call caching which can vastly speed up multiple runs (which is especially useful during development/debugging).

First start a mysql docker container to host the cromwell database:

    docker run -p 3306:3306 --name cromwell_mysql -e MYSQL_ROOT_PASSWORD=root_pass \
    -e MYSQL_DATABASE=cromwell -e MYSQL_USER=cromwell -e MYSQL_PASSWORD=cromwell_pass \
    -d mysql/mysql-server:5.5

A config file (e.g. `cromwell.local.conf`) specifies Cromwell options, including the database connection and usernames/passwords. Specify the config file when running Cromwell like this:

	java -Dconfig.file=cromwell.local.conf \
	-jar /usr/local/Cellar/cromwell/81/libexec/cromwell.jar \
	run -i tests/small-region-capture-micro-c/small_rcmc.json microc.wdl

Or, (if using Cromwell installed by homebrew on a Mac):

	export JAVA_OPTS=-Dconfig.file=cromwell.local.conf 
	cromwell run -i tests/small-region-capture-micro-c/small_rcmc.json microc.wdl


## Building or pulling docker images

Docker images are built automatically using GitHub Actions and GCP Cloud Build on every push to GitHub. You can pull the latest images from the container registry:

```
docker pull us-central1-docker.pkg.dev/aryeelab/docker/microc
docker pull us-central1-docker.pkg.dev/aryeelab/docker/juicer
```

You can see the status of image builds (completion, success, logs etc.) on the [GitHub Actions](https://github.com/aryeelab/topology-tools/actions) page.

You can alternatively build images locally:

```
cd Docker/microc
docker build -t us-central1-docker.pkg.dev/aryeelab/docker/microc .
```



## Configuring Workload Identity Federation for Github <-> GCP auth

	gcloud iam workload-identity-pools create "github-actions-pool" \
	  --project="aryeelab" \
	  --location="global" \
	  --display-name="Github Actions Pool"

	gcloud iam workload-identity-pools providers create-oidc "github-provider" \
	  --project="aryeelab" \
	  --location="global" \
	  --workload-identity-pool="github-actions-pool" \
	  --display-name="Github provider" \
	  --attribute-mapping="google.subject=assertion.sub,attribute.actor=assertion.actor,attribute.aud=assertion.aud" \
	  --issuer-uri="https://token.actions.githubusercontent.com"

	gcloud iam service-accounts add-iam-policy-binding "github-actions@aryeelab.iam.gserviceaccount.com" \
	  --project="aryeelab" \
	  --role="roles/iam.workloadIdentityUser" \
	  --member="principalSet://iam.googleapis.com/projects/303574531351/locations/global/workloadIdentityPools/github-actions-pool/*"