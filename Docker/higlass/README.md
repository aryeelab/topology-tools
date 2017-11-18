
# Build higlass image that contains a copy of all mcool files in the PWD
docker build -t aryeelab/higlass-colon .

# Save the docker image to a tgz
docker save aryeelab/higlass-colon | gzip > higlass-colon.tar.gz

# Load the tgz docker image
gunzip -c higlass-colon.tar.gz | docker load

# Run the higlass image
docker run --detach --publish 8989:80 --name higlass-container aryeelab/higlass-colon
sleep 3

# Ingest the mcool files that are baked into the image
docker exec higlass-container sh ingest.sh


