#!/bin/bash

for MCOOL in /mcool/*.mcool; do
    BASENAME=$(basename $MCOOL)
    SAMPLE_ID=${BASENAME%.*} # Remove extension
    echo $SAMPLE_ID
    python higlass-server/manage.py ingest_tileset --filename $MCOOL --name $SAMPLE_ID --datatype matrix --filetype cooler 
done
