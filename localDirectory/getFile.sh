#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <RCF_absolute_path> <output_dir>"
    exit 1
fi

echo "Downloading .root file from RCF: $1 and saved as $2/AnalysisOutput.root"
mkdir $2
#download the .root file from RCF using sftp
sftp mvranovsk@sftp.sdcc.bnl.gov:$1 $2/AnalysisOutput.root

echo "Successfully downloaded .root file. Running script getData.C to save the downloaded canvases to directory $2/"

root -l -b -q "getData.C(\"$2\")"

echo "Successfully created all the canvases and saved them. Goodbye!"



