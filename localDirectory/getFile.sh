#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <TString Argument>"
    exit 1
fi

echo "Downloading .root file from RCF with tag: $1 and saved as /$1/AnalysisOutput.root"
mkdir $1
#download the .root file from RCF using sftp
sftp mvranovsk@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/mvranovsk/Run17_P20ic/$1/merged/AnalysisOutput.root ./$1/AnalysisOutput.root

echo "Successfully downloaded .root file. Running script getData.C to save the downloaded canvases to directory $1/"

root -l -b -q "getData.C(\"$1\")"

echo "Successfully created all the canvases and saved them. Goodbye!"


