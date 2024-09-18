#!/bin/bash


get() {
    sftp mvranovsk@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/mvranovsk/Run17_P20ic/$1/merged/AnalysisOutput.root ./$1/AnalysisOutput.root
}

run() {
	
	ROOT_MACRO_ARGUMENT="$1"
	#echo "$ROOT_MACRO_ARGUMENT"
	root -l -q "getData.C(\"$ROOT_MACRO_ARGUMENT\")"

}

runFit() {
	ROOT_MACRO_ARGUMENT1="$1"
	#ROOT_MACRO_ARGUMENT2="fitPeakK0sSignalFit"
	#ROOT_MACRO_ARGUMENT3="invMassLambdaSignalFit"
	root -q -l "FittingWithRooFit.C(\"$ROOT_MACRO_ARGUMENT1\")"
}


runEff() {
	ROOT_MACRO_ARGUMENT1="$1"
	root -q -l "efficiency.C(\"$ROOT_MACRO_ARGUMENT1\")"
}

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <TString Argument>"
    exit 1
fi

echo "Downloading .root file from RCF with tag: $1 and saved as /$1/AnalysisOutput.root"
mkdir $1
get $1
echo "Successfully downloaded .root file. Running script getData.C to save the downloaded canvases to directory $1/"
run $1
echo "Successfully created all the canvases and saved them. Now it is time for RooFit."
runFit $1
#echo "Created all plots using RooFit. Now it is time for efficiencies."
#runEff $1

