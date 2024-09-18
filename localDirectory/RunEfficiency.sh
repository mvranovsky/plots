#!/bin/bash


runFit() {
	
	ROOT_MACRO_ARGUMENT1="$1"
    root -b -q -l "efficiency.C(\"$ROOT_MACRO_ARGUMENT1\")"

}

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <TString Argument>"
    exit 1
fi


mkdir "eff_$1"
runFit $1