#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <TString Argument>"
    exit 1
fi

./RunPlot.py $1 1> output.out 2>output.err &
