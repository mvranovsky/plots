#!/usr/bin/python

#--------------------------------------------------------------
# macro to run program plots
#
# usage:
# ./RunPlot.py tag 
# ./RunPlot.py tag inputsource
#
#
#--------------------------------------------------------------

import os
from subprocess import Popen
import re
import sys

#_____________________________________________________________________________

#name of config file from command line argument
args = sys.argv
inputSource = ""
tag = ""
inputDirectory = ""
args.pop(0) # cut first input from terminal = ./SubmitPlugin
if len(args) == 2:
    tag = args.pop(0) # read third input from terminal = tag (outputDir)
    inputSource = args.pop(0) # read second input from terminal = inputSource
elif len(args) == 1:
    tag = args.pop(0) # read third input from terminal = tag (outputDir)
    inputSource = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + tag + "/merged/StRP_production.list"
    #inputSource = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + tag + "/merged/StRP_production_0000.root" #set only to the first file, if running whole, switch to upper line->list

    inputDirectory = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + tag + "/merged"
else:
    print "Wrong input arguments"
    exit()

#run command
outputFile = inputDirectory + "/AnalysisOutput.root"
runCmd = "./PlotsManager " + inputSource + " " + outputFile
#print runCmd
Popen(runCmd, shell=True).communicate()
