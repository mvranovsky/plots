#!/usr/bin/python3


import sys
from subprocess import Popen, PIPE
import os

if __name__ == "__main__":

    args = sys.argv
    args.pop(0)
    if len(args) != 1:
        print("Usage: ProbRetainEvent.py <Production_name>")
        sys.exit(1)

    prod = args[0]

    # Check if the production name is valid
    top = os.path.join("/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" ,prod )

    if not os.path.exists(top):
        print("Production name does not exist")
        sys.exit(1)
    
    files = os.listdir(os.path.join(top, "*.txt"))

    if len(files) == 0:
        print("No files found in the production directory")
        sys.exit(1)

    # define a TTree and a TFile 
    # Create a TFile to store the results

    outData = ""
    for file in files:
        file_path = os.path.join(top, file)
        with open(file_path, "r") as f:
            lines = f.readlines()
        
        outData += lines[0] + "\n"
    
    with open("ProbRetainEvent.txt", "w") as f:
        f.write(outData)
    
    cmd = "root -l -b -q ProbRetainEvent.C(\"ProbRetainEvent.txt\")"
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
    


        



