#!/usr/bin/python 

import os
from subprocess import Popen, PIPE
import sys


# function which will collect information about average tracks in TOF, BEMC, TPC and save them to a single file


if __name__ == "__main__":

    prodDir = sys.argv[1]

    dir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + prodDir

    cmd = "ls " + dir + "/*.txt"
    print cmd
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()

    out = out[0].split("\n")
    data = ""
    i = 1
    for file_path in out[:-1]:
        print file_path

        with open(file_path, "r") as file:
            lines = file.readlines()
        data += str(i) + " "
        i += 1
        for line in lines:
            data += line.replace("\n", " ")
        
        data += "\n"
    
    print data
    with open("average_tracks.txt", "w") as file:
        file.write(data)
    
        
        



        
    



