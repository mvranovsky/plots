#!/usr/bin/python

# simple function which runs over values of average tracks in BEMC, TOF and clusters of BEMC, skips those that are out of defined range

if __name__ == "__main__":

    lowerLimTracks = 0.05
    upperLimTracks = 0.15
    lowerLimClusters = 1.0
    upperLimClusters = 1.6
    lowerLimVertices = 0.2
    upperLimVertices = 0.3

    # daq list
    with open("daqlist_nevt.txt", "r") as file:
        daq_all = file.readlines()
    
    with open("runNumberToBin.txt", "r") as file:
        run_number_to_bin = file.readlines()
    
    daq_chosen = []
    nEvents = 0
    for run in run_number_to_bin:
        run = run.split()
        runNumber = run[1]
        tracksBemc = float(run[3])
        clustersBemc = float(run[4])
        vertices = float(run[6])
        print "Run number: ", runNumber , " vertices: ", vertices
        if tracksBemc > upperLimTracks or tracksBemc < lowerLimTracks:
            continue
        if clustersBemc > upperLimClusters or clustersBemc < lowerLimClusters:
            continue
        if vertices > upperLimVertices or vertices < lowerLimVertices:
            continue
        # no condition on tof


        for line in daq_all:
            # check if the line contains the run number, if yes save it daq_chosen
            if runNumber in line:
                nEvents += int(line.split()[-1])
                daq_chosen.append(str(line))
                break
    
    with open("daqlist_embedding.txt", "w") as file:
        for line in daq_chosen:
            file.write(line) 

    print "Finished choosing daq files suitable for embedding. "
    print "New daq list: daqlist_embedding.txt"
    print "Number of events: ", nEvents




