#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <bitset>
#include "TString.h"


// if declaring variables here, make sure they are declared as const, otherwise the program will not compile

const bool DEBUGMODE = true; // true = print debug messages

// description shown in top right corner of every plot
//const TString plotDescription = "STAR Internal";
const TString plotDescription = "THIS THESIS";

// process which i am studying (will be shown below description)
const TString process = "Au + Au #rightarrow #phi (KK)";

// energy of collision at which I am studying the process (will be shown below process)
const TString energyOfCollision = "#sqrt{s_{NN}} = 200 GeV";


// if comparing results to an MC sample, declare the name of the file (with full path) and the name of the tree here
const TString MCFileName = "/path/to/MCfile.root";
const TString MCTreeName = "nameOfTreeForMCSimulations";


// settings for plot design
const int widthTypical = 1000;
const int heightTypical = 800;
const TString yAxisTitle = "counts"; //entries/counts
const int textAlign = 12;
const int textFont = 42; 
const Double_t textSize = 0.03;
const Double_t labelSize = 0.05;
const Double_t ticksSize = 0.05;
const int markerStyle = 20;	
const int lineStyle = 1;
const int lineColor = kBlue;
const int color = kBlue;
const int lineWidth = 2;


// settings for various analyses




#endif //ifndef CONFIG_H