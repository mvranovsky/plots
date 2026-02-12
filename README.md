# plots
Program located at remote computing facility of Brookhaven National Laboratory for experiment STAR. It is used to process single or multiple root files and create plots and histograms.



things to include:
- introduction + outline
- objects: classes, pointers
- plots: how do different people handle that? specific code for each plot? control plots? 
- show schematic view of the program
- introduce how each part looks like
- toy example 
- adding class to CMakeLists.txt
- FitJPsi class -> Fit class
- CrossSectionMaker
- bichselG10?



veci na prerobenie v kode ked skompilujem:
- Plot::loadInvMassHist
- Plot::getCondition
- setTH2Style, TH1General, TH2General na templates



//#ifdef __CLING__
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Plot+;
#pragma link C++ class PlotAnalysis+;
#pragma link C++ class Fit+;

#endif
