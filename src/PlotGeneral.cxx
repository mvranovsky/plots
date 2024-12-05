#include "../include/PlotGeneral.h"


PlotGeneral::PlotGeneral(TFile *mOutFile, const string mInputList, const char* filePath): Plot(mOutFile, mInputList, filePath){}

void PlotGeneral::Make(){


	hists1D = GetAllTH1D();
	if(hists1D.size() == 0){
		cerr << "Couldn't load 1D histograms from file. Leaving..." << endl;
		return;
	}
	for (int i = 0; i < hists1D.size(); ++i){
		if(!hists1D[i].first){
			cerr << "Couldn't load histogram " << hists1D[i].second << ". Leaving..." << endl;
			return;
		}
		TH1DGeneral(hists1D[i].second, hists1D[i].first);
		cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
	}




	hists2D = GetAllTH2F();
	if(hists2D.size() == 0){
		cerr << "Couldn't load 2D histograms from file. Leaving..." << endl;
		return;
	}
	for (int i = 0; i < hists2D.size(); ++i){
		if(!hists2D[i].first){
			cerr << "Couldn't load histogram " << hists2D[i].second << ". Leaving..." << endl;
			return;
		}
		TH2FGeneral(hists2D[i].second, hists2D[i].first);
		cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
	}




	cout << "All histograms successfully saved to canvases..." << endl;
	cout << "The output file is saved: " << outputPosition << endl;

}

void PlotGeneral::Init(){
	//define the output file which will store all the canvases
	outFile = new TFile(outputPosition, "recreate");

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaV0Tree);

    if(!tree){
    	cerr << "Couldn't open tree with data. Returning." << endl;
    	return;
    }
}



