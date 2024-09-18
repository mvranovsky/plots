#include <TFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TIterator.h>
#include <iostream>
#include <sstream> // For std::stringstream
#include <cstring>  // Include this for std::strcmp


using namespace std;



void compareHists(const char* name1, const char* name2){


	TFile *file1 = TFile::Open(name1 + TString("/AnalysisOutput.root"));
	TFile *file2 = TFile::Open(name2 + TString("/AnalysisOutput.root"));

	if(!file1 || !file2 ){
		cout << "Couldn't open one of the files. Leaving..." << endl;
		return;
	}

	TH1D* hist1 = (TH1D*)file1->Get("invMassK0sFit2");

    TH1D* hist2 = (TH1D*)file2->Get("invMassK0sFit2");

    TH1D *hist3 = new TH1D("hist3", "hist3",45, 0.45, 0.54 );
    if(!hist1 || !hist2){
    	cout << "Couldn't load one of the histograms. Leaving..." << endl;
    	return;
    }

    hist3->Add(hist1, 1);
    hist3->Add(hist2, -1);

    TCanvas *c = new TCanvas("comparison", "comparison", 800, 600);

    hist3->Draw("same");

    c->SaveAs("comparison.pdf");
 
}