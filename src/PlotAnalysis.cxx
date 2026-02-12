#include "PlotAnalysis.h"

PlotAnalysis::PlotAnalysis(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotAnalysis::PlotAnalysis(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotAnalysis::Make(){

    //EXERCISE
    // save all histograms to canvases and save them to the output file


    // load invariant mass histogram with name hSub


    //use Fit class to fit the peak in the invariant mass histogram and extract the yield of phi mesons


    // save the final canvas to the output file
    // output file should contain only final filled canvases ready to be saved as .pdf/.png files


}

void PlotAnalysis::Finish(){
    
    if(outFile) outFile->Close();
    if(histFile) histFile->Close();
    if(MCFile) MCFile->Close();
    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE) cout << "The output file is saved: " << outputPosition << endl;

}

void PlotAnalysis::Init(){
    if(DEBUGMODE) cout << "Initializing PlotAnalysis..." << endl;

    // setup the drawing style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // open the output file
	if(!outFile || outFile->IsZombie() ){
		cerr << "ERROR in PlotAnalysis::Init(): Couldn't open output file. Leaving..." << endl;
        return;
	}

    outFile->cd();
    // load file with histograms
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "ERROR in PlotAnalysis::Init(): Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
    /*
	ConnectInputTree(inputPosition, nameOfTree, tree);
    
    if(!tree || tree->IsZombie() || tree->GetEntries() == 0){
    	cerr << "ERROR in PlotAnalysis::Init(): Couldn't open tree with data. Returning." << endl;
    	return;
    }
    */

    // if background tree is needed, load it here the same way as the main one


    if(DEBUGMODE) cout << "Initialization finished successfully." << endl;

}






















/*
    // save the extracted histograms to the output file
    saveHistograms();
    
    // load inv mass histogram
    TH1D* invMass = (TH1D*) histFile->Get("hSub");
    
    if(!invMass || invMass->IsZombie()){
        cerr << "ERROR in PlotAnalysis::Make(): Could not load invariant mass histogram. Leaving..." << endl;
        return;
    }


    invMass->SetName("invariantMassFit");
    invMass->SetTitle("");

    // fit peak with Fit class
    Fit *fit = new Fit(invMass,"Gauss Poly1");
    fit->setFitRangeLow(1.005);
    fit->setFitRangeHigh(1.06);

    // set initial parameters for the fit
    fit->setMean(1.02);
    fit->setSigma(0.01);
    fit->setPolynomial( 0, 10);

    fit->fitPeak();
    fit->writeFitResult();

    // save the output histogram
    TCanvas *C = fit->getCanvas();
    C->SetName("invariantMassFitCanvas");
    DrawSTARTag();
    C->Update();
    outFile->cd();
    C->Write();
*/
