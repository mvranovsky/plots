#include "PlotAnalysis.h"

PlotAnalysis::PlotAnalysis(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotAnalysis::PlotAnalysis(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotAnalysis::Make(){
    
    cout << "Number of objects in the histogram file: " << histFile->GetListOfKeys()->GetEntries() << endl;

    for(int i = 0; i < histFile->GetListOfKeys()->GetEntries(); ++i){
        TObject *obj = (TObject*)histFile->GetListOfKeys()->At(i);
        cout << "Found object: " << obj->GetName() << endl;
        if(obj->InheritsFrom(TH1::Class())){
            cout << "This is a 1D histogram." << endl;
        }
    }
    // save the extracted histograms to the output file
    saveHistograms();
    
    // load inv mass histogram
    TH1D* invMass = (TH1D*) histFile->Get("hSub;2");

    if(!invMass || invMass->IsZombie()){
        cerr << "ERROR in PlotAnalysis::Make(): Could not load invariant mass histogram. Leaving..." << endl;
        return;
    }

    // fit peak with Fit class
    Fit *fit = new Fit(invMass,"Gauss + Poly1");

    // set initial parameters for the fit
    fit->setMean(1.02);
    fit->setSigma(0.05);
    fit->setPolynomial(10, 1, -0.5);

    fit->fitPeak();
    fit->writeFitResult();

    // save the output histogram
    fit->saveCanvas(outFile,"", "");

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


