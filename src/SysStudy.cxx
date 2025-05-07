#include "SysStudy.h"

SysStudy::SysStudy(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void SysStudy::Make(){


    // save all the histograms to canvases into outfile
    handleHistograms();
    

    outFile->Close();
    histFile->Close();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void SysStudy::Init(){
	//define the output file which will store all the canvases
	outFile = unique_ptr<TFile>(new TFile( outputPosition, "recreate") );

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}

    histFile = unique_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, true);


    if(!tree || !bcgTree){
    	cerr << "Couldn't open tree or background tree with data. Returning." << endl;
    	return;
    }

}







TString SysStudy::convertToString(double val) {

    ostringstream streamA;
    streamA << fixed << setprecision(1) << val;
    TString formattedA = streamA.str();

    return formattedA;
}

int SysStudy::makeInt(double val) {
	int result = val;
	return result;
}





