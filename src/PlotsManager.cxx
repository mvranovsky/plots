#include "../include/PlotsManager.h"

int main(int argc, char *argv[]){ 

	cout<<"Beginning to make plots..."<<endl;
	if (argc != 4 ){
		cerr << "ERROR: main: Incorrect number of arguments. Usage ./PlotsManager <analysis_tag> <input_source> <output_path> "<<endl;
		return 1;
	}

    if(DEBUGMODE) cout << "Extracting arguments..." << endl;
    const char* analysis = argv[1];       // analysis tag, for now just Analysis, can add stuff like Embedding, Zerobias, TofEfficiency based on which separate class is run
	const char* inputPosition = argv[2];  // input source can be a single root file or a .list file where each line is an absolute path to a .root file
	const char* outputPosition = argv[3]; // absolute path to output position
    const char* histogramFilePath = "histFile.root";


    
    if( strstr(analysis, "Analysis") ){
    	if(DEBUGMODE)  << "Creating plots for Analysis..." << endl;
    	mPlot = new PlotAnalysis(inputPosition, outputPosition);
    }else if( strstr(analysis, "CrossSection")){  // separate class which can run multiple analyses and calculate CS for example
        if(DEBUGMODE)  cout << "Calculating cross section..." << endl;
        runCrossSection(outputPosition);
    }
    /*                                              
    else if( strstr(inputPosition, "Embedding") ){  // example of another analysis
        cout << "Creating plots for embedding..." << endl;
        mPlot = new PlotEmbedding(inputPosition, outputPosition);
    }
    */

    // extract histograms and merge them (if it is a list of separate root files)
    if (!connectHists(inputPosition, histogramFilePath)){
        cerr << "ERROR: main: Couldn't connect all histograms." << endl;
        return 1;
    }
    if(DEBUGMODE) cout << "Successfully added all histograms together..." << endl;

    mPlot->Init();

    mPlot->Make();

    mPlot->Finish();

    if(DEBUGMODE) cout << "Ending analysis. All plots created. Goodbye..." << endl;

	return 0;

}


void runCrossSection(const char* outputFile) {. // for running multiple analyses
    if(DEBUGMODE) cout << "Running CrossSectionMaker..." << endl;

    // paths to main analysis, embedding, good run list analysis, zerobias study
    TString anaDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaJPsi_noRP_sysStudy_2.12.25/merged/StRP_production.list";
    TString embedDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/EmbeddingJPsi_sysStudy/merged/StRP_production.list";
    TString goodRunDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaGoodRun_2.12.25/merged/StRP_production.list";
    TString zeroBiasDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaZeroBias_26.8.25/merged/StRP_production.list";

    // extract histograms from each separate root file    

    if(!connectHists(anaDir, "histFile.root", false))  return;

    if(!connectHists(embedDir, "histFile.root", true))  return;

    if(!connectHists(goodRunDir, "histFile.root", true))  return;

    if(!connectHists(zeroBiasDir, "histFile.root", true))  return;

    // run the cross section analysis
    CrossSectionMaker *crossSectionMaker = new CrossSectionMaker(anaDir, embedDir, goodRunDir, zeroBiasDir, TString(outputFile));
    crossSectionMaker->Make();
    delete crossSectionMaker;

    if(DEBUGMODE) cout << "Finished running CrossSectionMaker..." << endl;
}

