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
    	if(DEBUGMODE)  cout << "Creating plots for Analysis..." << endl;
    	mPlot = new PlotAnalysis(inputPosition, outputPosition);
    }
    /*// example of another analysis                                              
    else if( strstr(inputPosition, "Embedding") ){  
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


