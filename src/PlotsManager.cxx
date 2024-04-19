
#include "../include/PlotsManager.h"
//argv[0] = inputSource
//argv[1] = outputFile
int main(int argc, char *argv[]){ 
	
	cout<<"Starting the analysis..."<<endl;
	if (argc != 3 ){
		cout << "Incorrect number of arguments. First argument is the input source, second is the outputfile position."<<endl;
		return 1;
	}
	const char* inputPosition = argv[1];
	const char* outputPosition = argv[2];

	if (!defineAnalysis()){
		cout << "Cannot define analysis from the possible choices." << endl;
		return 1;
	}
	
   	//create and open output file
   	mOutFile = CreateOutputFile("AnalysisOutput.root"); 
   	if(!mOutFile) {
      	cout << "Can not open output file." << endl; 
      	return 1;
   	}
   	cout << "Output file created..." << endl;


    if (!ConnectInputHist(argc, argv)){
    	cout << "Couldn't connect all histograms." << endl;
    	return 1;
    }
    cout << "Successfully added all histograms together..." << endl;


    if( runAnaV0 ){
    	cout << "Creating plots from AnaV0..." << endl;
    	mPlot = new PlotAnaV0(mOutFile, inputPosition, outputPosition);
    }

    mPlot->Init();
    mPlot->Make();


    Clear();

    cout << "Ending analysis. All plots created. Goodbye..." << endl;


	return 0;
}

