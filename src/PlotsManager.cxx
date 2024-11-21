
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


   	//create and open output file
   	mOutFile = CreateOutputFile("AnalysisOutput.root"); 
   	if(!mOutFile) {
      	cout << "Cannot open output file." << endl; 
      	return 1;
   	}
   	cout << "Output file created..." << endl;



    if (!connectHists(argc, argv)){
    	cout << "Couldn't connect all histograms." << endl;
    	return 1;
    }
    cout << "Successfully added all histograms together..." << endl;


    if( runAnaV0 ){
    	cout << "Creating plots from AnaV0..." << endl;
    	mPlot = new PlotAnaV0(mOutFile, inputPosition, outputPosition);
    }else if( runAnaV0SingleState ){
    	cout << "Creating plots from AnaV0 Single state..." << endl;
    	mPlot = new PlotV0SingleState(mOutFile, inputPosition, outputPosition);
    } else if( runTofEff ){
    	cout << "Creating plots from TofEff..." << endl;
    	mPlot = new PlotTofEff(mOutFile, inputPosition, outputPosition);
    } else if( runTofEffMult ){
        cout << "Creating plots from TofEffMult..." << endl;
        mPlot = new PlotTofEffMult(mOutFile, inputPosition, outputPosition);
    }else{
    	cout << "No plots to run. Leaving..." << endl;
    	Clear();
    	return 1;
    }

    mPlot->Init();
    mPlot->Make();


    Clear();

    cout << "Ending analysis. All plots created. Goodbye..." << endl;


	return 0;
}

