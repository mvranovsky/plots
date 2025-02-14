
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
    const char* histogramFilePath = "histFile.root";


    if (!connectHists(inputPosition, histogramFilePath)){
    	cout << "Couldn't connect all histograms." << endl;
    	return 1;
    }
    cout << "Successfully added all histograms together..." << endl;


    if( strstr(inputPosition, "AnaV0Mult") ){
    	cout << "Creating plots from AnaV0Mult..." << endl;
    	mPlot = new PlotAnaV0Mult(inputPosition, outputPosition);
    }else if( strstr(inputPosition, "AnaV0" ) ){
    	cout << "Creating plots from AnaV0..." << endl;
    	mPlot = new PlotAnaV0(inputPosition, outputPosition);
    } else if( strstr(inputPosition, "TofEffMult") ){
        cout << "Creating plots from TofEffMult..." << endl;
        mPlot = new PlotTofEffMult(inputPosition, outputPosition);
    } else if( strstr(inputPosition, "TofEff") ){
    	cout << "Creating plots from TofEff..." << endl;
    	mPlot = new PlotTofEff(inputPosition, outputPosition);
    } else if( strstr(inputPosition, "AnaJPsi") || strstr(inputPosition, "AnaJPSI") ){
        cout << "Creating plots from AnaJPsi..." << endl;
        mPlot = new PlotAnaJPsi(inputPosition, outputPosition);
    }else{
    	cout << "No plots to run. Leaving..." << endl;
    	Clear();
    	return 1;
    }

    //mPlot->Init();
    //mPlot->Make();



    cout << "Ending analysis. All plots created. Goodbye..." << endl;

    Clear();

	return 0;

}

