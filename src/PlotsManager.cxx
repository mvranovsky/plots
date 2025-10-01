#include "../include/PlotsManager.h"
//argv[0] = inputSource
//argv[1] = outputFile
int main(int argc, char *argv[]){ 
	


	cout<<"Starting making plots..."<<endl;
	if (argc != 3 ){
		cout << "Incorrect number of arguments. First argument is the input source, second is the outputfile position."<<endl;
		return 1;
	}
	const char* inputPosition = argv[1];
	const char* outputPosition = argv[2];
    const char* histogramFilePath = "histFile.root";

    if(strstr(argv[1], "-runCS")){
        runCrossSection(outputPosition);
        return 0;
    }
    
    
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
    } else if( strstr(inputPosition, "AnaJPsi") || strstr(inputPosition, "SysStudy") ){
        cout << "Creating plots from AnaJPsi..." << endl;
        mPlot = new PlotAnaJPsi(inputPosition, outputPosition);
    } else if( strstr(inputPosition, "EmbeddingJPsi") || strstr(inputPosition, "SysStudyEmbedding") ){
        cout << "Creating plots from EmbeddingJPsi..." << endl;
        mPlot = new PlotEmbeddingJPsi(inputPosition, outputPosition);
    }else if( strstr(inputPosition, "AnaGoodRun") ){
        cout << "Creating plots from AnaGoodRun" << endl;
        mPlot = new PlotGoodRun(inputPosition, outputPosition);
    }else if( strstr(inputPosition, "AnaZeroBias") ){
        cout << "Creating plots from AnaZeroBias..." << endl;
        mPlot = new PlotZeroBias(inputPosition, outputPosition);
    }else if( strstr(inputPosition, "BemcEfficiency") ){
        cout << "Creating plots from BemcEfficiency..." << endl;
        mPlot = new PlotBemcEfficiency(inputPosition, outputPosition);
    }else{
    	cout << "No plots to run. Leaving..." << endl;
    	return 1;
    }
    
    if (!connectHists(inputPosition, histogramFilePath)){
        cout << "Couldn't connect all histograms." << endl;
        return 1;
    }
    cout << "Successfully added all histograms together..." << endl;

    
    mPlot->Init();
    mPlot->Make();
    mPlot->Finish();

    cout << "Ending analysis. All plots created. Goodbye..." << endl;

	return 0;

}

