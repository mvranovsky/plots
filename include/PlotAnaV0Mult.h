#ifndef PlotAnaV0Mult_h
#define PlotAnaV0Mult_h


#include "Util.h"
#include "Plot.h"
#include "../../star-upcDst/work/include/RunDef.h"


using namespace std;
using namespace UTIL;

class PlotAnaV0Mult : public Plot {
	public:
		//accessible from outside the class
		PlotAnaV0Mult( const string mInputList, const char* filePath);
      	~PlotAnaV0Mult(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

    	void invMassK0s(Double_t minRange, Double_t maxRange, int nBins);
    	void invMassLambda(Double_t minRange, Double_t maxRange, int nBins);
    	void nSigmaCorr();
    	void histGeneral(TString nameOfHist, TString cutDescription);
    	void fitK0s(Double_t minRange, Double_t maxRange, int nBins);
      	void fitLambda(Double_t minRange, Double_t maxRange, int nBins);
		void general2f(TString nameOfHist , TString cutDescription);

		



};


#endif