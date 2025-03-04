#ifndef PlotAnaJPsi_h
#define PlotAnaJPsi_h


#include "Util.h"
#include "Plot.h"
#include "../../star-upcDst/work/include/RunDef.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class PlotAnaJPsi : public Plot {
	public:
		//accessible from outside the class
		PlotAnaJPsi(const string mInputList, const char* filePath);
      	~PlotAnaJPsi(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

    	void invMassJPsi(int numBins, Double_t minRange, Double_t maxRange, TString nameOfSignal, TString nameOfBcg);
    	void nSigmaCorr();
		double GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data);
		TString convertToString(double val);
		int makeInt(double val);
		void roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange);
		void processTreeForPt(TString nameOfTree);


};


#endif