#ifndef PlotV0SingleState_h
#define PlotV0SingleState_h


#include "Util.h"
#include "Plot.h"


using namespace std;
using namespace UTIL;

class PlotV0SingleState : public Plot {
	public:
		//accessible from outside the class
		PlotV0SingleState(TFile *mOutFile, const string mInputList, const char* filePath);
      	~PlotV0SingleState(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

    	void invMassK0s(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF);
    	void invMassLambda(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF);
    	void nSigmaCorr();
    	void fitK0s(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF);
      	void fitLambda(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF);

		
		vector<pair<TH1D*,TString>> hists1D;
		vector<pair<TH2F*,TString>> hists2D;

};


#endif