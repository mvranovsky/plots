#ifndef PlotEmbeddingJPsi_h
#define PlotEmbeddingJPsi_h

#include "Util.h"
#include "Plot.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class PlotEmbeddingJPsi : public Plot {
	public:
		//accessible from outside the class
		PlotEmbeddingJPsi(const string mInputList, const char* filePath);
      	~PlotEmbeddingJPsi(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

    	void nSigmaCorr();
		double GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data);
		TString convertToString(double val);
		int makeInt(double val);
		void roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange, bool noBcg, bool plotDataPeak = false);
		void processTreeForPt(TString nameOfTree);
		void JPsiControlPlots();
		void nSigmaCorrPlot(int particles, bool justJPsi);  // 1 = electron, 2 = pion, 3 = kaon, 4 = proton
		void reconstructionEfficiency(int SWITCH, TString nameOfStarlightFile, TString nameOfOutput);  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
		void controlPlotsComparison(bool justJPsi);
		bool plot2Dists(TH1 *hData, TH1* hEmb, TString outName);


		Double_t lowLimInvMass, topLimInvMass;

};


#endif