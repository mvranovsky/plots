#ifndef PlotTofEffMult_h
#define PlotTofEffMult_h


#include "Util.h"
#include "Plot.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class PlotTofEffMult : public Plot {
	public:
		//accessible from outside the class
		PlotTofEffMult(TFile *mOutFile, const string mInputList, const char* filePath);
      	~PlotTofEffMult(){};
      
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
      	void effPlot(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF,bool isTag, TString name);
		double GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data);
		TString convertToString(double val);
      	void efficiency(int switcher, int runSeparatePions);
      	void efficiencyPhi();
		void efficiencyPt();
		void efficiencyVz();
 		vector<pair<int,double>> effEtaFit(double etaMin, double etaMax, Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2], Double_t nGuess[2], Double_t alphaGuess[2]);
 		vector<pair<int,double>> effPhiFit(double etaMin, double etaMax, Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2]);
 		vector<pair<int,double>> effPtFit(double etaMin, double etaMax, Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2]);
 		vector<pair<int,double>> effVzFit(double vzMin, double vzMax, Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2]);
		int makeInt(double val);
		void plotBcgFit(TH1D*& hist1, TH1D*& hist2, double Min, double Max, TString variable);
		vector<pair<int,double>> effFit(int Switch ,double Min, double Max,Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2], int runSeparatePions);
		void efficiencyFromHists(TH1D* hist1, TH1D* hist2);
		vector<TString> getConditions(int Switch, int idx, double Min, double Max);
		TH1D* getHist(int Switch, double Min, double Max, bool is2TOF, int runSeparatePions, bool runFirstOnly);

		TH1D *hMean1, *hSigma1, *hN1, *hAlpha1;
		TH1D *hMean2, *hSigma2, *hN2, *hAlpha2;

		TH1D *polEta10, *polEta11, *polEta20, *polEta21;
		TH1D *polPhi10, *polPhi11, *polPhi20, *polPhi21;
		TH1D *polPt10, *polPt11, *polPt20, *polPt21;
		TH1D *polVz10, *polVz11, *polVz20, *polVz21;

		TFile *MCAna;

		vector<TH1D*> histsParFit;

		// wider binning
		Double_t effEtaAmpGuess1[9] = {12000, 17000, 9000, 8000, 8000, 8000, 7000, 8000, 8000};
		Double_t effEtaAmpGuess2[9] = {6000, 10000,5500, 6000, 6000, 7000, 6000, 7000, 7000};
		Double_t effEtaPol0Guess1[9] = {37000, 30000,15000, 1000, 4000, 4500, 4000, 5000, 5000};
		Double_t effEtaPol1Guess1[9] = {48000, -14000,10000, 5000, -5000, -4000, -7000, -5000, -7000,};
		Double_t effEtaPol0Guess2[9] = {20000, 18000,46000, 19000, 50000, 1000, 1000, 1000, 1000};
		Double_t effEtaPol1Guess2[9] = {14000, -15000,-48000, -20000, -50000, -1000, -1000, -1000, -1000};
		Double_t effEtaNGuess1[9] = {4, 4,4, 4, 4, 4, 4, 4, 4};
		Double_t effEtaAlphaGuess1[9] = {1.8, 1.8, 1.6, 1.5,1.9, 1.1, 1., 1, 1};
		Double_t effEtaNGuess2[9] = {4, 4,4, 4., 4, 4, 4., 4,4};
		Double_t effEtaAlphaGuess2[9] = {1.1,1,1.1, 1.1, 1.2, 1.2, 1., 1, 1};

		Double_t effPhiAmpGuess1[6] = {30000, 30000, 30000, 30000, 30000, 30000};
		Double_t effPhiAmpGuess2[6] = { 20000, 20000, 20000, 20000,20000, 20000};
		Double_t effPhiPol0Guess1[6] = { 20000, 20000, 20000, 20000, 20000, 20000};
		Double_t effPhiPol1Guess1[6] = {-10000, -10000, -10000, -10000, -10000, -10000};
		Double_t effPhiPol0Guess2[6] = { 1000, 1000, 1000, 1000,1000, 1000};
		Double_t effPhiPol1Guess2[6] = { 5000, 5000, 5000, 5000,5000, 5000};

		Double_t effPtAmpGuess1[6] = {140000, 40000, 40000, 30000, 30000, 30000};
		Double_t effPtAmpGuess2[6] = {80000, 20000, 20000, 5000,5000, 1000};
		Double_t effPtPol0Guess1[6] = { 20000, 90000, 20000, 5000,10000, 10000};
		Double_t effPtPol1Guess1[6] = { -20000, -90000, -5000, 1000,5000, 5000};
		Double_t effPtPol0Guess2[6] = { 30000, 30000, 1000, 1000,1000, 1000};
		Double_t effPtPol1Guess2[6] = { -5000, -5000, 5000, 5000,5000, 5000};

		Double_t effVzAmpGuess1[12] = { 15000, 15000, 15000, 15000, 15000, 15000, 15000, 15000,15000, 15000,15000, 15000};
		Double_t effVzAmpGuess2[12] = {4000, 4000, 4000, 4000, 4000, 4000, 4000, 4000,4000, 4000,4000, 4000};
		Double_t effVzPol0Guess1[12] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,10000, 10000,10000, 10000};
		Double_t effVzPol1Guess1[12] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,10000, 10000,10000, 10000};
		Double_t effVzPol0Guess2[12] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,10000, 10000,10000, 10000};
		Double_t effVzPol1Guess2[12] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,10000, 10000,10000, 10000};

		/*
		Double_t etaMC1[9] = {702, 986, 995, 990, 765, 950, 970, 1017, 748};
		Double_t etaMC1Err[9] = {55, 60, 55, 51, 72, 52, 57, 63, 57};
		Double_t etaMC2[9] = {392, 603, 650, 547, 383, 573, 622, 699, 443};
		Double_t etaMC2Err[9] = {22, 26, 28, 26, 21, 26, 28, 29, 23};

		Double_t VzMC1[6] = {949, 1390, 1642, 1768, 1442, 932};
		Double_t VzMC1Err[6] = {60, 72, 80, 80, 73, 61};
		Double_t VzMC2[6] = {522, 850, 996, 1078, 926, 540};
		Double_t VzMC2Err[6] = {25, 32, 34, 36, 32, 26};

		Double_t pTMC1[5] = {5380, 2099, 506, 130, 23};
		Double_t pTMC1Err[5] = {138, 85, 51, 31, 24};
		Double_t pTMC2[5] = {3133, 1352, 338, 71, 17};
		Double_t pTMC2Err[5] = {60, 40, 21, 10, 6};

		Double_t phiMC1[6] = {1217, 1603, 1060, 981, 1911, 1351};
		Double_t phiMC1[6] = {70, 73, 64, 70, 79, 72};
		Double_t phiMC1[6] = {724, 980, 595, 688, 1160, 765};
		Double_t phiMC1[6] = {29, 34, 27, 29, 37, 30};

		*/
		// default settings for binning and for fits for efficiency dists
		Double_t lowRange = 0.45000;
		Double_t topRange = 0.54000;
		Int_t nBins = 45;
		Double_t outerBorder[2] = {13, 37};
		Double_t innerBorder[2] = {19, 31}; 

		vector<pair<TH1D*,TString>> hists1D;
		vector<pair<TH2F*,TString>> hists2D;

};


#endif