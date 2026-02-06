#ifndef PlotEmbeddingJPsi_h
#define PlotEmbeddingJPsi_h

#include "Util.h"
#include "Plot.h"
#include "Config.h"
#include "FitJPsi.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

struct binning {
	int nBins;
	double lowLim;
	double topLim;
};

struct Value {
    double val;
    double errTop;
    double errLow;
};

class PlotEmbeddingJPsi : public Plot {
	public:
		//accessible from outside the class
		PlotEmbeddingJPsi(const string mInputList, const char* filePath);
		PlotEmbeddingJPsi(const string mInputList, shared_ptr<TFile> file);
      	~PlotEmbeddingJPsi(){};

		Value fitPeakJPsi(TH1D* invMass = nullptr);
		void fitContinuum();
	
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;

		void createBemcEfficiencyPlots();
		bool runStudy(int VAR, TString condition);
		void runSysStudy();
		void peakFittingStudy(); // function where i will be changing binning of the invMass histogram and see how it affects the fit results
		void plotContinuum(int nBins, double low, double top);

		double getEfficiency(int VARIABLE, int VARIATION) const {
			return (1.0*yieldResults.at(mUtil->nameOfVariable(VARIABLE))[VARIATION]/mGoodStarlightEntries);
		}

		double getEfficiencyFinal(){ return mEfficiencyFinal; }
		double getEfficiencyErrFinal() { return mEfficiencyErrFinal; }
		TGraphAsymmErrors* reconstructionEfficiency(int SWITCH, TString nameOfOutput = "", bool plotFitFunc = false, bool sumLastBins = false);  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters

		TGraphAsymmErrors* plotEfficiency(vector<TH1D*> h1,vector<TH1D*> h2, TString nameOfOutput, vector<double> range , TString xAxisDescription,TString yAxisDescription = "reconstruction efficiency", TString dir = "recoEffPlots" , bool plotFitFunc = false, bool sumLastBins = false); // plots efficiency graph from two histograms{
		TGraphAsymmErrors* plotEfficiency(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TString nameOfOutput, vector<double> range, TString xAxisDescription,TString yAxisDescription,  TString dir = "recoEffPlots" , bool plotFitFunc = false, bool sumLastBins = false);
		TH1D* scanDir(TDirectory *&dir, TString nameOfHist);
		bool isDir(TString dirPath);

		void changeBinning(int SWITCH, int nBins, double low, double top);
		bool isJPsiEmbedding(TH1D* h);
		bool isRunSysStudy();

    private:
    	//store variables, which can be accessed only from within this class
		void controlPlotsComparison(bool justJPsi);
		bool plot2Dists(TH1 *hData, TH1* hEmb, TString outName);


		Double_t lowLimInvMass = 0.0, topLimInvMass = 7.0;

		TString bemcEffDir = "BemcEffPlots";



		Util *mUtil;


		double peakFittingMean = 0.0;
		double peakFittingSigma = 0.0;
		double peakFittingMeanError = 0.0;
		double peakFittingSigmaError = 0.0;

		double mEfficiencyFinal;
		double mEfficiencyErrFinal;
		int mGoodStarlightEntries;

		map<TString, vector<int>> yieldResults;
		int nBins = 40; // number of bins for invariant mass histogram
		double lowerLim = 2.5; // lower limit for invariant mass histogram
		double upperLim = 3.5; // upper limit for invariant mass histogram

		vector<binning> recoBins = {
			{40, 2.5, 3.5}, // inv mass fitting
			{20, -1, 1},  //pair rapidity
			{40, -1, 1},  //eta of daughters
			{40, -3.14, 3.14}, // phi of daughters
			{30, 0, 1.5}, // pT of Jpsi
			{40, 0.5, 2.5}  // pT of daughters
		};

};


#endif