#ifndef PlotEmbeddingJPsi_h
#define PlotEmbeddingJPsi_h

#include "Util.h"
#include "Plot.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"
#include "FitJPsi.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class PlotEmbeddingJPsi : public Plot {
	public:
		//accessible from outside the class
		PlotEmbeddingJPsi(const string mInputList, const char* filePath);
		PlotEmbeddingJPsi(const string mInputList, unique_ptr<TFile> &file);
      	~PlotEmbeddingJPsi(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


		bool runStudy(int VAR, TString condition);
		void runSysStudy();
		void peakFittingStudy(); // function where i will be changing binning of the invMass histogram and see how it affects the fit results


		double getEfficiency(int VARIABLE, int VARIATION) const {
			return (yieldResults.at(mUtil->nameOfVariable(VARIABLE))[VARIATION]/starlightTree->GetEntries());
		}

		double getEfficiencyFinal(){ return mEfficiencyFinal; }
		double getEfficiencyErrFinal() { return mEfficiencyErrFinal; }


    private:
    	//store variables, which can be accessed only from within this class
		void reconstructionEfficiency(int SWITCH, TString nameOfOutput);  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
		void controlPlotsComparison(bool justJPsi);
		bool plot2Dists(TH1 *hData, TH1* hEmb, TString outName);


		void  saveSysStudyYieldsHists();
		Double_t lowLimInvMass, topLimInvMass;


		unique_ptr<TFile> starlightFile;
		TTree *starlightTree;

		Util *mUtil;
		// variables for sys study
		TH1 *hSysStudyLoose, *hSysStudyTight;

		double peakFittingMean = 0.0;
		double peakFittingSigma = 0.0;
		double peakFittingMeanError = 0.0;
		double peakFittingSigmaError = 0.0;

		double mEfficiencyFinal;
		double mEfficiencyErrFinal;

		map<TString, vector<int>> yieldResults;
		int nBins = 40; // number of bins for invariant mass histogram
		double lowerLim = 2.5; // lower limit for invariant mass histogram
		double upperLim = 3.5; // upper limit for invariant mass histogram
};


#endif