#ifndef PlotAnaJPsi_h
#define PlotAnaJPsi_h


#include "Util.h"
#include "Plot.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"
#include "FitJPsi.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class PlotAnaJPsi : public Plot {
	public:
		//accessible from outside the class
		PlotAnaJPsi(const string mInputList, const char* filePath);
		PlotAnaJPsi(const string mInputList, unique_ptr<TFile> &file);
      	~PlotAnaJPsi(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;

		// methods about Vz study
		void vertexZStudy(); // function to calculate systematic error of Vz position
		double getVzStudyFitEff() const { return VzStudyFitEff; }
		double getVzStudyHistEff() const { return VzStudyHistEff; }
		double getVzStudyAverageEff() const { return (VzStudyFitEff + VzStudyHistEff) / 2.0; }
		double getVzStudySysError() const { return abs(VzStudyFitEff - VzStudyHistEff)/2.0; }

		// methods about peak fitting study
		double getPeakFittingMean() const { return peakFittingMean; }
		double getPeakFittingSigma() const { return peakFittingSigma; }
		double getPeakFittingMeanError() const { return peakFittingMeanError; }
		double getPeakFittingSigmaError() const { return peakFittingSigmaError; }


		bool runStudy(int VAR, TString condition);
		void runSysStudy();
		void peakFittingStudy(); // function where i will be changing binning of the invMass histogram and see how it affects the fit results

		double getYieldFinal(){ return mYieldFinal; }
		double getYieldErrFinal(){ return mYieldErrFinal; }

		double getYield(int VARIABLE, int VARIATION) const {return yieldResults.at(mUtil->nameOfVariable(VARIABLE))[VARIATION];}

		bool isRomanPots(){return !noRomanPots; }
    private:
		// control plots for comparison between data and embedding. This function just saves data histograms to a TFile
		void controlPlotsComparison(bool justJPsi);
		void saveSysStudyYieldsHists(); // saves the yields of systematic study to histograms

		Double_t lowLimInvMass, topLimInvMass;
		FitJPsi* fit;

		Util *mUtil;

		// results of fitting
		double mYieldFinal, mYieldErrFinal;

		// graphs for vertex z study
		TGraph *VzStudyMeanGraph, *VzStudySigmaGraph, *VzStudyEfficiencyFitGraph, *VzStudyEfficiencyHistGraph;
		double VzStudyFitEff = 0.0;
		double VzStudyHistEff = 0.0;
		

		// variables for sys study
		TH1D *hSysStudyLoose, *hSysStudyTight;

		double peakFittingMean = 0.0;
		double peakFittingSigma = 0.0;
		double peakFittingMeanError = 0.0;
		double peakFittingSigmaError = 0.0;


		map<TString, vector<int>> yieldResults;
		int nBins = 40; // number of bins for invariant mass histogram
		double lowerLim = 2.5; // lower limit for invariant mass histogram
		double upperLim = 3.5; // upper limit for invariant mass histogram
		
};


#endif