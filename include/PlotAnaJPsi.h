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
		PlotAnaJPsi(const string mInputList, shared_ptr<TFile> file);
      	~PlotAnaJPsi(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;

		// methods about Vz study
		void vertexZStudy(); // function to calculate systematic error of Vz position
		double getVzStudyFitEff() const { return VzStudyFitEff; }
		double getVzStudyHistEff() const { return VzStudyHistEff; }
		double getVzStudyAverageEff() const { return (VzStudyFitEff + VzStudyHistEff) / 2.0; }
		double getVzStudySysError() const { return abs(VzStudyFitEff - VzStudyHistEff)/2.0; }

		// methods about peak fitting study
		double getPeakFittingMean() const { return mPeakFittingMean; }
		double getPeakFittingSigma() const { return mPeakFittingSigma; }
		double getPeakFittingMeanError() const { return mPeakFittingMeanError; }
		double getPeakFittingSigmaError() const { return mPeakFittingSigmaError; }

		void setPeakFittingMean(double val) { mPeakFittingMean = val; }
		void setPeakFittingSigma(double val) { mPeakFittingSigma = val; }
		void setPeakFittingMeanError(double val) { mPeakFittingMeanError = val; }
		void setPeakFittingSigmaError(double val) { mPeakFittingSigmaError = val; }

		void drawAndFitHist(TH1D *h, TString outName, double min, double max);
		bool runStudy(int VAR, TString condition);
		void runSysStudy();
		void peakFittingStudy(); // function where i will be changing binning of the invMass histogram and see how it affects the fit results

		double getYieldFinal(){ return mYieldFinal; }
		double getYieldErrFinal(){ return mYieldErrFinal; }
		void setYieldFinal(double val) {mYieldFinal = val; }
		void setYieldErrFinal(double val) {mYieldErrFinal = val; }

		double getYield(int VARIABLE, int VARIATION) const {return yieldResults.at(mUtil->nameOfVariable(VARIABLE))[VARIATION];}


		double getRapidityYield(int bin) const { return mRapidityYields[bin]; }
		double getRapidityYieldError(int bin) const { return mRapidityYieldsErrors[bin]; }
		pair<double,double> getRapidityBin(int bin) const { return mRapidityBins[bin]; }
		int getNRapidityBins() const { return mRapidityBins.size(); }

		bool isRomanPots(){return (!noRomanPots); }

    private:
		// control plots for comparison between data and embedding. This function just saves data histograms to a TFile
		void controlPlotsComparison(bool justJPsi);
		void saveSysStudyYieldsHists(); // saves the yields of systematic study to histograms
		void plotRapidityDependence(int nRapBins);

		void plotContinuum(int nBins, double low, double top);

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

		double mPeakFittingMean = 0.0;
		double mPeakFittingSigma = 0.0;
		double mPeakFittingMeanError = 0.0;
		double mPeakFittingSigmaError = 0.0;


		map<TString, vector<int>> yieldResults;
		int nBins = 40; // number of bins for invariant mass histogram
		double lowerLim = 2.5; // lower limit for invariant mass histogram
		double upperLim = 3.5; // upper limit for invariant mass histogram


		vector<double> mRapidityYields, mRapidityYieldsErrors;
		vector<pair<double,double>> mRapidityBins;
		
};


#endif