#ifndef SysStudy_h
#define SysStudy_h

#include "Util.h"
#include "Plot.h"
#include "RunDef.h"
#include "ProbRetainEvent.h"
#include "FitJPsi.h"


using namespace RooFit;
using namespace std;
using namespace UTIL;

class SysStudy : public Plot {
	public:
		//accessible from outside the class
		SysStudy(const string mInputList, const char* filePath);
		SysStudy(const string mInputList, unique_ptr<TFile> &file);
      	~SysStudy(){};
      
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

		double getYield(int VARIABLE, int VARIATION) const {
			return yieldResults.at(varNames[VARIABLE])[VARIATION];
		}


		void setEmbeddingTrue() { isEmbedding = true; }
    private:
    	//store variables, which can be accessed only from within this class

		double GetAverageYValue(TH1 *hist);
		int studyDiscriminator(TString canvasName);
		void saveSysStudyYieldsHists(); // saves the yields of systematic study to histograms

		vector<pair<double,double>> yRange = { {0.9, 1.1} , {0.8, 1.2}, {0.9, 1.1}, {0.8, 1.2} }; //nHitsFit, nHitsDedx, dcaZInCm, dcaXYInCm
		bool isEmbedding = false;


		// graphs for vertex z study
		TGraph *VzStudyMeanGraph, *VzStudySigmaGraph, *VzStudyEfficiencyFitGraph, *VzStudyEfficiencyHistGraph;
		double VzStudyFitEff = 0.0;
		double VzStudyHistEff = 0.0;
		

		// variables for sys study
		TH1D *hSysStudyLoose, *hSysStudyTight;

		enum {
			NHITSFIT = 0,
			NHITSDEDX,
			DCAZINCM,
			DCAXYINCM,
			PID,
			VERTEXZ,
			ETA,
			nVariables
		};

		TString nameOfHist[nVariables] = {
			"N^{Fit}_{hits}",
			"N^{Dedx}_{hits}",
			"DCA_{Z} [cm]",
			"DCA_{XY} [cm]",
			"PID",
			"V_{Z} [cm]",
			"#eta"
		};

		TString varNames[nVariables] = {
			"nHitsFit",
			"nHitsDedx",
			"dcaZInCm",
			"dcaXYInCm",
			"PID",
			"vertexZInCm",
			"etaHadron"
		};

		

		enum {NOMINAL = 0,TIGHT, LOOSE, nVariations};

		TString variations[nVariations] = {
			"nominal",
			"tight",
			"loose"
		};



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