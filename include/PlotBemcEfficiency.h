#ifndef PlotBemcEfficiency_h
#define PlotBemcEfficiency_h


#include "Plot.h"
#include "RunDef.h"
#include "Libraries.h"
#include "Util.h"

using namespace RooFit;
using namespace std;



class PlotBemcEfficiency : public Plot {
	public:
		//accessible from outside the class
		PlotBemcEfficiency(const string mInputList, const char* filePath);
		PlotBemcEfficiency(const string mInputList, shared_ptr<TFile> file);
      	~PlotBemcEfficiency(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;

		inline void setEfficiencyGraph(TGraphAsymmErrors* g) { mEfficiencyGraph = g; }

		inline TGraphAsymmErrors* getEfficiencyGraph() const { return mEfficiencyGraph; }	

		void setEpsilon(double val) { mEpsilon = val; }
		void setEpsilonError(double val) { mEpsilonError = val; }
		void setN(double val) { mN = val; }
		void setNError(double val) { mNError = val; }
		void setPtThreshold(double val) { mPtThreshold = val; }
		void setPtThresholdError(double val) { mPtThresholdError = val; }
		void setSigma(double val) { mSigma = val; }
		void setSigmaError(double val) { mSigmaError = val; }

		double getEpsilon() const { return mEpsilon; }
		double getEpsilonError() const { return mEpsilonError; }
		double getN() const { return mN; }
		double getNError() const { return mNError; }
		double getPtThreshold() const { return mPtThreshold; }
		double getPtThresholdError() const { return mPtThresholdError; }
		double getSigma() const { return mSigma; }
		double getSigmaError() const { return mSigmaError; }

	private:

		TGraphAsymmErrors* bemcEfficiencyPt(int nBins, double low, double top);

		double mEpsilon, mEpsilonError, mN, mNError, mPtThreshold, mPtThresholdError, mSigma, mSigmaError;
		TGraphAsymmErrors* mEfficiencyGraph;
};


#endif