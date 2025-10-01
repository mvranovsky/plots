#ifndef PlotZeroBias_h
#define PlotZeroBias_h


#include "Plot.h"
#include "RunDef.h"
#include "Libraries.h"
#include "Util.h"
#include "ProbRetainEvent.h"
#include "FitJPsi.h"

using namespace RooFit;
using namespace std;



class PlotZeroBias : public Plot {
	public:
		//accessible from outside the class
		PlotZeroBias(const string mInputList, const char* filePath);
		PlotZeroBias(const string mInputList, shared_ptr<TFile> file);
      	~PlotZeroBias(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;

		
		double getTriggerEfficiency() const {return mTriggerEfficiency;}
		double getTriggerEfficiencyError() const {return mTriggerEffError;}
		double getTriggerEfficiencySysError() const {return mTriggerEffSysError; }

		void setNumerator(double val) { mNumerator = val; }
		void setNumeratorError(double val) { mNumeratorError = val; }
		void setDenominator(double val) { mDenominator = val; }
		void setDenominatorError(double val) { mDenominatorError = val; }

		double getNumerator() const { return mNumerator; }
		double getDenominator() const { return mDenominator; }
		double getNumeratorError() const {return mNumeratorError; }
		double getDenominatorError() const {return mDenominatorError; }
		
	private:
		void calculateTriggerEfficiency();
		double triggerEfficiency(TString cond);
		void invMassDependence();
		void invMassDependence2();




		double mNumerator, mDenominator, mNumeratorError, mDenominatorError;
		double mTriggerEffError;
		double mTriggerEfficiency;
		double mTriggerEffSysError;
};


#endif