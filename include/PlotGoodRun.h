#ifndef PlotGoodRun_h
#define PlotGoodRun_h


#include "Plot.h"
#include "RunDef.h"
#include "Libraries.h"
#include "Util.h"
#include "ProbRetainEvent.h"

using namespace RooFit;
using namespace std;



class PlotGoodRun : public Plot {
	public:
		//accessible from outside the class
		PlotGoodRun(const string mInputList, const char* filePath);
		PlotGoodRun(const string mInputList, unique_ptr<TFile> &file);
      	~PlotGoodRun(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;

		double getTriggerEfficiency() const { return triggerEfficiency; }
		double getTriggerEfficiencyBeforeCuts() const { return triggerEfficiencyBeforeCuts; }
		vector<int> getGoodRunList() const { return goodRunList; }
		vector<int> getRpGoodRunList() const { return rpGoodRunList; }
		


    private:
    	//store variables, which can be accessed only from within this class
		pair<double,double> FitEtaDistributions(TH1D* hist, int SWITCH);
		Double_t GetMean(TH1* hist);
		void createGoodRunList(vector<int> goodRunList, TString outName);
		vector<int> loadBadEtaPhiRuns(TString nameOfFile);
		//void createAverageTrackHists(vector<pair<array<int, 2>, array<double, 10>>> runNumberMap);

		bool inRangeEtaBemc(double average);
		bool inRangePhiBemc(double average);
		vector<int> getPREList();
		void setBranchAddresses();

		TTree *ZBTree;

		TH1* hGoodRunListFlow;
		int mRunNumber;
		int atLeast1JPsiTrigger, RPsClose, nEventsZBVetoAll, nEventsZBVetoPassed, nEventsTEAll, nEventsTEPassed, nEventsLumiFile, nEventsJPsi;
		double luminosity, luminosityError, instLumi;
		double vetoEfficiency, topologyEfficiency;
		double nTracksBEMC, nClustersBEMC, nTracksTPC, nTracksTOF, nVertices;
		double tpcEtaAverage, bemcEtaAverage, tpcPhiAverage, bemcPhiAverage;
	
	    vector<int> goodRunList, rpGoodRunList;

		//trigger efficiency variables
		TGraph *hTrigEff;
		double triggerEfficiency = 0.0;
		double triggerEfficiencyBeforeCuts= 0.0;


		// average tracks plots
		void createAverageTracksPlots();
		pair<double,double> bemcFitPars, tpcFitPars, bemcPhiFitPars;
};


#endif