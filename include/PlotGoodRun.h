#ifndef PlotGoodRun_h
#define PlotGoodRun_h


#include "Plot.h"
#include "Config.h"
#include "Libraries.h"
#include "Util.h"
#include "ProbRetainEvent.h"

using namespace RooFit;
using namespace std;



class PlotGoodRun : public Plot {
	public:
		//accessible from outside the class
		PlotGoodRun(const string mInputList, const char* filePath);
		PlotGoodRun(const string mInputList, shared_ptr<TFile> file);
      	~PlotGoodRun(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;

		double getTriggerEfficiency() const { return triggerEfficiency; }
		double getTriggerEfficiencyBeforeCuts() const { return triggerEfficiencyBeforeCuts; }
		vector<int> getGoodRunList() const { return goodRunList; }
		vector<int> getRpGoodRunList() const { return rpGoodRunList; }
		
		void createGoodRunList(vector<int> goodRunList, TString outName);
		vector<int> loadBadEtaPhiRuns(TString nameOfFile);
		void printBadRunsToFile(vector<int> badRuns, TString outName);
		void loadGoodRunList(TString fileName, TString option);
		
		ProbRetainEvent* mProbRetainEvent;
		
	private:
		pair<double,double> FitEtaDistributions(TH1D* hist, int SWITCH);
		Double_t GetMean(TH1* hist);
    	//store variables, which can be accessed only from within this class
		//void createAverageTrackHists(vector<pair<array<int, 2>, array<double, 10>>> runNumberMap);

		bool inRangeEtaBemc(double average);
		bool inRangePhiBemc(double average);
		vector<int> getPREList();
		void setBranchAddresses();

		TTree *ZBTree;

		TH1* hGoodRunListFlow, *hNEventsJPsiFlow;
		int mRunNumber;
		int atLeast1JPsiTrigger, RPsClose, nEventsZBVetoAll, nEventsZBVetoPassed, nEventsTEAll, nEventsTEPassed, nEventsLumiFile;
		double nEventsJPsi;
		double luminosity, luminosityError, instLumi;
		double vetoEfficiency, topologyEfficiency;
		double nTracksBEMC, nClustersBEMC, nTracksTPC, nTracksTOF, nVertices;
		double tpcEtaAverage, bemcEtaAverage, tpcPhiAverage, bemcPhiAverage;
		double nEventsJPsiTrigger0, nEventsJPsiTrigger1, nEventsJPsiTrigger2;
	    vector<int> goodRunList, rpGoodRunList, goodRunTrigger0, goodRunTrigger1, goodRunTrigger2;

		//trigger efficiency variables
		TGraph *hTrigEff;
		double triggerEfficiency = 0.0;
		double triggerEfficiencyBeforeCuts= 0.0;


		// average tracks plots
		void createAverageTracksPlots();
		pair<double,double> bemcFitPars, tpcFitPars, bemcPhiFitPars;
};


#endif