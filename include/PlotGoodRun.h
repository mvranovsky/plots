#ifndef PlotGoodRun_h
#define PlotGoodRun_h


#include "Plot.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"
#include "Libraries.h"
#include "Util.h"

using namespace RooFit;
using namespace std;

class PlotGoodRun : public Plot {
	public:
		//accessible from outside the class
		PlotGoodRun(const string mInputList, const char* filePath);
      	~PlotGoodRun(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class
		pair<double,double> FitEtaDistributions(TH1D* hist, int SWITCH);
		Double_t GetMean(TH1* hist);
		void AnalysisBEMC(vector<int> runNumberList);
		void createGoodRunList(vector<int> goodRunList);
		vector<int> loadProbRetainEventList();
		vector<int> loadBadEtaPhiRuns(TString nameOfFile);

		void createAverageTrackHists(vector<pair<array<int, 2>, array<double, 10>>> runNumberMap);


		Int_t mIsRpOk , mIsJPsiTrigger1, mIsJPsiTrigger2, mIsJPsiTrigger3, nTracksBemc, nTracksTof, nClustersBemc, runNumber;
		vector<Double_t> mTpcTrack_phi, mTpcTrack_eta, mBemcTrack_eta, mBemcTrack_phi, mTpcNHitsFit, mTpcNHitsDEdx, mTpcNSigmaElectron, mTpcNSigmaPion, mTpcNSigmaProton, mTpcNSigmaKaon;
	    // set branch addresses to these variables:
		UInt_t mFillNumber,mEventNumber, mBunchCrossId, mBunchCrossId7bit;
		UShort_t mTofMult;
		UInt_t mNVertecies, mNGoodTpcTrks;
		TH1* hGoodRunListFlow;



		int mRunNumber, atLeast1JPsiTrigger, RPsClose, nEventsAll, nEventsPassed;
		double luminosity, luminosityError, nTracksBEMC, nClustersBEMC, nTracksTPC, nTracksTOF, nVertices;
		double tpcEtaAverage, bemcEtaAverage, tpcPhiAverage, bemcPhiAverage;
	


};


#endif