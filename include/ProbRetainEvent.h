#ifndef ProbRetainEvent_h
#define ProbRetainEvent_h


#include "Libraries.h"
#include "Config.h"

using namespace std;


struct runInfo {
    int runNumber;
    double luminosity;
    double probRetainEvent;
    double instLumi;
    int nEventsZBVetoAll;
    int nEventsZBVetoPassed;
    int nEventsZBTofMultPassed;
    int nEventsZBBbcEPassed;
    int nEventsZBBbcWPassed;
    double nEventsJPsi;
    int nVerticesZB0;
    int nVerticesZB1;
    int nVerticesZB2More;
};


class ProbRetainEvent {
    
    public:
        ProbRetainEvent(const std::shared_ptr<TFile>& outFile, TTree* tree);
        ~ProbRetainEvent();

        void Make(TString fitOption);  

        void defineFitOption(TString option);

        vector<runInfo> getData();
        double getLuminosity(const vector<int>& runs, double A, double B, double pileUpA, double pileUpB);
        vector<int> filterProbRetainEvent(const vector<runInfo>& data, double mean, double sigma, double A, double B);
        vector<int> getGoodRunList() const { return goodRuns; }
		void DrawSTARInternal(double xl = 0.8, double yl = 0.89, double xr = 0.93, double yr = 0.93);

        void pileUpCorrection(vector<runInfo> &data);

        double getA() const { return a; }
        double getB() const { return b; }
        double getAError() const { return aErr; }
        double getBError() const { return bErr; }

        void setPileUpCorrectionA(double val) { pileUpA = val; }
        void setPileUpCorrectionB(double val) { pileUpB = val; }
        void setPileUpCorrectionAError(double val) { pileUpAErr = val; }
        void setPileUpCorrectionBError(double val) { pileUpBErr = val; }

        double getPileUpCorrectionA() const { return pileUpA; }
        double getPileUpCorrectionB() const { return pileUpB; }
        double getPileUpCorrectionAError() const { return pileUpAErr; }
        double getPileUpCorrectionBError() const { return pileUpBErr; }

        void plotTofMultVeto(TString option);
        void plotBbcEVeto(TString option);
        void plotBbcWVeto(TString option);      
        void runCustomList(vector<int> runList, TString option);  
    
        private:
        void readLumiFile(bool isZB );
        TString dir;
   		void saveTopologyEfficiency(vector<runInfo> &dat);
        TGraphErrors* fillGraph(vector<runInfo> &data);
        vector<double> fitGaussian(TH1D *h, vector<runInfo> data);
        double exponential(double x, double A, double B);  // function to calculate the exponential
        double linear(double x, double A, double B);
        vector<runInfo> data;  // vector to store run information
        shared_ptr<TFile> outFile;
        TTree* tree;    
        vector<int> goodRuns;  // vector to store good runs

        TString usedFitOption;
        double a = 0;
        double b = 0;
        double aErr = 0;
        double bErr = 0;

        double pileUpA = 0;
        double pileUpB = 0;
        double pileUpAErr = 0;
        double pileUpBErr = 0;

        bool tofMultVeto = false, bbcEVeto = false, bbcWVeto = false;

        map<int,double> mInstLumiPerRun;
        map<int,int> mNEventsLumiFile;
        map<int,double> mLumiPerRun;

        vector<int> customList;
        bool runningCustomList = false;

        // fit option variables
        TString fitFunc, fitDescription, fitFuncLegend;
};


#endif
