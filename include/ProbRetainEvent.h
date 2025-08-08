#ifndef ProbRetainEvent_h
#define ProbRetainEvent_h


#include "Libraries.h"


using namespace std;


struct runInfo {
    int runNumber;
    double luminosity;
    double probRetainEvent;
    double instLumi;
    int nEventsZBVetoAll;
    int nEventsZBVetoPassed;
    int TEAll;
    int TEPassed;
};


class ProbRetainEvent {
    
    public:
        ProbRetainEvent(const std::unique_ptr<TFile>& outFile, TTree* tree);
        ~ProbRetainEvent();

        void Make();  

        vector<runInfo> getData();
        double calculateSysError(vector<int> runs = {});
        double getLuminosity(const vector<int>& runs, double A, double B);
        vector<int> filterProbRetainEvent(const vector<runInfo>& data, double mean, double sigma, double A, double B);
        vector<int> getGoodRunList() const { return goodRuns; }

        double getA() const { return a; }
        double getB() const { return b; }
        double getAError() const { return aErr; }
        double getBError() const { return bErr; }
        
        // popridavat funkcie pre ziskanie luminozity pre dany run, korekcia na luminozitu... 
        
    private:

        TString dir;
   		void saveTopologyEfficiency(vector<runInfo> &dat);
        TGraph* fillGraph(vector<runInfo> data);
        vector<double> fitGaussian(TH1D *h, vector<runInfo> data);
        double exponential(double x, double A, double B);  // function to calculate the exponential
        map<int, double> RunNumToLumi;  // map to store run number and luminosity
        map<int, double> RunNumToInstLumi;  // map to store run number and instantaneous luminosity
        vector<runInfo> data;  // vector to store run information
        TFile* outFile;
        TTree* tree;    
        vector<int> goodRuns;  // vector to store good runs


        // values from the exponential fit
        double a = 0;
        double b = 0;
        double aErr = 0;
        double bErr = 0;

};


#endif
