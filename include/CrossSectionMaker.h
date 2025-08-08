#ifndef CrossSectionMaker_h
#define CrossSectionMaker_h

#include "Plot.h"
#include "Libraries.h"
#include "Util.h"
#include "ProbRetainEvent.h"
#include "PlotEmbeddingJPsi.h"
#include <iostream>
#include <vector>
#include "PlotAnaJPsi.h"
#include "PlotGoodRun.h"


using namespace std;
using namespace UTIL;

class CrossSectionMaker  {

    public:
        // Constructor
        CrossSectionMaker(const TString AnaDir, const TString EmbedddingDir, const TString GoodRunDir, TString outputFile);
        ~CrossSectionMaker();

        void Make();

        double CalculateCrossSection();

        void CalculateSysError(int VAR);

        double getRPEfficiency() const { return isRomanPots ? mRPEfficiency : 1.0; }
        double getTopologyEfficiency() const { return 0.42; } // hard coded for now, we will change it later
    private:

        //vector<double> corrrectedYield();
        unique_ptr<TFile> outFile;

        Util *mUtil;
        ProbRetainEvent *mProbRetainEvent;
        PlotEmbeddingJPsi *mEmbedding;
        PlotAnaJPsi *mAna;
        PlotGoodRun *mGoodRun;

        TString mAnaDir, mEmbeddingDir, mGoodRunDir, mOutFileName;
        
        unique_ptr<TFile> mOutFile;

        TTree *mGoodRunTree = nullptr;  // tree to store good runs

        double mRPEfficiency = 0.85; // so far hard coded, we will change it later

        bool isRomanPots;

        map<int, double> mSysErrorTight;
        map<int, double> mSysErrorLoose;


};

#endif //CrossSectionMaker_h

