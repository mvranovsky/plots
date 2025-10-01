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
#include "PlotZeroBias.h"


using namespace std;
using namespace UTIL;

class CrossSectionMaker  {

    public:
        // Constructor
        CrossSectionMaker(const TString AnaDir, const TString EmbedddingDir, const TString GoodRunDir, const TString ZeroBiasDir, TString outputFile);
        ~CrossSectionMaker();

        void Make();

        // Cross section functions
        double CalculateCrossSection();
        double getYieldCorrection(int spectrumBins, double spectrumLow, double spectrumTop, TString condition);
        void setNominalCrossSection(double val) { mNominalCrossSection = val; }
        double getNominalCrossSection() const { return mNominalCrossSection; }

        // plot of cross section dependent on rapidity
        void rapidityDependencePlot(int nBins, double low, double top,vector<double> plotEdges);
    
        // hard coded for now
        double getRPEfficiency() const { return isRomanPots ? mRPEfficiency : 1.0; }
        double getRPEffSysError() const { return isRomanPots ? 0.01 : 0.0; } // hard coded for now, we will change it later
        double getVertexRecoEfficiency() const { return 0.9; } // hard coded for now, we will change it later
        double getVertexRecoEffSysError() const { return 0.01; }

        double crossSection(double Y,double eff, double lumi, double deltaRap = 2);
        double efficiency(vector<double> effs);
        double luminosity();
        TGraph* getMcNormalizedPlot(int nBins, double low, double top);

        // yield correction saving
        void setYieldCorrectionCanvas(TCanvas* can) { yieldCorrectionCanvas = can; }
        void setYieldCorrectionSpectrum(TH1D* h) { hYieldCorrection = h; }
        void setYieldCorrectionGraph(TGraphAsymmErrors* g) { gYieldCorrection = g; }

        TH1D* getYieldCorrectionSpectrum() const { return hYieldCorrection; }
        TGraphAsymmErrors* getYieldCorrectionGraph() const { return gYieldCorrection; }
        TCanvas* getYieldCorrectionCanvas() const { return yieldCorrectionCanvas; }

        void saveYieldCorrection(TString dir);

        void setYieldCorrFinal(double val){ mCorrYieldFinal = val;}
        void setYieldCorrErrTopFinal(double val){ mCorrYieldErrTopFinal = val;}
        void setYieldCorrErrLowFinal(double val){ mCorrYieldErrLowFinal = val;}

        double getYieldCorrFinal() const {return mCorrYieldFinal;}
        double getYieldCorrErrTopFinal() const {return mCorrYieldErrTopFinal;}
        double getYieldCorrErrLowFinal() const {return mCorrYieldErrLowFinal;}
        
        
        
    private:
        
        // computation and summing of all systematic errors
        pair<double,double> SysError();
        pair<double,double> CalculateSysErrorCuts(int VAR);
        pair<double,double> CalculateRPEffSysError();
        pair<double,double> CalculateVertexRecoSysError();
        pair<double,double> CalculateYieldExtractionSysError();
        pair<double,double> CalculateBemcSysError();
        pair<double,double> CalculateLumiSysError();
        pair<double,double> CalculateTriggerVetoSysError();
        pair<double,double> CalculateTriggerTopologySysError();






        //vector<double> corrrectedYield();
        shared_ptr<TFile> outFile;

        Util *mUtil;
        ProbRetainEvent *mProbRetainEvent;
        PlotEmbeddingJPsi *mEmbedding;
        PlotAnaJPsi *mAna;
        PlotGoodRun *mGoodRun;
        PlotZeroBias *mZeroBias;

        TString mAnaDir, mEmbeddingDir, mGoodRunDir, mZeroBiasDir, mOutFileName;

        shared_ptr<TFile> mOutFile;

        TTree *mGoodRunTree = nullptr;  // tree to store good runs

        double mRPEfficiency = 0.85; // so far hard coded, we will change it later

        bool isRomanPots = false;

        map<int, double> mSysErrorTight;
        map<int, double> mSysErrorLoose;

        // yield correction variables
        double mCorrectedYield = 0, mCorrYieldErrTop = 0, mCorrYieldErrLow = 0, mRawYield = 0;
        double mNominalCrossSection;

        double mCorrYieldFinal = 0, mCorrYieldErrTopFinal = 0, mCorrYieldErrLowFinal = 0;


        TCanvas *yieldCorrectionCanvas = nullptr;
        TH1D* hYieldCorrection = nullptr;
        TGraphAsymmErrors* gYieldCorrection = nullptr;

};

#endif //CrossSectionMaker_h

