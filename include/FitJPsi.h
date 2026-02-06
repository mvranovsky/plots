#ifndef FitJPsi_h
#define FitJPsi_h

#include "Libraries.h"
#include <memory>
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "TVectorD.h"
#include "RooFormulaVar.h"
#include "Config.h"


using namespace std;
using namespace RooFit;

// have to do this, because we are in ROOT 5.34 and setCovarianceMatrix is protected
class RooFitResultPublic : public RooFitResult {
public:
    using RooFitResult::setCovarianceMatrix;


    RooFitResultPublic(const char* name, const char* title)
        : RooFitResult(name, title) {}
};


class FitJPsi {
    public:

        TH1D *hist;
        TString bcg;

        FitJPsi(TH1D *&h, TString b);  // bcg takes arguments "poly1", "poly2", "poly3"
        ~FitJPsi();

        void fitPeak();

        void writeContinuumResult(double left = 0.2, double bottom = 0.35, double right = 0.45, double top = 0.62);
        void writeFitResult();
        void writeSignalResult(bool &fitSignal2S,double left, double bottom, double right, double top);

        RooAddPdf* buildModel(bool &fitBcg, bool &fitSignal, bool &fitSignal2s);

        void setAlphaLoose(bool value) { mAlphaLoose = value; }
        void setNLoose(bool value) { mNLoose = value; }
        bool isAlphaLoose() const { return mAlphaLoose; }
        bool isNLoose() const { return mNLoose; }

        bool isEmbedding() {return mIsEmbedding; }
        void setEmbeddingDescription(TString desc) { mEmbeddingDescription = desc; }
        void setIsEmbedding(bool value) { mIsEmbedding = value; }

        TCanvas* getCanvas() const { c->Update(); return c; }
        void saveCanvas(TFile *&file, TString name = "", TString dir = "");
        void saveCanvas(shared_ptr<TFile> &file, TString name = "", TString dir = "");
        RooFitResult* getFitResult() const { return fitResult; }

        // variables from JPsi fit
        double getMean() const { if(cb || gauss) return cbmean->getVal();  return -999; }
        double getMeanError() const { if(cb || gauss) return cbmean->getError(); return -999; }
        double getSigma() const { if(cb || gauss) return cbsigma->getVal(); return -999; }
        double getSigmaError() const { if(cb || gauss) return cbsigma->getError(); return -999; }
        double getAlpha() const { if(cb) return cbalpha->getVal(); return -999;}
        double getN() const { if(cb) return cbn->getVal(); return -999; }
        double getNSig() const {if(cb || gauss) return nsig->getVal(); return -999; }
        double getNBkg() const {if(bkg) return nbkg->getVal(); return -999; }
        double getLowLimitFit() const { if(cb || gauss) return (cbmean->getVal() -4.0 * cbsigma->getVal()); else return -999; }
        double getHighLimitFit() const { if(cb || gauss) return (cbmean->getVal() + 3.0 * cbsigma->getVal()); else return -999; }
        double getChiSquare() const { return frame->chiSquare("model", "data"); }
        double getChiSquareNDF() const { return (frame->chiSquare("model", "data")/NDF); }

        // variables from JPsi fit integration
        double getYield() const { if(cb || gauss) return netYield; else return -999; }  // extracted yield
        double getErrYield() const { if(cb || gauss) return errNetYield; else return -999; } // error of the extracted yield
        double getYieldBackground() const { if(cb || gauss) return yieldBackground; else return -999; } // yield of background in the range (mu - 4sigma, mu + 3sigma)
        double getErrBackground() const { if(cb || gauss) return errBackground; else return -999; }
        double getSOverB() const { if(cb || gauss) return sOverB; else return -999; }  // ratio of signal yield over background in the range (mu - 4sigma, mu + 3sigma)
        double getErrSOverB() const { if(cb || gauss) return errSOverB; else return -999; }

        // variables from psi(2S) fit
        double getMean2S() const { if(cb2S || gauss2S) return mean2S->getVal();  return -999; }
        double getMean2SError() const { if(cb2S || gauss2S) return mean2S->getError(); return -999; }
        double getSigma2S() const { if(cb2S || gauss2S) return sigma2S->getVal(); return -999; }
        double getSigma2SError() const { if(cb2S || gauss2S) return sigma2S->getError(); return -999; }
        double getNSig2S() const {if(cb2S || gauss2S) return nsig2S->getVal(); return -999; }

        // setters for fit range
        void defineFitRange();
        void setFitRangeLow(double val) { fitRangeLow = val; }
        void setFitRangeHigh(double val) { fitRangeHigh = val; }

        // for fitting additional function on the canvas
        void showAdditionalContinuumFunction(bool value) { showAdditionalContinuum = value; }
        void addContinuumFunction(double c_1, double c_2, double c_3, RooFitResultPublic *&fitResCont);  // for comparison between embedding and data, call only after fitPeak()
        void setVisualizeError(bool value){ visualizeError = value; }
        void saveRooFitResult(TString filename);
        void loadRooFitResult(TString filename);

        // for plotting additional CB peak
        void showDataPeak(double mean = 3.086, double sigma = 0.051) {mShowDataPeak = true; mDataMean = mean; mDataSigma = sigma; }
        void loadDataPeak(double mean, double sigma);


        void setLegendPosition(double x1, double y1, double x2, double y2){
            legPosX1 = x1;
            legPosY1 = y1;
            legPosX2 = x2;
            legPosY2 = y2;
        }

        void printCovarianceMatrix();

    private:

        void integrate(TString bcg);
        void drawLines();

        TCanvas *c;
        RooFitResult* fitResult;
        RooRealVar* x;
        RooDataHist* dh;
        RooPlot* frame;
        TH1D *hSignal;
        int NDF;

        vector<RooAbsPdf*> models;
        vector<TString> modelNames;

        RooRealVar *a0, *a1, *a2;
        RooGenericPdf *bkg;
        
        RooRealVar *cbmean, *cbsigma, *cbalpha, *cbn;
        RooCBShape *cb;
        RooGaussian *gauss;
        RooRealVar *nbkg, *nsig, *nsig2S;
        RooAddPdf *model;

        RooRealVar *c1, *c2, *c3, *c4, *c5;

        //variables to show additional peak from data
        double mDataMean = 3.086;
        double mDataSigma = 0.051;
        bool mShowDataPeak = false;
        RooRealVar *datamean, *datasigma, *dataalpha, *datan;
        RooCBShape *datacb;
        RooAddPdf *datamodel;

        bool mIsEmbedding = false;
        TString mEmbeddingDescription = "";

        RooRealVar *mean2S, *sigma2S, *cbalpha2S, *cbn2S;
        RooCBShape *cb2S;
        RooGaussian *gauss2S;

        TLegend *leg1;

        // variables to define whether to fit crystal ball with set n, alpha or that they are free
        bool mAlphaLoose = false;
        bool mNLoose = false;
        bool visualizeError = false;
        bool showAdditionalContinuum = false;

        double fitRangeLow, fitRangeHigh;


        bool fitBcg = false;
        bool fitSignal = false, fitSignal2S = false;
        
        // integration variables
        double netYield = 0, errNetYield = 0;
        double yieldSignal, yieldBackground;
        double errSignal;
        double errBackground = 0.0;
        double sOverB, errSOverB;

        double legPosX1 = 0.2;
        double legPosY1 = 0.68;
        double legPosX2 = 0.45;
        double legPosY2 = 0.88;

};

#endif