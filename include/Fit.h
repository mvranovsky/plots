#ifndef Fit_h
#define Fit_h

#include "Libraries.h"
#include "Config.h"

using namespace std;
using namespace RooFit;


class Fit {
    public:

        TH1D *hist;
        TString bcg;

        Fit(TH1D *&h, TString b);  // bcg takes arguments "poly1", "poly2", "poly3"
        ~Fit();

        void fitPeak();

        // methods to write results on the canvas
        void writeFitResult();
        void writeSignalResult(double left, double bottom, double right, double top);
        
        // methods for handling canvas
        TCanvas* getCanvas() const { c->Update(); return c; }
        void saveCanvas(TFile *&file, TString name = "", TString dir = "");
        void saveCanvas(shared_ptr<TFile> &file, TString name = "", TString dir = "");
        
        // variables from peak fit
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
        
        // variables from peak fit integration
        double getYield() const { if(cb || gauss) return netYield; else return -999; }  // extracted yield
        double getErrYield() const { if(cb || gauss) return errNetYield; else return -999; } // error of the extracted yield
        double getYieldBackground() const { if(cb || gauss) return yieldBackground; else return -999; } // yield of background in the range (mu - 4sigma, mu + 3sigma)
        double getErrBackground() const { if(cb || gauss) return errBackground; else return -999; }
        double getSOverB() const { if(cb || gauss) return sOverB; else return -999; }  // ratio of signal yield over background in the range (mu - 4sigma, mu + 3sigma)
        double getErrSOverB() const { if(cb || gauss) return errSOverB; else return -999; }
        
        // methods for crystal ball fit parameters n, alpha
        void setAlphaLoose(bool value) { mAlphaLoose = value; }
        void setNLoose(bool value) { mNLoose = value; }
        bool isAlphaLoose() const { return mAlphaLoose; }
        bool isNLoose() const { return mNLoose; }

        //initial guesses for fit parameters
        void setPolynomial(double c0, double c1 = 0, double c2 = 0) {mPolynomialPars[0] = c0; mPolynomialPars[1] = c1; mPolynomialPars[2] = c2;}
        void setMean(double mean) { mMeanPar = mean; }
        void setSigma(double sigma) { mSigmaPar = sigma; }

        // setters for fit range
        void defineFitRange();
        void setFitRangeLow(double val) { fitRangeLow = val; }
        void setFitRangeHigh(double val) { fitRangeHigh = val; }
        

        
        void setLegendPosition(double x1, double y1, double x2, double y2){
            legPosX1 = x1;
            legPosY1 = y1;
            legPosX2 = x2;
            legPosY2 = y2;
        }
        
        
    private:
        // method which decides what type of functions will be used
        RooAddPdf* buildModel();
        
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
        RooPolynomial *bkg;
        
        RooRealVar *cbmean, *cbsigma, *cbalpha, *cbn;
        RooCBShape *cb;
        RooGaussian *gauss;
        RooRealVar *nbkg, *nsig, *nsig2S;
        RooAddPdf *model;

        RooRealVar *c1, *c2, *c3, *c4, *c5;

        TLegend *leg1;

        // variables to define whether to fit crystal ball with set n, alpha or that they are free
        bool mAlphaLoose = false;
        bool mNLoose = false;
        bool visualizeError = false;
        Double_t mPolynomialPars[3] = {-2, 2, 1};
        Double_t mMeanPar = 3.1;
        Double_t mSigmaPar = 0.05;

        double fitRangeLow, fitRangeHigh;

        bool fitBcg = false;
        bool fitSignal = false;
        
        // integration variables
        double netYield = 0, errNetYield = 0;
        double yieldSignal, yieldBackground;
        double errSignal;
        double errBackground = 0.0;
        double sOverB, errSOverB;

        // position of legend
        double legPosX1 = 0.2;
        double legPosY1 = 0.68;
        double legPosX2 = 0.45;
        double legPosY2 = 0.88;

};

#endif