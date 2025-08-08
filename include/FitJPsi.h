#ifndef FitJPsi_h
#define FitJPsi_h

#include "Libraries.h"
#include <memory>

using namespace std;
using namespace RooFit;


class FitJPsi {
    public:

        TH1D *hist;
        TString bcg;

        FitJPsi(TH1D *&h, TString b);  // bcg takes arguments "poly1", "poly2", "poly3"
        ~FitJPsi();

        void fitPeak();


        TCanvas* getCanvas() const { return c; }
        void saveCanvas(TFile *&file, TString name = "", TString dir = "");
        void saveCanvas(unique_ptr<TFile> &file, TString name = "", TString dir = "");
        RooFitResult* getFitResult() const { return fitResult; }
        double getMean() const { return cbmean->getVal(); }
        double getMeanError() const { return cbmean->getError(); }
        double getSigma() const { return cbsigma->getVal(); }
        double getSigmaError() const { return cbsigma->getError(); }
        double getAlpha() const { return cbalpha->getVal(); }
        double getN() const { return cbn->getVal(); }
        double getNSig() const { return nsig->getVal(); }
        double getNBkg() const { return nbkg->getVal(); }
        double getLowLimitFit() const { return cbmean->getVal() -4.0 * cbsigma->getVal(); }
        double getHighLimitFit() const { return cbmean->getVal() + 3.0 * cbsigma->getVal(); }
        double getChiSquare() const { return frame->chiSquare("model", "data"); }
        double getChiSquareNDF() const { return frame->chiSquare("model", "data")/NDF; }

        double getYield() const { return netYield; }
        double getErrYield() const { return errNetYield; }
        double getYieldSignal() const { return yieldSignal; }
        double getYieldBackground() const { return yieldBackground; }
        double getErrSignal() const { return errSignal; }
        double getErrBackground() const { return errBackground; }
        double getSOverB() const { return sOverB; }
        double getErrSOverB() const { return errSOverB; }

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

        RooRealVar *a0, *a1, *a2;
        RooPolynomial *bkg;

        RooRealVar *cbmean, *cbsigma, *cbalpha, *cbn;
        RooCBShape *cb;
        RooRealVar *nbkg, *nsig;
        RooAddPdf *model;


        bool fitBcg = false;
        // integration variables
        double netYield, errNetYield;
        double yieldSignal, yieldBackground;
        double errSignal;
        double errBackground = 0.0;
        double sOverB, errSOverB;

};

#endif