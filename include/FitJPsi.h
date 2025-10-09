#ifndef FitJPsi_h
#define FitJPsi_h

#include "Libraries.h"
#include <memory>
#include "RooGenericPdf.h"

using namespace std;
using namespace RooFit;


class FitJPsi {
    public:

        TH1D *hist;
        TString bcg;

        FitJPsi(TH1D *&h, TString b);  // bcg takes arguments "poly1", "poly2", "poly3"
        ~FitJPsi();

        void fitPeak();
        void fitContinuum();

        double getCorrectedYield(TH1D *spectrum, TGraphAsymmErrors* graph, double averageCorrection);
        double getCorrectedYieldErrTop() {return correctedYieldErrTop; }
        double getCorrectedYieldErrLow() {return correctedYieldErrLow; }
        void setAlphaLoose(bool value) { mAlphaLoose = value; }
        void setNLoose(bool value) { mNLoose = value; }
        bool isAlphaLoose() const { return mAlphaLoose; }
        bool isNLoose() const { return mNLoose; }

        TCanvas* getCanvas() const { return c; }
        void saveCanvas(TFile *&file, TString name = "", TString dir = "");
        void saveCanvas(shared_ptr<TFile> &file, TString name = "", TString dir = "");
        RooFitResult* getFitResult() const { return fitResult; }
        double getMean() const { return cbmean->getVal(); }
        double getMeanError() const { return cbmean->getError(); }
        double getSigma() const { return cbsigma->getVal(); }
        double getSigmaError() const { return cbsigma->getError(); }
        double getAlpha() const { return cbalpha->getVal(); }
        double getN() const { return cbn->getVal(); }
        double getNSig() const { return nsig->getVal(); }
        double getNBkg() const { return nbkg->getVal(); }
        double getLowLimitFit() const { return (cbmean->getVal() -4.0 * cbsigma->getVal()); }
        double getHighLimitFit() const { return (cbmean->getVal() + 3.0 * cbsigma->getVal()); }
        double getChiSquare() const { return frame->chiSquare("model", "data"); }
        double getChiSquareNDF() const { return (frame->chiSquare("model", "data")/NDF); }

        double getYield() const { return netYield; }
        double getErrYield() const { return errNetYield; }
        double getYieldSignal() const { return yieldSignal; }
        double getYieldBackground() const { return yieldBackground; }
        double getErrSignal() const { return errSignal; }
        double getErrBackground() const { return errBackground; }
        double getSOverB() const { return sOverB; }
        double getErrSOverB() const { return errSOverB; }


        void setContinuumMaximum(double max) { maximumContinuum = max; }
        void setContinuumMinimum(double min) { minimumContinuum = min; }
        double getContinuumMaximum() const { return maximumContinuum; }
        double getContinuumMinimum() const { return minimumContinuum; }
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

        // variables to define whether to fit crystal ball with set n, alpha or that they are free
        bool mAlphaLoose = false;
        bool mNLoose = false;


        bool fitBcg = false;
        // integration variables
        double netYield = 0, errNetYield = 0;
        double yieldSignal, yieldBackground;
        double errSignal;
        double errBackground = 0.0;
        double sOverB, errSOverB;

        // yield correction
        double correctedYield = 0;
        double correctedYieldErrTop = 0;
        double correctedYieldErrLow = 0;

        double minimumContinuum = 1.5;
        double maximumContinuum = 2.5;

};

#endif