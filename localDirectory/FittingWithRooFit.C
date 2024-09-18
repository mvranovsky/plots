#include <iostream>
#include <cmath>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>


using namespace RooFit;
using namespace std;


void fit(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);
void fitLambda(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);
double integrateGaussian(RooRealVar& x, RooGenericPdf& gauss, RooRealVar& mean, RooRealVar& sigma, double binWidth);
double GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data);
void fitWith2Gauss(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);
void crystalBall(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);
void crystalBallLambda(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);


void FittingWithRooFit(const string input) {
//const TString& inputTag, const TString& histK0sName, const TString& histLambdaName
	TString inputTag = input;
	TString histK0sName1 = "invMassK0sFit1";
	TString histLambdaName1 = "invMassLambdaFit1"; //change accordingly
    TString histK0sName2 = "invMassK0sFit2";
    TString histLambdaName2 = "invMassLambdaFit2"; 
	TString inputFile = inputTag + TString("/AnalysisOutput.root");

   	TFile* file = new TFile(inputFile);
   	if(!file){
   		cerr << "Couldn't open one of the input files. Leaving." << endl;
   		return;
   	}
   	cout << "Opened input file." << endl;

   	//load the correct histograms
   	TH1D* histK0s1TOF = (TH1D*)file->Get(histK0sName1);
   	TH1D* histLambda1TOF = (TH1D*)file->Get(histLambdaName1);

    TH1D* histK0s2TOF = (TH1D*)file->Get(histK0sName2);
    TH1D* histLambda2TOF = (TH1D*)file->Get(histLambdaName2);

   	if(!histK0s1TOF || !histK0s2TOF ){
   		cerr << "Could not load one of the histograms. Leaving..." << endl;
   		return;
   	}

   	cout << "Fitting K0s with 1 TOF track with RooFit...";
   	//fit(histK0s1TOF, "RooFitK0s1", inputTag, false);
    crystalBall(histK0s1TOF, "RooFitCBK0s1", inputTag, false);

    cout << "Fitting K0s with 2 TOF tracks with RooFit...";
    //fit(histK0s2TOF, "RooFitK0s2", inputTag, true);
    crystalBall(histK0s2TOF, "RooFitCBK0s2", inputTag, true);

   	cout << "Fitting Lambda with 1 TOF track with RooFit...";
   	//fitLambda(histLambda1TOF, "RooFitLambda1", inputTag, false);
    crystalBallLambda(histLambda1TOF, "RooFitCBLambda1", inputTag, false);

    cout << "Fitting Lambda with 2 TOF tracks with RooFit...";
    //fitLambda(histLambda2TOF, "RooFitLambda2", inputTag, true);
    crystalBallLambda(histLambda2TOF, "RooFitCBLambda2", inputTag, true);


   	cout << "All plots created, Goodbye."<< endl;
    file->Close(); // Close the ROOT files

}

void fit(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF) {

    TString xAxisDescription;
    double textSize = 0.04;
    double a, b, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background, smallpeak;
    double meanAlternate, sigmaAlternate, ampAlternate;
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
        a = -61111;
        b = 55500;
        mean = 0.50;
        sigma = 0.01;
        amp = 4000;
        meanAlternate = 0.49;
        sigmaAlternate = 0.1;
        ampAlternate = 500;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 9000;
        yRangeBottom = 0;
        signal = 100000;
        background = 3000000;
        smallpeak = 20000;
        
    }else{

        a = -30111;
        b = 55500;
        mean = 0.50;
        sigma = 0.005;
        amp = 11000;
        meanAlternate = 0.49;
        sigmaAlternate = 0.1;
        ampAlternate = 2000;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 40000;
        yRangeBottom = 0;
        signal = 150000;
        background = 20000000;
        smallpeak = 200000;
    }

    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0, 1);
    RooGaussian gauss1("gauss1", "Gaussian 1", x, mean1, sigma1);

    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -100000, 1000000);
    RooRealVar coef1("coef1", "linear term", a, -1000000, 1000000);
    //RooRealVar coef2("coef2", "quadratic term", 250000, 0, 1000000);
    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2


    // Coefficients for PDFs
    RooRealVar nsig("nsig", "signal ", signal, 1000, 1000000);
    RooRealVar nbkg("nbkg", "Background",background, 1000, 600000000);


    RooAddPdf model("model", "Gaussians plus polynomial", RooArgList(gauss1, poly), RooArgList(nsig, nbkg)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, Save());

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    frame->GetXaxis()->SetTitle(xAxisDescription);
    frame->GetYaxis()->SetTitleOffset(1.5);
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(poly), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gauss1), LineStyle(kDashed), LineColor(kBlue));
    //frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);


    frame->Draw();

    // create a legend with results of fit
    TPaveText* text = new TPaveText(0.15, 0.85, 0.4, 0.65,"brNDC");
    text -> SetTextSize(textSize - 0.01);
    text -> SetTextAlign(11);
    text -> SetFillColor(0);
    text -> SetTextFont(62);
    int yield = ceil(nsig.getVal());
    double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield = " + to_string(yield);
    TString polCoef0 = "Pol0 = " + to_string(coef0.getVal());
    TString polCoef1 = "Pol1 = " + to_string(coef1.getVal());
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(polCoef0);
    text -> AddText(polCoef1);
    text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> Draw("same");

    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(textSize);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");

    TLegend *leg = new TLegend(0.7, 0.8 ,0.9, 0.55 );
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(textSize-0.01);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(frame->findObject("data"), "Data", "lep");
    leg -> AddEntry(frame->findObject("model"), "Total fit", "l");
    leg -> AddEntry(frame->findObject("gauss1"), "Gauss 1", "l");
    leg -> AddEntry(frame->findObject("gauss2"), "Gauss 2", "l");
    leg -> Draw("same");



    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file

}

void crystalBall(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF){

    TString xAxisDescription;
    double textSize = 0.04;
    double a, b,d, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background, smallpeak, n_guess, alpha_guess;
    //double meanAlternate, sigmaAlternate, ampAlternate;
    double legendPos[4];
    double resultPos[4];
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
        a = -61111;
        b = 55500;
        d = 5000;
        mean = 0.497;
        sigma = 0.004;
        amp = 4000;
        //meanAlternate = 0.49;
        //sigmaAlternate = 0.1;
        //ampAlternate = 500;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 9000;
        yRangeBottom = 0;
        signal = 140000;
        background = 6000000;
        n_guess = 4;
        alpha_guess = 1.5;

        resultPos[0] = 0.15;
        resultPos[1] = 0.85;
        resultPos[2] = 0.4;
        resultPos[3] = 0.65;
        legendPos[0] = 0.65;
        legendPos[1] = 0.72;
        legendPos[2] = 0.85;
        legendPos[3] = 0.55;

    }else{

        a = -30111;
        b = 55500;
        d = 5000;
        mean = 0.497;
        sigma = 0.0045;
        //amp = 11000;
        amp = 5000;

        //meanAlternate = 0.49;
        //sigmaAlternate = 0.1;
        //ampAlternate = 2000;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 35000;
        yRangeBottom = 0;
        signal = 200000;
        background = 40000000;
        n_guess = 4;
        alpha_guess = 1.5;
        resultPos[0] = 0.15;
        resultPos[1] = 0.6;
        resultPos[2] = 0.4;
        resultPos[3] = 0.3;
        legendPos[0] = 0.65;
        legendPos[1] = 0.6;
        legendPos[2] = 0.85;
        legendPos[3] = 0.4;
    }

    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , 0.495,0.51); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0.003, 0.005);
    RooRealVar alpha("alpha", "Tail parameter", alpha_guess, 1, 10);
    RooRealVar n("n", "Power parameter", n_guess, 3, 15);

    RooCBShape CBfunc("cb", "Crystal Ball PDF", x, mean1, sigma1, alpha, n);
    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -1000000, 10000000);
    RooRealVar coef1("coef1", "linear term", a, -10000000, 10000000);
    //RooRealVar coef2("coef2", "quadratic term", d, -1000000, 1000000);

    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2


    // Coefficients for PDFs
    RooRealVar nsig("nsig", "signal ", signal, 1000, 500000);
    RooRealVar nbkg("nbkg", "Background",background, 1000, 6e+10);


    RooAddPdf model("model", "Gaussians plus polynomial", RooArgList(CBfunc, poly), RooArgList(nsig, nbkg)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, PrintLevel(-1));

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    frame->GetXaxis()->SetTitle(xAxisDescription);
    frame->GetYaxis()->SetTitleOffset(1.5);
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(poly), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(CBfunc), LineStyle(kDashed), LineColor(kBlue));
    //frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);
    frame->Draw();


    // create a legend with results of fit
    TPaveText* text = new TPaveText(resultPos[0],resultPos[1], resultPos[2], resultPos[3],"brNDC");
    text -> SetTextSize(textSize - 0.01);
    text -> SetTextAlign(11);
    text -> SetFillColor(0);
    text -> SetTextFont(62);
    int yield = ceil(nsig.getVal());
    int yieldErr = ceil(nsig.getError());
    double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString nFinal = "n = " + to_string(n.getVal()) + " +/- " + to_string(n.getError());
    TString alphaFinal = "#alpha = " + to_string(alpha.getVal()) + " +/- " + to_string(alpha.getError());
    TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield = " + to_string(yield) + " +/- " + to_string(yieldErr);
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(nFinal);
    text -> AddText(alphaFinal);
    text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> Draw("same");

    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(textSize);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");

    int j = frame->numItems();
    cout << "numItems: " << j << endl;
    TString names[j];
    TObject* obj[j];
    cout << "Frame objects:\n";
    j = 0;
    for (int i=0; i<frame->numItems(); i++) {
        TString obj_name=frame->nameOf(i); 
        if (obj_name=="") 
            continue;
        cout << Form("%d. '%s'\n",i,obj_name.Data());
        names[j] = Form("%d. '%s'\n",i,obj_name.Data());
        obj[j] = frame->findObject(obj_name.Data());
        ++j;
    }


    TLegend *leg = new TLegend(legendPos[0], legendPos[1], legendPos[2], legendPos[3]);
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(textSize-0.01);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);
    leg ->AddEntry(obj[0], "Data", "lep");
    leg ->AddEntry(obj[1], "Fit function", "l");
    leg ->AddEntry(obj[2], "Poly1", "l");
    leg ->AddEntry(obj[3], "Crystal ball", "l");


    leg -> Draw("same");
    

    TPaveText* text3 = new TPaveText( 0.6, 0.8, 0.85, 0.8,"brNDC");
    text3 -> SetTextSize(textSize-0.01);    
    text3 -> SetFillColor(0);
    text3 -> SetTextFont(62);
    TString tracksText;
    if(is2TOF){
        tracksText = "2 ToF tracks";
    }else{
        tracksText = "at least 1 ToF track";
    }
    text3 -> AddText(tracksText);
    text3 -> Draw("same");


    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file

}


void fitLambda(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF) {

    TString xAxisDescription;
    double textSize = 0.04;
    double a, b, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background;
    double meanAlternate, sigmaAlternate, ampAlternate;
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
    a = 50000;
    b = -54500;
    mean = 1.115;
    sigma = 0.001;
    amp = 700;
    xAxisDescription = "m_{p #pi} [GeV/c^{2}]";
    yRangeTop = 1200;
    yRangeBottom = 0;
    signal = 5000;
    background = 5000;
    } else {
    a = 40000;
    b = -40500;
    mean = 1.115;
    sigma = 0.001;
    amp = 2000;
    xAxisDescription = "m_{p #pi} [GeV/c^{2}]";
    yRangeTop = 1200;
    yRangeBottom = 0;
    signal = 10000;
    background = 500000;
    }
    

    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0, 0.1);
    RooRealVar amplitude1("amplitude1", "Amplitude1 of Gaussian", amp, 0, 10000);
    RooGaussian gauss1("gauss1", "Gaussian 1", x, mean1, sigma1);

    RooGenericPdf gaussWithAmp1("gaussWithAmp1", "Amplitude * Gauss", "amplitude1 * gauss1", RooArgList(amplitude1, gauss1));

    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -100000 , 100000);
    RooRealVar coef1("coef1", "linear term", a, -100000, 100000);
    //RooRealVar coef2("coef2", "quadratic term", 250000, 0, 1000000);
    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2

    RooRealVar nsig("nsig", "signal ", signal, 1000, 50000);
    RooRealVar nbkg("nbkg", "Background",background, 2000, 100000);
    // Coefficients for PDFs
    RooRealVar gaussFrac1("gaussFrac1", "fraction of Gauss 1", 0.54, 0, 1);
    RooAddPdf model("model", "Gaussians plus polynomial", RooArgList(gauss1, poly), RooArgList(nsig, nbkg)); //,poly   gaussFrac2



    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, Save());

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    frame->GetXaxis()->SetTitle(xAxisDescription);
    frame->GetYaxis()->SetTitleOffset(1.5);
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(poly), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gauss1), LineStyle(kDashed), LineColor(kBlue));
    //model.plotOn(frame, Components(gaussWithAmp2), LineStyle(kDashed), LineColor(kGreen));
    //frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);


    frame->Draw();

    // create a legend with results of fit
    TPaveText* text = new TPaveText(0.15, 0.85, 0.4, 0.65,"brNDC");
    text -> SetTextSize(textSize - 0.01);
    text -> SetTextAlign(11);
    text -> SetFillColor(0);
    text -> SetTextFont(62);
    int yield = ceil(nsig.getVal());
    double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield = " + to_string(yield);
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> Draw("same");

    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(textSize);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");

    TPaveText* text3 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text3 -> SetTextSize(textSize);
    text3 -> SetFillColor(0);
    text3 -> SetTextFont(62);
    TString tracksText;
    if(is2TOF){
        tracksText = "2 ToF tracks";
    }else{
        tracksText = "at least 1 ToF track";
    }
    text3 -> AddText(tracksText);
    text3 -> Draw("same");


    TLegend *leg = new TLegend(0.7, 0.8 ,0.9, 0.55 );
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(textSize-0.01);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(frame->findObject("data"), "Data", "lep");
    leg -> AddEntry(frame->findObject("model"), "Total fit", "l");
    leg -> AddEntry(frame->findObject("gauss1"), "Gauss 1", "l");
    //leg -> AddEntry(frame->findObject("gauss2"), "Gauss 2", "l");
    leg -> Draw("same");



    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file

}



void crystalBallLambda(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF){
    TString xAxisDescription;
    double textSize = 0.04;
    double a, b, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background, n_guess, alpha_guess;
    double meanAlternate, sigmaAlternate, ampAlternate;
    double legendPos[4];
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
        a = 50000;
        b = -54500;
        mean = 1.115;
        sigma = 0.001;
        amp = 700;
        xAxisDescription = "m_{p #pi} [GeV/c^{2}]";
        yRangeTop = 1200;
        yRangeBottom = 0;
        signal = 5000;
        background = 5000;
        n_guess = 1;
        alpha_guess = 1;
        legendPos[0] = 0.8;
        legendPos[1] = 0.75;
        legendPos[2] = 0.6;
        legendPos[3] = 0.6;

    } else {
        a = 40000;
        b = -40500;
        mean = 1.115;
        sigma = 0.001;
        amp = 2000;
        xAxisDescription = "m_{p #pi} [GeV/c^{2}]";
        yRangeTop = 1200;
        yRangeBottom = 0;
        signal = 10000;
        background = 500000;
        n_guess = 1;
        alpha_guess = 1;
        legendPos[0] = 0.6;
        legendPos[1] = 0.75;
        legendPos[2] = 0.8;
        legendPos[3] = 0.6;
    }
    
    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0, 1);
    RooRealVar alpha("alpha", "Tail parameter", alpha_guess, 0, 10);
    RooRealVar n("n", "Power parameter", n_guess, 0, 1e+6);

    RooCBShape CBfunc("cb", "Crystal Ball PDF", x, mean1, sigma1, alpha, n);
    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -100000, 1000000);
    RooRealVar coef1("coef1", "linear term", a, -1000000, 1000000);
    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2


    // Coefficients for PDFs
    RooRealVar nsig("nsig", "signal ", signal, 1000, 1000000);
    RooRealVar nbkg("nbkg", "Background",background, 1000, 600000000);


    RooAddPdf model("model", "crystal ball plus polynomial", RooArgList(CBfunc, poly), RooArgList(nsig, nbkg)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, Save());

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    frame->GetXaxis()->SetTitle(xAxisDescription);
    frame->GetYaxis()->SetTitleOffset(1.5);
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(poly), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(CBfunc), LineStyle(kDashed), LineColor(kBlue));
    //frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);
    frame->Draw();

    // create a legend with results of fit
    TPaveText* text = new TPaveText(0.15, 0.85, 0.4, 0.65,"brNDC");
    text -> SetTextSize(textSize - 0.01);
    text -> SetTextAlign(11);
    text -> SetFillColor(0);
    text -> SetTextFont(62);
    int yield = ceil(nsig.getVal());
    int yieldErr = ceil(nsig.getError());
    double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString nFinal = "n = " + to_string(n.getVal()) + " +/- " + to_string(n.getError());
    TString alphaFinal = "#alpha = " + to_string(alpha.getVal()) + " +/- " + to_string(alpha.getError());
    TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield = " + to_string(yield) + " +/- " + to_string(yieldErr);
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(nFinal);
    text -> AddText(alphaFinal);
    text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> Draw("same");

    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(textSize);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");


    int j = frame->numItems();
    cout << "numItems: " << j << endl;
    TString names[j];
    TObject* obj[j];
    cout << "Frame objects:\n";
    j = 0;
    for (int i=0; i<frame->numItems(); i++) {
        TString obj_name=frame->nameOf(i); 
        if (obj_name=="") 
            continue;
        cout << Form("%d. '%s'\n",i,obj_name.Data());
        names[j] = Form("%d. '%s'\n",i,obj_name.Data());
        obj[j] = frame->findObject(obj_name.Data());
        ++j;
    }

    TLegend *leg = new TLegend(legendPos[0], legendPos[1], legendPos[2], legendPos[3] );
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(textSize-0.01);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg ->AddEntry(obj[0], "Data", "lep");
    leg ->AddEntry(obj[1], "Fit function", "l");
    leg ->AddEntry(obj[2], "Poly1", "l");
    leg ->AddEntry(obj[3], "Crystal ball", "l");
    leg -> Draw("same");

    TPaveText* text3 = new TPaveText( 0.6, 0.79, 0.85, 0.79,"brNDC");
    text3 -> SetTextSize(textSize-0.01);    
    text3 -> SetFillColor(0);
    text3 -> SetTextFont(62);
    TString tracksText;
    if(is2TOF){
        tracksText = "2 ToF tracks";
    }else{
        tracksText = "at least 1 ToF track";
    }
    text3 -> AddText(tracksText);
    text3 -> Draw("same");



    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file

}



double GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data){
    double chiSquare = frame->chiSquare();  // This calculates chi-square per degree of freedom

    int nBins = frame->GetNbinsX();  // Number of bins used in the histogram
    int nParams = model.getParameters(data)->getSize();  // Number of floating parameters in the model
    int ndf = nBins - nParams;

    return chiSquare/ndf;
}


double integrateGaussian(RooRealVar& x, RooGenericPdf& gauss, RooRealVar& mean, RooRealVar& sigma, double binWidth) {
    // Define the range for integration as mean ± 3sigma
    double meanVal = mean.getVal();
    double sigmaVal = sigma.getVal();
    double lowRange = meanVal - 3 * sigmaVal;
    double highRange = meanVal + 3 * sigmaVal;

    // Set the range on the variable
    x.setRange("peakRange", lowRange, highRange);


    // Integrate the Gaussian over the specified range
    RooAbsReal* integral = gauss.createIntegral(x, RooFit::NormSet(x), RooFit::Range("peakRange"));
    double result = integral->getVal(); // This gives the fraction of the total Gaussian within the range

    // To get the total number of events under the peak, multiply by the number of expected events
    // Assuming 'yield' is the RooRealVar representing the number of events (amplitude of Gaussian)
    // This should be defined in your model if you are using a normalized PDF for fitting
    double totalYield = result / binWidth;

    cout << "--------------------------------------------------" << endl;
    cout << "meanVal: " << meanVal << endl;
    cout << "sigmaVal: " << sigmaVal << endl;
    cout << "lowRange: " << lowRange << endl;
    cout << "highRange: " << highRange << endl;
    cout << "result: " << result << endl;
    cout << "binWidth: " << binWidth << endl;
    cout << "totalYield: " << totalYield << endl;
    cout << "--------------------------------------------------" << endl;

    delete integral; // Clean up to avoid memory leaks
    return totalYield;
}

void fitWith2Gauss(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF) {

    TString xAxisDescription;
    double textSize = 0.04;
    double a, b, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background, smallpeak;
    double meanAlternate, sigmaAlternate, ampAlternate;
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
        a = -61111;
        b = 55500;
        mean = 0.50;
        sigma = 0.01;
        amp = 4000;
        meanAlternate = 0.49;
        sigmaAlternate = 0.1;
        ampAlternate = 500;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 9000;
        yRangeBottom = 0;
        signal = 100000;
        background = 3000000;
        smallpeak = 20000;
        
    }else{

        a = -30111;
        b = 55500;
        mean = 0.50;
        sigma = 0.005;
        amp = 11000;
        meanAlternate = 0.49;
        sigmaAlternate = 0.1;
        ampAlternate = 2000;
        xAxisDescription = "m_{#pi^{+} #pi^{-}} [GeV/c^{2}]";
        yRangeTop = 40000;
        yRangeBottom = 0;
        signal = 150000;
        background = 20000000;
        smallpeak = 200000;
    }

    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0, 1);
    RooRealVar amplitude1("amplitude1", "Amplitude1 of Gaussian", amp, 0, 20000);
    RooGaussian gauss1("gauss1", "Gaussian 1", x, mean1, sigma1);
    
    RooGenericPdf gaussWithAmp1("gaussWithAmp1", "Amplitude * Gauss", "amplitude1 * gauss1", RooArgList(amplitude1, gauss1));


    RooRealVar mean2("mean2", "mean of Gaussian 2", meanAlternate,hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // adjust these ranges and initial values as needed
    RooRealVar sigma2("sigma2", "width of Gaussian 2", sigmaAlternate,0,1);
    RooRealVar amplitude2("amplitude2", "amplitude of Gaussian 2", ampAlternate, 0, 5000);
    RooGaussian gauss2("gauss2", "Gaussian 2", x, mean2, sigma2);

    RooGenericPdf gaussWithAmp2("gaussWithAmp2", "Amplitude * Gauss", "amplitude2 * gauss2", RooArgList(amplitude2, gauss2));
    
    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -100000, 1000000);
    RooRealVar coef1("coef1", "linear term", a, -1000000, 1000000);
    //RooRealVar coef2("coef2", "quadratic term", 250000, 0, 1000000);
    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2


    // Coefficients for PDFs
    RooRealVar gaussFrac1("gaussFrac1", "fraction of Gauss 1", 0.3, 0, 1);
    RooRealVar gaussFrac2("gaussFrac2", "fraction of Gauss 2", 0.2, 0, 1);
    RooAddPdf model("model", "Gaussians plus polynomial", RooArgList(gaussWithAmp1, gaussWithAmp2, poly), RooArgList(gaussFrac1,gaussFrac2)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, Save());

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    frame->GetXaxis()->SetTitle(xAxisDescription);
    frame->GetYaxis()->SetTitleOffset(1.5);
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(poly), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(gauss1), LineStyle(kDashed), LineColor(kBlue));
    //model.plotOn(frame, Components(gaussWithAmp2), LineStyle(kDashed), LineColor(kGreen));
    frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);


    frame->Draw();

    // create a legend with results of fit
    TPaveText* text = new TPaveText(0.15, 0.85, 0.4, 0.65,"brNDC");
    text -> SetTextSize(textSize - 0.01);
    text -> SetTextAlign(11);
    text -> SetFillColor(0);
    text -> SetTextFont(62);
    int yield = integrateGaussian(x, gaussWithAmp1, mean1, sigma1, hist->GetBinWidth(3));
    double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield = " + to_string(yield);
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> Draw("same");

    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(textSize);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");

    TLegend *leg = new TLegend(0.7, 0.8 ,0.9, 0.55 );
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(textSize-0.01);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(frame->findObject("data"), "Data", "lep");
    leg -> AddEntry(frame->findObject("model"), "Total fit", "l");
    leg -> AddEntry(frame->findObject("gauss1"), "Gauss 1", "l");
    leg -> AddEntry(frame->findObject("gauss2"), "Gauss 2", "l");
    leg -> Draw("same");



    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file

}



// staci 1 state per event-> hodit si do recTree bool ToF match, mat dcery rozdelene na tag/probe
// 
// pseudorap, phi, pT distribucie do plotov
// 