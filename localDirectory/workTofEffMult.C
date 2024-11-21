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

TFile *file;
TTree *tree;


Double_t outerBorder[2] = {13, 37};
Double_t innerBorder[2] = {19, 31}; 

TH1D* getHistK0s(bool is2TOF, bool onlyFirst );
void fitK0s();
void crystalBall(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF);

void workTofEffMult(){

	file = new TFile("TofEffMult/StRP_production_0000.root");

	if(!file){
		cout << "Did not get file." << endl;
		return;
	}

	tree = (TTree*)file->Get("recTreeTofEffMult");

	if(!tree){
		cout << "Did not get tree" << endl;
		return;
	}


	cout << "Obtaining data from tree into a histogram..." << endl;

	TH1D* hist1 = getHistK0s(false, false);
	TH1D* hist2 = getHistK0s(true, false);

	if(hist1 == nullptr || hist2 == nullptr)
		return;

	cout << "Now fitting histograms..." << endl;

	crystalBall(hist1, "invMass1All", "TofEffMult", false);
	crystalBall(hist2, "invMass2All", "TofEffMult", true);



	cout << "Obtaining data from tree into a histogram for only the first state..." << endl;

	hist1 = getHistK0s(false, true);
	hist2 = getHistK0s(true, true);

	if(hist1 == nullptr || hist2 == nullptr)
		return;

	cout << "Now fitting histograms..." << endl;

	crystalBall(hist1, "invMass1justFirst", "TofEffMult", false);
	crystalBall(hist2, "invMass2justFirst", "TofEffMult", true);


}



TH1D* getHistK0s(bool is2TOF, bool onlyFirst ) {

	int numBins = 45;
	Double_t minRange = 0.45;
	Double_t maxRange = 0.54;

    TH1D *signalFinal = new TH1D("fitPeakK0sSignalFit", "fitPeakK0sSignalFit", numBins, minRange, maxRange);

    signalFinal->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    
    TString cmd, condition;

    //draw the pairs invMasses which satisfy condition
    for (int i = 0; i < 5; ++i){


	   	tree->Draw(TString::Format("invMass%d>>hist(%d, %f, %f)", i,numBins, minRange,maxRange), TString::Format("pairID%d == 0 && tofHit%d > 0", i, 2*i +1) );
		signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
        cout << "Number of entries: " << signalFinal->GetEntries() << endl;
	   	tree->Draw(TString::Format("invMass%d>>hist(%d, %f, %f)", i, numBins, minRange,maxRange), TString::Format("pairID%d == 0 && tofHit%d > 0", i, 2*i) );
	   	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
		

		if(!is2TOF){
	   		tree->Draw(TString::Format("invMass%d>>hist(%d, %f, %f)",i, numBins, minRange,maxRange),TString::Format("pairID%d == 0 && tofHit%d < 0", i, 2*i +1) );
	   		signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
		}

		if(onlyFirst && i == 0){
			return signalFinal;
		}

    }

    if(signalFinal->GetEntries() == 0){
    	cout << "Empty signalFinal. " << endl;
    	return nullptr;
    }

    return signalFinal;

}



void crystalBall(TH1D*& hist, const TString& outName, const TString& inputTag, bool is2TOF){

    TString xAxisDescription;
    double textSize = 0.04;
    double a, b,d, mean, sigma, amp, yRangeTop, yRangeBottom, signal, background, smallpeak, n_guess, alpha_guess;
    //double meanAlternate, sigmaAlternate, ampAlternate;
    double legendPos[4];
    double resultPos[4], yAxisText[2];
    //define all variables which are different for K0s and Lambda
    if(is2TOF){
        a = -61111;
        b = 5500;
        d = -200000;
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

        resultPos[0] = 0.12;
        resultPos[1] = 0.85;
        resultPos[2] = 0.4;
        resultPos[3] = 0.65;
        legendPos[0] = 0.6;
        legendPos[1] = 0.72;
        legendPos[2] = 0.8;
        legendPos[3] = 0.6;
        yAxisText[0] = 0.74;
        yAxisText[1] = 0.8;

    }else{

        a = -30111;
        b = 55500;
        d = -10000;
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
        resultPos[0] = 0.12;
        resultPos[1] = 0.6;
        resultPos[2] = 0.4;
        resultPos[3] = 0.25;
        legendPos[0] = 0.6;
        legendPos[1] = 0.6;
        legendPos[2] = 0.8;
        legendPos[3] = 0.45;
        yAxisText[0] = 0.74;
        yAxisText[1] = 0.8;
    }

    // Define the observable
    RooRealVar x("x", "Observable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Convert TH1D to RooDataHist
    RooDataHist data("data", "dataset with x", x, Import(*hist));

    // Define the Gaussian PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", mean , 0.495,0.51); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", sigma , 0.003, 0.005);
    RooRealVar alpha("alpha", "Tail parameter", alpha_guess, 1, 10);
    RooRealVar n("n", "Power parameter", n_guess, 3, 11);

    RooCBShape CBfunc("cb", "Crystal Ball PDF", x, mean1, sigma1, alpha, n);
    // Define the polynomial PDF for the background
    RooRealVar coef0("coef0", "constant term", b, -1000000, 10000000);
    RooRealVar coef1("coef1", "linear term", a, -10000000, 10000000);
    //RooRealVar coef2("coef2", "quadratic term", d, -1000000, 1000000);

    RooPolynomial poly("poly", "Background Polynomial", x, RooArgList(coef0, coef1)); // , coef2


    // Coefficients for PDFs
    RooRealVar nsig("nsig", "signal ", signal, 1000, 2000000);
    RooRealVar nbkg("nbkg", "Background",background, 1000, 6e+10);


    RooAddPdf model("model", "Gaussians plus polynomial", RooArgList(CBfunc, poly), RooArgList(nsig, nbkg)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data, PrintLevel(-1));


    double sumSignal = 0;
    double sumSignalError = 0;
    double sumBcg = 0;
    double sumBcgError = 0;
    // integrate signal when approximating linear rise of background

    for (int iBin = outerBorder[0]; iBin < outerBorder[1] ; ++iBin){ // integrating background: 0.46-0.48 and 0.52-0.54, signal 0.48-0.52
        
        if(iBin >= innerBorder[0]  && iBin < innerBorder[1] ){
            sumSignal += hist->GetBinContent(iBin);
            sumSignalError += pow(hist->GetBinError(iBin),2);
            cout << "just added bin number " << iBin << " to signal" << endl;
        }else if(iBin < innerBorder[0] || iBin >= innerBorder[1] ){
            sumBcg += hist->GetBinContent(iBin);
            sumBcgError += pow(hist->GetBinError(iBin),2);
            cout << "just added bin number " << iBin << " to bcg" << endl;
        }
    }


    int sum = sumSignal - sumBcg;
    int sumError = sqrt(sumSignalError) + sqrt(sumBcgError);


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
    //double goodyFit = GoodnessOfFit(frame, model, data);
    TString meanFinal = "#mu_{1} = " + to_string(mean1.getVal()) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString nFinal = "n = " + to_string(n.getVal()) + " +/- " + to_string(n.getError());
    TString alphaFinal = "#alpha = " + to_string(alpha.getVal()) + " +/- " + to_string(alpha.getError());
    //TString fitFinal = "#chi^{2}/NDF = " + to_string(goodyFit);
    TString yieldText = "Yield (fitting) = " + to_string(yield) + " +/- " + to_string(yieldErr);
    TString integrationText = "Yield (integration) = " + to_string(sum) + " +/- " + to_string(sumError);
    text -> AddText(meanFinal);
    text -> AddText(sigmaFinal);
    text -> AddText(nFinal);
    text -> AddText(alphaFinal);
    //text -> AddText(fitFinal);
    text -> AddText(yieldText);
    text -> AddText(integrationText);
    //text -> AddText("a_0 = " + TString(to_string(coef0.getVal() )) + " +/- " + TString(to_string(coef0.getError())) );
    //text -> AddText("a_1 = " + TString(to_string(coef1.getVal() )) + " +/- " + TString(to_string(coef1.getError())) );
    //text -> AddText( "a_2 = " + TString(to_string(coef2.getVal() )) + " +/- " + TString(to_string(coef2.getError())) );
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
    

    TPaveText* text3 = new TPaveText( legendPos[0], yAxisText[0], legendPos[2], yAxisText[1],"brNDC");
    text3 -> SetTextSize(textSize-0.01);    
    text3 -> SetFillColor(0);
    text3 -> SetTextAlign(11);  
    text3 -> SetTextFont(62);
    TString tracksText;
    //if(globalTracks){
    //    tracksText = "using global tracks";
    //}else{
    //    tracksText = "using primary tracks only";
    //}
    text3 -> AddText(tracksText);
    if(is2TOF){
        tracksText = "2 ToF tracks";
    }else{
        tracksText = "at least 1 ToF track";
    }
    text3 -> AddText(tracksText);
    text3 -> Draw("same");


    c->SaveAs(inputTag + "/" + outName + ".pdf");  // Save the canvas to a file
    c->Close();
}