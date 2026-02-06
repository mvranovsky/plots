#include "PlotZeroBias.h"


PlotZeroBias::PlotZeroBias(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotZeroBias::PlotZeroBias(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotZeroBias::Make() {


    if(DEBUGMODE) cout << "Starting PlotZeroBias::Make()..." << endl;

    //cout << "Now calculating trigger efficiency for varying pT..." << endl;
    calculateTriggerEfficiency("pT");

    if(DEBUGMODE) cout << "Now calculating trigger efficiency for varying invmass..." << endl;
    calculateTriggerEfficiency("mee");

    if(DEBUGMODE) cout << "Now creating invMass dependence plot..." << endl;
    invMassDependence2();

    if(DEBUGMODE) cout << "Finished PlotZeroBias::Make()..." << endl;
    handleHistograms(nameOfAnaZeroBiasDir, "ZBData");

}


void PlotZeroBias::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();

    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE) cout << "The output file is saved: " << outputPosition << endl;
}

void PlotZeroBias::Init(){
	//define the output file which will store all the canvases
    if(DEBUGMODE) cout << "Starting PlotZeroBias::Init()..." << endl;
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	if(!outFile || outFile->IsZombie() ){
		cerr << "ERROR in PlotZeroBias::Init(): Couldn't open output file. Leaving..." << endl;
        return;
	}

    outFile->cd();
    outFile->mkdir(nameOfAnaZeroBiasDir);
    outFile->cd();
    
    if(DEBUGMODE) cout << "Opening histogram file in PlotZeroBias::Init()..." << endl;
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );
    
    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "ERROR in PlotZeroBias::Init(): Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
    if(DEBUGMODE) cout << "Connecting input tree for zerobias..." << endl;
	ConnectInputTree(inputPosition, nameOfAnaZeroBiasTree, tree, bcgTree);
    if(!tree){
    	cerr << "ERROR in PlotZeroBias::Init(): Couldn't open tree with data. Returning." << endl;
    	return;
    }

    mTriggerEffError = 0.0;

    if(DEBUGMODE) cout << "Finished PlotZeroBias::Init()..." << endl;
}

void PlotZeroBias::calculateTriggerEfficiency(TString var){

    if(DEBUGMODE) cout << "Calculating trigger efficiency from ZeroBias data for variable " << var << endl;

    TH1 *hist = new TH1D("h_" + var, "trigger efficiency; #epsilon_{topology} [%]; counts", 40,0,0.8);

    
    for(int i = 0; i < 20; i++){
        double low = getLow(var, i);
        for(int j = 0; j < 20; j++){

            double top = getTop(var, j);
            if(top <= low){
                cerr << "ERROR in PlotZeroBias::calculateTriggerEfficiency(): This should not be happening. Leaving." << endl;
                return;
            }
            TString cond = "";
            if(var.Contains("pT")){
                cond = Form("pTInGev0 > %f && pTInGev0 < %f && pTInGev1 > %f && pTInGev1 < %f", low, top, low, top); // limit to pT < 1.8 GeV/c
            }
            else if(var.Contains("mee")){
                cond = Form("invMass > %f && invMass < %f", low, top); // limit to mee < 3.3 GeV/c^2
            }
            double eff = triggerEfficiency(cond);
            hist->Fill(eff);

        }
    }


    CreateCanvas(&canvas, TString("TriggerEfficiencySysStudy_") + var, 1200, 800);
    SetHistStyle(hist, kBlue, markerStyleTypical);
    hist->Draw("hist");
    hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);

    // fit the distribution with a Gaussian
    TF1 *gauss = new TF1("gauss", "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    gauss->SetLineColor(kRed);
    gauss->SetLineWidth(2);
    hist->Fit(gauss, "R", "", 0, 1.0);
    gauss->Draw("same");

    //print the fit results onto the canvas
    CreateText(Form("#mu = %.4f +/- %.4f", gauss->GetParameter(1), gauss->GetParError(1)),textFont, 0.72, 0.75, 0.9, 0.72);
    CreateText(Form("#sigma = %.4f +/- %.4f", gauss->GetParameter(2), gauss->GetParError(2)), textFont, 0.72, 0.72, 0.9, 0.69);

    DrawSTARInternal(0.77, 0.89, 0.9, 0.93);

    
    CreateLegend(&legend, 0.2, 0.8, 0.45, 0.73);
    legend->AddEntry(hist, "Trigger topology efficiency", "lep");
    legend->AddEntry(gauss, "Gaussian Fit", "l");
    legend->Draw("same");

    mTriggerEfficiency = gauss->GetParameter(1);
    mTriggerEffError = gauss->GetParError(1);
    mTriggerEffSysError = gauss->GetParameter(2);
    
    if(DEBUGMODE) cout << "Trigger topology efficiency from zerobias data: " << mTriggerEfficiency << " +/- " << mTriggerEffError << endl;

    canvas->Update();
    outFile->cd();
    outFile->cd(nameOfAnaZeroBiasDir);
    canvas->SetName(TString("triggerEfficiencyGaussFit_") + var);
    canvas->Write();
    canvas->Clear();

    return;
}

double PlotZeroBias::getLow(TString var, int i){

    if(var.Contains("pT")){
        return 0.8 + i*0.02;
    }
    else if(var.Contains("mee")){
        return 2.0 + i*0.04;
    }
    else{
        cerr << "ERROR in PlotZeroBias::getLow(): Invalid variable name. Returning -1." << endl;
        return -1;
    }
}

double PlotZeroBias::getTop(TString var, int i){

    if(var == "pT"){
        return 1.8 + i*0.02;
    }
    else if(var == "mee"){
        return 3.3 + i*0.04;
    }else{
        cerr << "ERROR in PlotZeroBias::getTop(): Invalid variable name. Returning -1." << endl;
        return -1;
    }
}
double PlotZeroBias::triggerEfficiency(TString cond){


    int nBins = 3;
    double lowerLim = -1.5;
    double upperLim = 1.5;

    //cout << "Plotting trigger efficiency with condition: " << cond << endl;
    tree->Draw(Form("isTopology>>hist(%d, %f, %f)", nBins, lowerLim, upperLim), cond);
    TH1D* hist = (TH1D*)gDirectory->Get("hist");


    if(!hist || hist->GetEntries() == 0){
        cerr << "ERROR in PlotZeroBias::triggerEfficiency(): Could not get histogram. Returning." << endl;
        return -1;
    }

    int minus = hist->GetBinContent(1);
    int plus = hist->GetBinContent(3);

    double minusErr = hist->GetBinError(1);
    double plusErr = hist->GetBinError(3);

    if(plus == 0 || minus == 0){
        cerr << "ERROR in PlotZeroBias::triggerEfficiency(): No entries in plus or minus bin. Returning." << endl;
        return -1;
    }

    // calculate error
    double errCombined = sqrt(pow(plusErr, 2) + pow(minusErr, 2));

    setNumerator(plus);
    setDenominator(plus + minus);
    setNumeratorError(plusErr);
    setDenominatorError(errCombined);


    return (1.0* plus/ (plus + minus));
}


void PlotZeroBias::invMassDependence(){

    CreateCanvas(&canvas, "invMassDependence");
    canvas->SetName("invMassDependence");
    SetGPad();

    TH1D *hNumerator = new TH1D("hNumerator", "Numerator", 40, 0, 4);
    TH1D *hDenominator = new TH1D("hDenominator", "Denominator", 40, 0, 4);

    Double_t invMass[40], eff[40], effErr[40], invMassErr[40];
    // Fill the graph with data points
    for (int j = 0; j < 40; j++) {
        double mass = j * 0.2;
        double efficiency = triggerEfficiency(Form("invMass > %f", mass));
        invMass[j] = mass;
        eff[j] = efficiency;
        effErr[j] = mTriggerEffError;
        invMassErr[j] = 0.0;  // Example error

        hNumerator->SetBinContent(j + 1, getNumerator());
        hNumerator->SetBinError(j + 1, getNumeratorError());
        hDenominator->SetBinContent(j + 1, getDenominator());
        hDenominator->SetBinError(j + 1, getDenominatorError());
    }

    TEfficiency *teff = new TEfficiency(*hNumerator, *hDenominator);
    TGraphAsymmErrors* g = teff->CreateGraph();
    g->SetName("invMassDependence");
    g->SetTitle("");
    SetTGraphStyle(g);
    g->GetXaxis()->SetTitle("m_{pairs} [GeV/c^{2}]");
    g->GetYaxis()->SetTitle("Trigger Efficiency");

    // set x limits to 0 and 4
    g->GetXaxis()->SetLimits(0, 4);
    g->SetMaximum(1.0);
    g->SetMinimum(0.0);

    g->Draw("ALP");
    DrawSTARInternal();


    outFile->cd();
    outFile->cd(nameOfAnaZeroBiasDir);
    canvas->SetName("TGraphAssymErrorEfficiency");
    canvas->Write();

    canvas->Clear();

    TGraphErrors *gr = new TGraphErrors(40, invMass, eff, invMassErr, effErr);
    gr->SetName("invMassDependence");
    gr->SetTitle("");
    SetTGraphStyle(gr);
    gr->GetXaxis()->SetTitle("m_{pairs} [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Trigger Efficiency");

    // set x limits to 0 and 4
    gr->GetXaxis()->SetLimits(0, 4);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.0);

    gr->Draw("ALP");
    DrawSTARInternal();


    canvas->SetName("TGraphErrorEfficiency");
    canvas->Write();

    return;
}



void PlotZeroBias::invMassDependence2(){

    CreateCanvas(&canvas, "invMassDependence");
    canvas->SetName("invMassDependence");
    SetGPad();

    Double_t invMass[5] = {2.85, 2.9, 3.0, 2.6 , 3.0};
    Double_t invMassErr[5] = {0.35, 0.4, 0.5, 0.6, 1.0};
    Double_t eff[5], effErr[5];
    for(int i = 0; i< 5; ++i){
        eff[i] = triggerEfficiency(Form("invMass > %f && invMass < %f", invMass[i]- invMassErr[i], invMass[i]+ invMassErr[i]));
        effErr[i] = mTriggerEffError;
    }
    

    TGraphErrors *gr = new TGraphErrors(5, invMass, eff, invMassErr, effErr);
    gr->SetName("invMassDependence");
    gr->SetTitle("");
    SetTGraphStyle(gr);
    gr->GetXaxis()->SetTitle("m_{pairs} [GeV/c^{2}]");
    gr->GetYaxis()->SetTitle("Trigger Efficiency");
    
    // set x limits to 0 and 4
    gr->GetXaxis()->SetLimits(2, 4);
    gr->SetMaximum(0.5);
    gr->SetMinimum(0.3);
    
    gr->Draw("AEP");
    DrawSTARInternal();
    
    
    outFile->cd();
    outFile->cd(nameOfAnaZeroBiasDir);
    canvas->SetName("invMassWindow");
    canvas->Write();

    return;
}