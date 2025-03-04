#include "../include/PlotAnaJPsi.h"


PlotAnaJPsi::PlotAnaJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotAnaJPsi::Make(){


	int nBins = 40;
    double lowerLim = 2;
    double upperLim = 4;
    roofit_invMassFit_Misa(nBins, lowerLim, upperLim);
    cout << "Finished Misa's plot" << endl;

    TString nameOfSignal = "hInvMassJPsi_beforeRP";
    TString nameOfBcg = "hInvMassJPsiBcg_beforeRP";
    invMassJPsi(nBins, lowerLim, upperLim,nameOfSignal, nameOfBcg );
    cout << "Finished JPsi plot" << endl;



    // save all the histograms to canvases into outfile
    handleHistograms();
    
    // get control plots and save them as histograms
    TH1D *hSignalControl = (TH1D*)histFile->Get("hInvMassJPsi");
    TH1D* hBcgControl = (TH1D*)histFile->Get("hInvMassJPsiBcg");

    if(hSignalControl && hBcgControl){
        outFile->cd();
        hSignalControl->Write("hInvMassJPsi_Control");
        hBcgControl->Write("hInvMassJPsiBcg_Control");
    }else{
        cout << "Could not get control invariant mass histograms." << endl;
    }


    // save pT plots
    TH1D *hPtMissing = (TH1D*)histFile->Get("hPtMissing");
    TH1D *hPtMissingBcg = (TH1D*)histFile->Get("hPtMissingBcg");

    TH1D *hPhotonAMomX = (TH1D*)histFile->Get("hPhotonAMomX");
    TH1D *hPhotonAMomY = (TH1D*)histFile->Get("hPhotonAMomY");
    TH1D *hPhotonAMomXBcg = (TH1D*)histFile->Get("hPhotonAMomXBcg");
    TH1D *hPhotonAMomYBcg = (TH1D*)histFile->Get("hPhotonAMomYBcg");

    TH1D *hPhotonBMomX = (TH1D*)histFile->Get("hPhotonBMomX");
    TH1D *hPhotonBMomY = (TH1D*)histFile->Get("hPhotonBMomY");
    TH1D *hPhotonBMomXBcg = (TH1D*)histFile->Get("hPhotonBMomXBcg");
    TH1D *hPhotonBMomYBcg = (TH1D*)histFile->Get("hPhotonBMomYBcg");

    outFile->cd();
    hPtMissing->Write(hPtMissing->GetName());
    hPtMissingBcg->Write(hPtMissingBcg->GetName());

    hPhotonAMomX->Write(hPhotonAMomX->GetName());
    hPhotonAMomY->Write(hPhotonAMomY->GetName());
    hPhotonAMomXBcg->Write(hPhotonAMomXBcg->GetName());
    hPhotonAMomYBcg->Write(hPhotonAMomYBcg->GetName());

    hPhotonBMomX->Write(hPhotonBMomX->GetName());
    hPhotonBMomY->Write(hPhotonBMomY->GetName());
    hPhotonBMomXBcg->Write(hPhotonBMomXBcg->GetName());
    hPhotonBMomYBcg->Write(hPhotonBMomYBcg->GetName());

    outFile->Close();
    histFile->Close();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void PlotAnaJPsi::Init(){
	//define the output file which will store all the canvases
	outFile = unique_ptr<TFile>(new TFile( outputPosition, "recreate") );

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}

    histFile = unique_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaJPsiTree, true);


    if(!tree || !bcgTree){
    	cerr << "Couldn't open tree or background tree with data. Returning." << endl;
    	return;
    }

}


void PlotAnaJPsi::invMassJPsi(int numBins, Double_t minRange, Double_t maxRange, TString nameOfSignal, TString nameOfBcg) {
    
    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("invMassJPsiCanvas"), "Fit Result", 1200, 800); 
    SetGPad();
    c->cd();
    TH1D *signal = new TH1D("invMassJPsi", "invMassJPsi", numBins, minRange, maxRange);
    TH1D *bcg = new TH1D("invMassJPsiBcg", "invMassJPsiBcg", numBins, minRange, maxRange);

    signal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    signal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    
    signal = (TH1D*)histFile->Get( nameOfSignal );
    bcg = (TH1D*)histFile->Get( nameOfBcg );

    if(!signal || !bcg || signal->GetEntries() == 0 || bcg->GetEntries() == 0 ){
        cout << "Could not load signal or background to histograms" << endl;
        return;
    }
    signal->Sumw2();
    bcg->Sumw2();

    signal->Add(bcg, -1);
    signal->SetTitle("");

    // Declare observable x
    RooRealVar x("x","m_{e^{+}e^{-}} [GeV/c^{2}]",2,4) ;

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*signal));
    RooPlot* frame = x.frame(Title("Uncorrected invariant mass of e^{+}e^{-} pairs")) ;
    dh.plotOn(frame,DataError(RooAbsData::SumW2),Name("data")); 

    //define function to fit background - polynomial 2 degree
    RooRealVar a0("a0","a0",1,-10,50); 
    RooRealVar a1("a1","a1",-1,-10,50);
    RooRealVar a2("a2","a2",0.5,-1,50);
    RooPolynomial bkg("bkg","bkg",x,RooArgList(a0,a1,a2));
  
    // define function to fit peak in data - crystal ball
    RooRealVar cbmean("mean","mean",3.1,2.8,3.2);
    RooRealVar cbsigma("sigma","sigma",5,0.01,10);
    //RooRealVar cbn("cb n","cb_n",6,0.1,10);
    RooRealVar cbn("cb n","cb_n",2,0.5,10);
    RooRealVar cbalpha("cb alpha","cb_alpha",2,0.5,10);
    RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;

    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.5,0.,1.);
    //RooAddPdf model("model","model",RooArgList(cb,bkg),bkgfrac);
    //RooAddPdf model("model","g+a",RooArgList(cb,bkg), bkgfrac) ;
    RooRealVar nsig("nsig","signal events",150,0,10000);
    RooRealVar nbkg("nbkg","signal background events",100,-100,10000);
    

    RooAddPdf model("model","model",RooArgList(cb,bkg),RooArgList(nsig,nbkg));

    //fit data with model
    RooFitResult* fitResult = model.fitTo(dh, RooFit::Save()); 


    //plotting
    model.plotOn(frame, Name("model"));
    //bkg.plotOn(frame);
    model.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed), Name("background")); 
    frame->Draw("hist E");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleOffset(0.8);
    gStyle->SetOptTitle(0);

    cerr << "Before integration" << endl;

    //integration - raw yield
    x.setRange("signal",2.8682,3.2462);
    double i_s, i_b;
    RooAbsReal *intPeakX = model.createIntegral(x,NormSet(x),Range("signal"));   //under peak of pol+gaus
    RooAbsReal *intBKGx = bkg.createIntegral(x,NormSet(x), Range("signal")) ;     //(2.908, 3.261); 
    i_s = intPeakX->getVal()*(nsig.getVal()+nbkg.getVal()); 
    i_b = intBKGx->getVal()*(nbkg.getVal()); 

    double nsig_err = nsig.getError();
    double nbkg_err = nbkg.getError();
    double yield = i_s - i_b;

    // Error propagation
    double i_s_factor = intPeakX->getVal();
    double i_b_factor = intBKGx->getVal();

    double err_yield = sqrt( pow(i_s_factor * nsig_err, 2) + pow((i_s_factor - i_b_factor) * nbkg_err, 2) );

    cout << "i_s = " << i_s << "\n";
    cout << "i_b = " << i_b << "\n";
    cout << "Yield = " << yield << " ± " << err_yield << "\n";



    //lines around mass +3-4 sigma
    auto l1 = new TLine(2.8682, 0, 2.8682, 100);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->Draw("same hist E");

    auto l2 = new TLine(3.2462, 0, 3.2462, 100);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->Draw("same hist E");

    DrawSTARpp510JPsi(0.7, 0.8, 0.9, 0.9);

    TPaveText *text = new TPaveText(0.65,0.42,0.8,0.7, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
    text->SetTextSize(0.03);
    text->SetFillColor(0);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    //text->AddText("p + p #rightarrow p #oplus J/#psi #oplus p");
    //text->AddText("Run17, #sqrt{s} = 510 GeV");
    //text->AddText("p_{T} range (0.6-1.5)");
    //text->AddText("#Delta #Phi > 1.6 && JPSIpT < 1.5"); 
    //text->AddText("After RP cuts");
    //text->AddText("This work");
    text->AddText("");
    text->AddText(TString("Without RP conditions"));
    text->AddText(Form("#mu = %.4f #pm %.4f",cbmean.getVal(),cbmean.getError()));
    text->AddText(Form("#sigma = %.4f #pm %.4f",cbsigma.getVal(),cbsigma.getError()));
    text->AddText(Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), fitResult->floatParsFinal().getSize(), frame->chiSquare()/fitResult->floatParsFinal().getSize() ) ); 
    text->AddText(Form("Raw yield = %.0f #pm %.0f", round(yield/10)*10,round(err_yield/10)*10 ));
    text->Draw("same hist E");

    cout << "The actual value and error of yield: " << nsig.getVal() << " +/- " << nsig.getError() << endl;


    TLegend *leg1 = new TLegend(0.25,0.68,0.4,0.9);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->AddEntry("data","Data", "LEP");
    leg1->AddEntry("model","Crystal Ball + Poly2","LP");
    leg1->AddEntry("background","Poly2", "LP");
    leg1->Draw("same hist E");


    outFile->cd();
    c->Write("invMassJPsi_beforeRP");

}

double PlotAnaJPsi::GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data){
    double chiSquare = frame->chiSquare();  // This calculates chi-square per degree of freedom

    int numBins = frame->GetNbinsX();  // Number of bins used in the histogram
    int nParams = model.getParameters(data)->getSize();  // Number of floating parameters in the model
    int ndf = numBins - nParams;

    return chiSquare/ndf;
}


TString PlotAnaJPsi::convertToString(double val) {

    ostringstream streamA;
    streamA << fixed << setprecision(1) << val;
    TString formattedA = streamA.str();

    return formattedA;
}

int PlotAnaJPsi::makeInt(double val) {
	int result = val;
	return result;
}



void PlotAnaJPsi::roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange){ 

    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("invMassJPsiCanvas_Misa"), "Fit Result", 1200, 800); 
    SetGPad();

    TH1D *hSignal = new TH1D("invMassJPsi_Misa", "invMassJPsi", numBins, minRange, maxRange);
    TH1D *hBcg = new TH1D("invMassJPsiBcg_Misa", "invMassJPsiBcg", numBins, minRange, maxRange);

    hSignal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    
    TString cmd, condition;
    tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
    hSignal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

    bcgTree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
    hBcg->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

    if(hSignal->GetEntries() == 0 || hBcg->GetEntries() == 0 ){
        cout << "Could not load signal or background to histograms" << endl;
        return;
    }

    hSignal->Add(hBcg, -1);
    hSignal->SetTitle("");
    // Declare observable x
    RooRealVar x("x","m_{e^{+}e^{-}} [GeV/c^{2}]",2,4) ;

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*hSignal));
    RooPlot* frame = x.frame(Title("Uncorrected invariant mass of e^{+}e^{-} pairs")) ;
    dh.plotOn(frame,DataError(RooAbsData::SumW2),Name("data")); 

    //define function to fit background - polynomial 2 degree
    RooRealVar a0("a0","a0",1,-1,20); 
    RooRealVar a1("a1","a1",-1,-10,10);
    RooRealVar a2("a2","a2",0.5,-1,1);
    RooPolynomial bkg("bkg","bkg",x,RooArgList(a0,a1,a2));
  
    // define function to fit peak in data - crystal ball
    RooRealVar cbmean("mean","mean",3.1,2.8,3.2);
    RooRealVar cbsigma("sigma","sigma",5,0.01,10);
    //RooRealVar cbn("cb n","cb_n",6,0.1,10);
    RooRealVar cbn("cb n","cb_n",2,1,10);
    RooRealVar cbalpha("cb alpha","cb_alpha",2,1,10);
    RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;



    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.5,0.,1.);
    //RooAddPdf model("model","model",RooArgList(cb,bkg),bkgfrac);
    //RooAddPdf model("model","g+a",RooArgList(cb,bkg), bkgfrac) ;
    RooRealVar nsig("nsig","signal events",150,0,10000);
    RooRealVar nbkg("nbkg","signal background events",100,-10000,10000);
    

    RooAddPdf model("model","model",RooArgList(cb,bkg),RooArgList(nsig,nbkg));

    //fit data with model
    RooFitResult* fitResult = model.fitTo(dh, RooFit::Save()); 

    //plotting
    model.plotOn(frame, Name("model"));
    //bkg.plotOn(frame);
    model.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed), Name("background")); 
    frame->Draw("hist E");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleOffset(0.8);
    gStyle->SetOptTitle(0);

    cerr << "Before integration" << endl;

    //integration - raw yield
    x.setRange("signal",2.8682,3.2462);
    double i_s, i_b;
    RooAbsReal *intPeakX = model.createIntegral(x,NormSet(x),Range("signal"));   //under peak of pol+gaus
    RooAbsReal *intBKGx = bkg.createIntegral(x,NormSet(x), Range("signal")) ;     //(2.908, 3.261); 
    i_s = intPeakX->getVal()*(nsig.getVal()+nbkg.getVal()); 
    i_b = intBKGx->getVal()*(nbkg.getVal()); 

    double nsig_err = nsig.getError();
    double nbkg_err = nbkg.getError();
    double yield = i_s - i_b;

    // Error propagation
    double i_s_factor = intPeakX->getVal();
    double i_b_factor = intBKGx->getVal();

    double err_yield = sqrt( pow(i_s_factor * nsig_err, 2) + pow((i_s_factor - i_b_factor) * nbkg_err, 2) );

    cout << "i_s = " << i_s << "\n";
    cout << "i_b = " << i_b << "\n";
    cout << "Yield = " << yield << " ± " << err_yield << "\n";



    //lines around mass +3-4 sigma
    auto l1 = new TLine(2.8682, -5, 2.8682, 17);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->Draw("same hist E");

    auto l2 = new TLine(3.2462, -5, 3.2462, 17);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->Draw("same hist E");

    DrawSTARpp510JPsi(0.7, 0.8, 0.9, 0.9);

    TPaveText *text = new TPaveText(0.65,0.42,0.8,0.7, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
    text->SetTextSize(0.03);
    text->SetFillColor(0);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    //text->AddText("p + p #rightarrow p #oplus J/#psi #oplus p");
    //text->AddText("Run17, #sqrt{s} = 510 GeV");
    //text->AddText("p_{T} range (0.6-1.5)");
    //text->AddText("#Delta #Phi > 1.6 && JPSIpT < 1.5"); 
    //text->AddText("After RP cuts");
    //text->AddText("This work");
    text->AddText("");
    text->AddText(TString("With RP conditions"));
    text->AddText(Form("#mu = %.4f #pm %.4f",cbmean.getVal(),cbmean.getError()));
    text->AddText(Form("#sigma = %.4f #pm %.4f",cbsigma.getVal(),cbsigma.getError()));
    text->AddText(Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), fitResult->floatParsFinal().getSize(), frame->chiSquare()/fitResult->floatParsFinal().getSize() ) ); 
    text->AddText(Form("Raw yield = %.0f #pm %.0f", round(yield/10)*10, round( (err_yield + 10)/10)*10 ));
    text->Draw("same hist E");



    TLegend *leg1 = new TLegend(0.25,0.68,0.4,0.9);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->AddEntry("data","Data", "LEP");
    leg1->AddEntry("model","Crystal Ball + Poly2","LP");
    leg1->AddEntry("background","Poly2", "LP");
    leg1->Draw("same hist E");


    outFile->cd();
    c->Write("invMassJPsi_afterRP");
}

/*
void PlotAnaJPsi::processTreeForPt(TString nameOfTree){

    int nEntries = tree->GetEntries();
    cout << "Total entries in tree: " << nEntries << std::endl;
    
    // Loop over all events
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        





    }
}
*/
