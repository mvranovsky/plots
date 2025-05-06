#include "PlotAnaJPsi.h"

PlotAnaJPsi::PlotAnaJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotAnaJPsi::Make(){


	int nBins = 40;
    double lowerLim = 2;
    double upperLim = 4;
    roofit_invMassFit_Misa(nBins, lowerLim, upperLim, true,false);  // true = subtract background, false = do not fit background with poly2
    cout << "Finished Misa's plot" << endl;

    nSigmaCorrPlot(1, false);  // 1 = electron, 2 = pion, 3 = kaon, 4 = proton
    
    // plots for reconstruction efficiency

    if(strstr(outputPosition,"EmbeddingJPsi")){
        reconstructionEfficiency(1, "slight_EtaCut1_100k.root", "recoEff_pairRap");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(2, "slight_EtaCut1_100k.root", "recoEff_eta");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(3, "slight_EtaCut1_100k.root", "recoEff_phi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(4, "slight_EtaCut1_100k.root", "recoEff_pTJPsi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(5, "slight_EtaCut1_100k.root", "recoEff_pT");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    }else{
        cout << "This is not embedding, this is data. Skipping reconstruction efficiency plots." << endl;
    }
    
    
    controlPlotsComparison(true); // true = JPsi, false = all
    
    
    // save all the histograms to canvases into outfile
    handleHistograms();
    

    /*
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
    */


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

void PlotAnaJPsi::JPsiControlPlots(){

    tree->Draw("phiHadron0 - phiHadron1>>hist(24, -3.14, 3.14)", "invMass > 2.8682 && invMass < 3.2462");
    TH1* phiDelta = (TH1*)gPad->GetPrimitive("hist");
    if(phiDelta && phiDelta->GetEntries() > 0){
        phiDelta->SetTitle("hDeltaPhi");
        phiDelta->GetXaxis()->SetTitle("#Delta #phi [rad]");
        phiDelta->GetYaxis()->SetTitle("counts");
        TH1General( TString("deltaPhiDaughters"), phiDelta);
    }

    tree->Draw("|etaHadron0 - etaBemc0|>>hist(40, 0,0.5)");
    tree->Draw("|etaHadron1 - etaBemc1|+>>hist");
    TH1* etaDeltaTpcBemc = (TH1*)gPad->GetPrimitive("hist");
    if(etaDeltaTpcBemc && etaDeltaTpcBemc->GetEntries() >0){
        etaDeltaTpcBemc->SetTitle("hDeltaEta");
        etaDeltaTpcBemc->GetXaxis()->SetTitle("#Delta #eta [-]");
        etaDeltaTpcBemc->GetYaxis()->SetTitle("counts");
        TH1General(TString("etaDeltaTB"),etaDeltaTpcBemc);
    }

    /*
    bcgTree->Draw("|etaHadron0 - etaBemc0|>>hist(40, 0,0.5)");
    bcgTree->Draw("|etaHadron1 - etaBemc1|+>>hist");
    TH1* etaDeltaTpcBemcBcg = (TH1*)gPad->GetPrimitive("hist");
    if(etaDeltaTpcBemcBcg && etaDeltaTpcBemcBcg->GetEntries() >0){
        etaDelta
        TH1General(TString("etaDeltaTBBcg"),etaDeltaTpcBemcBcg);
    }

    
    // pridat korelacny plot pT pre J/Psi a pT pre proton v RP
    tree->Draw("pt:(ptRpEast + ptRpWest)>>hist(20,0,2,20,0,2)");
    TH2* pTcorr = (TH2*)gPad->GetPrimitive("hist");
    if(pTcorr && pTcorr->GetEntries() >0){
        TH2General(TString("pTJPsiRPPcorrPlot"), pTcorr);
    }
    */


}


bool PlotAnaJPsi::plot2Dists(TH1 *hData, TH1* hEmb, TString outName){


    if(!hData || !hEmb){
        cout << "Could not get histograms. Returning." << endl;
        return false;
    }

    CreateCanvas(&canvas, outName, 1200, 800);
    SetGPad();
    SetHistStyle(hData, kBlack, markerStyleTypical);
    SetHistStyle(hEmb, kRed, markerStyleTypical+1);
    
    hEmb->Draw();
    hData->Draw("same");
    DrawSTARpp510JPsi(0.6,0.85,0.93,0.93, 0.01);

    CreateLegend(&legend, 0.15, 0.8, 0.35, 0.9);
    legend->AddEntry(hData, "Data", "LEP");
    legend->AddEntry(hEmb, "Embedding", "LEP");
    legend->Draw("same");

    outFile->cd();
    canvas->Write(outName);
    canvas->Close();
    return true;


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


void PlotAnaJPsi::nSigmaCorrPlot(int particles, bool justJPsi){
    // 1 = electron, 2 = pion, 3 = kaon, 4 = proton
    TString condition = "";
    if(justJPsi){
        condition = Form("invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
    }

    if(particles == 1){
        tree->Draw("nSigmaTPCelectronMinus:nSigmaTPCelectronPlus>>histnsig(40,-4,4,40,-4,4)", condition);
        tree->Draw("chiSquareelectron>>chiPlot(50,0,50)", condition );
    }else if(particles == 2){
        tree->Draw("nSigmaTPCpionPlus:nSigmaTPCpionMinus>>histnsig(40,-4,4,40,-4,4)", condition );
        tree->Draw("chiSquarekion>>chiPlot(50,0,50)", condition );
    }else if(particles == 3){
        tree->Draw("nSigmaTPCkaonPlus:nSigmaTPCkaonMinus>>histnsig(40,-4,4,40,-4,4)", condition );
        tree->Draw("chiSquarekaon>>chiPlot(50,0,50)", condition );
    }else if(particles == 4){
        tree->Draw("nSigmaTPCprotonPlus:nSigmaTPCprotonMinus>>histnsig(40,-4,4,40,-4,4)", condition );
        tree->Draw("chiSquareproton>>chiPlot(50,0,50)", condition );
    }else{
        cout << "Unknown particle type. Returning." << endl;
        return;
    }
    TH2* nSigma = dynamic_cast<TH2*>(gDirectory->FindObject("histnsig"));
    if(nSigma && nSigma->GetEntries() > 0){
        nSigma->GetXaxis()->SetTitle("n#sigma_{e} [-]");
        nSigma->GetYaxis()->SetTitle("n#sigma_{e} [-]");
        gStyle->SetOptStat("euo"); // Enable underflow and overflow bins
        TH2General(TString("nSigmaJPsi"), nSigma);
        cout << "Obtained and saved nSigma histogram." << endl;
    }else{
        cout << "Could not get nSigma histogram." << endl;
    }
    
    TH1* chiSquare = dynamic_cast<TH1*>(gDirectory->FindObject("chiPlot"));
    if(chiSquare && chiSquare->GetEntries() > 0){
        chiSquare->GetXaxis()->SetTitle("#chi^{2}_{ee} [-]");
        chiSquare->GetYaxis()->SetTitle("counts");
        TH1General(TString("chiSquareJPsi"), chiSquare);
        cout << "Obtained and saved chiSquare histogram." << endl;
    }else{
        cout << "Could not get chiSquare histogram." << endl;
    }
    
}



void PlotAnaJPsi::roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange,bool subtractBcg ,bool noBcgFit){ 

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
    if(subtractBcg){
        hSignal->Add(hBcg, -1);
    }
    // Declare observable x
    RooRealVar x("x","m_{e^{+}e^{-}} [GeV/c^{2}]",2,4) ;

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*hSignal));
    RooPlot* frame = x.frame(Title("")) ;
    dh.plotOn(frame,DataError(RooAbsData::SumW2),Name("data")); 

    //define function to fit background - polynomial 2 degree
    RooRealVar a0("a0","a0",1,-20,20); 
    RooRealVar a1("a1","a1",-1,-20,20);
    RooRealVar a2("a2","a2",0.5,-20,20);
    RooPolynomial bkg("bkg","bkg",x,RooArgList(a0,a1,a2));
  
    // define function to fit peak in data - crystal ball
    RooRealVar cbmean("mean","mean",3.0908,2.8,3.2);
    RooRealVar cbsigma("sigma","sigma",0.0543,0.0,0.065);
    RooRealVar cbn("cb n","cb_n",1.41,1.410,1.411);
    RooRealVar cbalpha("cb alpha","cb_alpha",1.94,1.940,1.941);
    //RooRealVar cbn("cb n","cb_n",2,0.1,10);
    //RooRealVar cbalpha("cb alpha","cb_alpha",2,0,10);
    RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;
    

    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.5,0.,1.);
    //RooAddPdf model("model","model",RooArgList(cb,bkg),bkgfrac);
    //RooAddPdf model("model","g+a",RooArgList(cb,bkg), bkgfrac) ;
    RooRealVar nsig("nsig","signal events",150,0,2000);
    RooRealVar nbkg("nbkg","signal background events",100,-10000,10000);
    

    RooAddPdf* model;
    if(noBcgFit){
        model = new RooAddPdf("model","model",RooArgList(cb),RooArgList(nsig));
    }else{
        model = new RooAddPdf("model","model",RooArgList(cb,bkg),RooArgList(nsig,nbkg));
    }
    //fit data with model
    RooFitResult* fitResult = model->fitTo(dh, RooFit::Save()); 

    //plotting
    model->plotOn(frame, Name("model"));
    //bkg.plotOn(frame);
    if(!noBcgFit){
        model->plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed), Name("background")); 
    }
    frame->Draw("hist E");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleOffset(0.8);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    // set the title of the canvas to ""
    c->SetTitle("");
    // set the title of RooFit frame to ""
    frame->SetTitle("");

    lowLimInvMass = cbmean.getVal() - 4*cbsigma.getVal();
    topLimInvMass = cbmean.getVal() + 3*cbsigma.getVal();
    //integration - raw yield
    x.setRange("signal",lowLimInvMass,topLimInvMass);
    double i_s, i_b;
    RooAbsReal *intPeakX = model->createIntegral(x,NormSet(x),Range("signal"));   //under peak of pol+gaus
    double nsig_err = nsig.getError();
    double nbkg_err = 0;
    double i_b_factor = 0;
    double i_s_factor = intPeakX->getVal();
    
    if(!noBcgFit){
        RooAbsReal *intBKGx = bkg.createIntegral(x,NormSet(x), Range("signal")) ;     //(2.908, 3.261); 
        nbkg_err = nbkg.getError(); 
        i_b_factor = intBKGx->getVal();
        i_b = intBKGx->getVal()*(nbkg.getVal());
        i_s = intPeakX->getVal()*(nsig.getVal()+nbkg.getVal());
    }else{
        i_b = 0;
        i_s = intPeakX->getVal()*(nsig.getVal());
    }

    double yield = i_s - i_b;
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
    if(!noBcgFit){
        leg1->AddEntry("background","Poly2", "LP");
    }
    leg1->Draw("same hist E");


    outFile->cd();
    c->Write("invMassJPsi");
}



void PlotAnaJPsi::reconstructionEfficiency(int SWITCH, TString nameOfStarlightFile, TString nameOfOutput){  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    // open tree from starlight


    TFile *starlightFile = new TFile(nameOfStarlightFile, "READ");
    if(!starlightFile || starlightFile->IsZombie() || !starlightFile->IsOpen()){
        cerr << "Could not get file with starlight data. Returning." << endl;
        return;
    }
    TTree *starlightTree = (TTree*)starlightFile->Get("T");
    if(!starlightTree){
        cerr << "Could not get starlight tree. Returning." << endl;
        return;
    }
    
    
    TH1 *h1, *h2;
    TString xAxisDescription;
    if(SWITCH == 1){ // pair rapidity
        tree->Draw("pairRapidity>>hist1(20,-1,1)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw("rapVM>>hist2(20,-1,1)");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "y_{J/#psi} [-]";
    }else if(SWITCH == 2){ // eta of daughters
        tree->Draw("etaHadron0>>hist1(40,-1,1)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("etaHadron1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw("ed1>>hist2(40,-1,1)");
        starlightTree->Draw("ed2>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#eta_{e} [-]";
    }else if(SWITCH == 3){ // phi of daughters
        tree->Draw("phiHadron0>>hist1(40,-3.14,3.14)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("phiHadron1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw("(phiatemchd1 - 3.14)>>hist2(40,-3.14,3.14)");
        starlightTree->Draw("(phiatemchd2 - 3.14)>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#phi [rad]";
        
    }else if(SWITCH == 4){ // pT of pair
        tree->Draw("pt>>hist1(40,0,2)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw("ptVM>>hist2(40,0,2)");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        
        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
        
    }else if(SWITCH == 5){ // pT of daughters
        tree->Draw("pTInGev0>>hist1(50,0.5,3)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("pTInGev1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw("ptd1>>hist2(50,0.5,3)");
        starlightTree->Draw("ptd2>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "p_{T}^{e} [GeV/c]";
    }else{
        cout << "Unknown switch. Returning." << endl;
        return;
    }
    if(h1 && h2 && h1->GetEntries() > 0 && h2->GetEntries() > 0){

        CreateCanvas(&canvas, nameOfOutput, 1200, 800);
        SetGPad();
        h1->Sumw2();
        h2->Sumw2();
        // instead of dividing the graphs, se TGraphAssymErrors and calculate the efficiency of reconstruction
        TEfficiency *eff = new TEfficiency(*h2, *h1);
        TGraphAsymmErrors *g = eff->CreateGraph();
        g->SetMarkerColor(kRed);
        g->SetName(nameOfOutput);
        g->SetTitle(nameOfOutput);
        g->GetXaxis()->SetTitle(xAxisDescription);
        g->GetYaxis()->SetTitle("reconstruction efficiency");
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kRed);
        g->SetLineColor(kRed); 
        g->GetYaxis()->SetRangeUser(0,1);
        g->Draw("AEP");

        DrawSTARpp510JPsi(0.7, 0.8, 0.9, 0.9);

        outFile->cd();
        canvas->Write(nameOfOutput);

    }
    else{
        cout << "Could not get histograms." << endl;
    }
    // close the starlight file
    starlightFile->Close();

    return;
}


void PlotAnaJPsi::controlPlotsComparison(bool justJPsi){

    if(!gPad){
        cerr << "gPad is not initialized. Returning." << endl;
        return;
    }

    TString condition = "";
    if(justJPsi){
        condition = Form("invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
    }

    TFile *controlFile;
    controlFile = new TFile("controlPlots.root", "RECREATE");

    tree->Draw("etaHadron0>>hist1(40,-1.5,1.5)", condition);
    tree->Draw("etaHadron1>>+hist1", condition);
    TH1D* etaHadron = (TH1D*)gPad->GetPrimitive("hist1");

    tree->Draw("phiHadron0>>hist2(40,-3.14,3.14)", condition);
    tree->Draw("phiHadron1>>+hist2", condition);
    TH1D* phiHadron = (TH1D*)gPad->GetPrimitive("hist2");

    tree->Draw("pTInGev0>>hist3(40,0.5,2.5)", condition);
    tree->Draw("pTInGev1>>+hist3", condition);
    TH1D* pTInGev = (TH1D*)gPad->GetPrimitive("hist3");

    
    tree->Draw("vertexZInCm>>hist4(200,-100,100)", condition);
    TH1D* vertexZ = (TH1D*)gPad->GetPrimitive("hist4");
    
    tree->Draw("dcaZInCm0>>hist5(30,-1.5,1.5)", condition);
    tree->Draw("dcaZInCm1>>+hist5", condition);
    TH1D* dcaZ = (TH1D*)gPad->GetPrimitive("hist5");
    
    tree->Draw("dcaXYInCm0>>hist6(20,0,2)", condition);
    tree->Draw("dcaXYInCm1>>+hist6", condition);
    TH1D* dcaXY = (TH1D*)gPad->GetPrimitive("hist6");
    

    tree->Draw("nHitsFit0>>hist7(50,0,50)", condition);
    tree->Draw("nHitsFit1>>+hist7", condition);
    TH1D* nHitsFit = (TH1D*)gPad->GetPrimitive("hist7");

    tree->Draw("nHitsDEdx0>>hist8(50,0,50)", condition);
    tree->Draw("nHitsDEdx1>>+hist8", condition);
    TH1D* nHitsDEdx = (TH1D*)gPad->GetPrimitive("hist8");
    

    tree->Draw("pt>>hist9(60,0,1.5)", condition);
    TH1D* pt = (TH1D*)gPad->GetPrimitive("hist9");

    controlFile->cd();
    if(etaHadron && etaHadron->GetEntries() > 0){
        etaHadron->Write("hEta");
    }else{
        cout << "Could not get etaHadron histogram." << endl;
    }

    if(phiHadron && phiHadron->GetEntries() > 0){
        phiHadron->Write("hPhi");
    }else{
        cout << "Could not get phiHadron histogram." << endl;
    }

    if(pTInGev && pTInGev->GetEntries() > 0){
        pTInGev->Write("hPt");
    }else{
        cout << "Could not get pTInGev histogram." << endl;
    }

    
    if(vertexZ && vertexZ->GetEntries() > 0){
        vertexZ->Write("hVertexZ");
    }else{
        cout << "Could not get vertexZ histogram." << endl;
    }
    
    if(dcaZ && dcaZ->GetEntries() > 0){
        dcaZ->Write("hDCAZ");
    }else{
        cout << "Could not get dcaZ histogram." << endl;
    }

    if(dcaXY && dcaXY->GetEntries() > 0){
        dcaXY->Write("hDCAXY");
    }else{
        cout << "Could not get dcaXY histogram." << endl;
    }
    

    if(nHitsFit && nHitsFit->GetEntries() > 0){
        nHitsFit->Write("hNHitsFit");
    }else{
        cout << "Could not get nHitsFit histogram." << endl;
    }

    if(nHitsDEdx && nHitsDEdx->GetEntries() > 0){
        nHitsDEdx->Write("hNHitsDEdx");
    }else{
        cout << "Could not get nHitsDEdx histogram." << endl;
    }

    if(pt && pt->GetEntries() > 0){
        pt->Write("hPtJPsi");
    }else{
        cout << "Could not get pt histogram." << endl;
    }

    cout << "All plots from data have been saved to controlPlots.root file" << endl;
    
    controlFile->Close();
}


