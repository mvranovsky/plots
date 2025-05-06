#include "PlotEmbeddingJPsi.h"

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotEmbeddingJPsi::Make(){



	int nBins = 40;
    double lowerLim = 2;
    double upperLim = 4;
    roofit_invMassFit_Misa(nBins, lowerLim, upperLim, false, true);  // true if fitting without background
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
    
    //JPsiControlPlots();
    
    controlPlotsComparison(true); // true = JPsi, false = all
    
    // save all the histograms to canvases into outfile
    handleHistograms();
    

    outFile->Close();
    histFile->Close();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void PlotEmbeddingJPsi::Init(){
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
	ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, true);


    if(!tree || !bcgTree){
    	cerr << "Couldn't open tree or background tree with data. Returning." << endl;
    	return;
    }

}

void PlotEmbeddingJPsi::JPsiControlPlots(){

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


void PlotEmbeddingJPsi::controlPlotsComparison(bool justJPsi){

    if(!gPad){
        cerr << "gPad is not initialized. Returning." << endl;
        return;
    }

    TString condition = "";
    if(justJPsi){
        condition = Form("invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
    }

    TFile *controlFile;

    controlFile = new TFile("controlPlots.root", "read");
    if(!controlFile || controlFile->IsZombie() || !controlFile->IsOpen()){
        cerr << "Could not open control plots file. Returning." << endl;
        return;
    }

    TH1* hEtaData = dynamic_cast<TH1*>(controlFile->Get("hEta"));
    tree->Draw("etaHadron0>>hist11(40,-1.5,1.5)", condition);
    tree->Draw("etaHadron1>>+hist11", condition);
    TH1* hEtaEmb = dynamic_cast<TH1*>(gPad->FindObject("hist11"));
    if(hEtaData && hEtaEmb){
        hEtaData->GetXaxis()->SetTitle("#eta [-]");
        hEtaEmb->GetXaxis()->SetTitle("#eta [-]");
        hEtaData->GetYaxis()->SetTitle("counts");
        hEtaEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hEtaData, hEtaEmb, "hEtaComparison");
    }else{
        cout << "Could not get eta histograms for comparison." << endl;
    }

    TH1* hPhiData = dynamic_cast<TH1*>(controlFile->Get("hPhi"));
    tree->Draw("phiHadron0>>hist12(40,-3.14,3.14)", condition);
    tree->Draw("phiHadron1>>+hist12", condition);
    TH1* hPhiEmb = dynamic_cast<TH1*>(gPad->FindObject("hist12"));
    if(hPhiData && hPhiEmb){
        hPhiData->GetXaxis()->SetTitle("#phi [rad]");
        hPhiEmb->GetXaxis()->SetTitle("#phi [rad]");
        hPhiData->GetYaxis()->SetTitle("counts");
        hPhiEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hPhiData, hPhiEmb, "hPhiComparison");
    }else{
        cout << "Could not get phi histograms for comparison." << endl;
    }

    TH1* hPtData = dynamic_cast<TH1*>(controlFile->Get("hPt"));
    tree->Draw("pTInGev0>>hist13(40,0.5,2.5)", condition);
    tree->Draw("pTInGev1>>+hist13", condition);
    TH1* hPtEmb = dynamic_cast<TH1*>(gPad->FindObject("hist13"));
    if(hPtData && hPtEmb){
        hPtData->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        hPtEmb->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        hPtData->GetYaxis()->SetTitle("counts");
        hPtEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hPtData, hPtEmb, "hPtComparison");
    }else{
        cout << "Could not get pT histograms for comparison." << endl;
    }
    
    TH1* hVertexZData = dynamic_cast<TH1*>(controlFile->Get("hVertexZ"));
    tree->Draw("vertexZInCm>>hist14(200,-100,100)", condition);
    TH1* hVertexZEmb = dynamic_cast<TH1*>(gPad->FindObject("hist14"));
    if(hVertexZData && hVertexZEmb){
        hVertexZData->GetXaxis()->SetTitle("z_{vertex} [cm]");
        hVertexZEmb->GetXaxis()->SetTitle("z_{vertex} [cm]");
        hVertexZData->GetYaxis()->SetTitle("counts");
        hVertexZEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hVertexZData, hVertexZEmb, "hVertexZComparison");
    }else{
        cout << "Could not get vertexZ histograms for comparison." << endl;
    }
    TH1 *hDCAZData = dynamic_cast<TH1*>(controlFile->Get("hDCAZ"));
    tree->Draw("dcaZInCm0>>hist15(30,-1.5,1.5)", condition);
    tree->Draw("dcaZInCm1>>+hist15", condition);
    TH1 *hDCAZEmb = dynamic_cast<TH1*>(gPad->FindObject("hist15"));
    if(hDCAZData && hDCAZEmb){
        hDCAZData->GetXaxis()->SetTitle("DCA_{z} [cm]");
        hDCAZEmb->GetXaxis()->SetTitle("DCA_{z} [cm]");
        hDCAZData->GetYaxis()->SetTitle("counts");
        hDCAZEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hDCAZData, hDCAZEmb, "hDCAZComparison");
    }else{
        cout << "Could not get dcaZ histograms for comparison." << endl;
        if(!hDCAZData){
            cout << "hDCAZData is null." << endl;
        }
        if(!hDCAZEmb){
            cout << "hDCAZEmb is null." << endl;
        }
    }
    TH1* hDCAXYData = dynamic_cast<TH1*>(controlFile->Get("hDCAXY"));
    tree->Draw("dcaXYInCm0>>hist16(20,0,2)", condition);
    tree->Draw("dcaXYInCm1>>+hist16", condition);
    TH1* hDCAXYEmb = dynamic_cast<TH1*>(gPad->FindObject("hist16"));
    if(hDCAXYData && hDCAXYEmb){
        hDCAXYData->GetXaxis()->SetTitle("DCA_{xy} [cm]");
        hDCAXYEmb->GetXaxis()->SetTitle("DCA_{xy} [cm]");
        hDCAXYData->GetYaxis()->SetTitle("counts");
        hDCAXYEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hDCAXYData, hDCAXYEmb, "hDCAXYComparison");
    }else{
        cout << "Could not get dcaXY histograms for comparison." << endl;
    }
    

    TH1* hNHitsFitData = dynamic_cast<TH1*>(controlFile->Get("hNHitsFit"));
    tree->Draw("nHitsFit0>>hist17(50,0,50)", condition);
    tree->Draw("nHitsFit1>>+hist17", condition);
    TH1* hNHitsFitEmb = dynamic_cast<TH1*>(gPad->FindObject("hist17"));
    if(hNHitsFitData && hNHitsFitEmb){
        hNHitsFitData->GetXaxis()->SetTitle("N_{Hits}^{Fit}");
        hNHitsFitEmb->GetXaxis()->SetTitle("N_{Hits}^{Fit}");
        hNHitsFitData->GetYaxis()->SetTitle("counts");
        hNHitsFitEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hNHitsFitData, hNHitsFitEmb, "hNHitsFitComparison");
    }else{
        cout << "Could not get nHitsFit histograms for comparison." << endl;
    }

    TH1* hNHitsDEdxData = dynamic_cast<TH1*>(controlFile->Get("hNHitsDEdx"));
    tree->Draw("nHitsDEdx0>>hist18(50,0,50)", condition);
    tree->Draw("nHitsDEdx1>>+hist18", condition);
    TH1* hNHitsDEdxEmb = dynamic_cast<TH1*>(gPad->FindObject("hist18"));
    if(hNHitsDEdxData && hNHitsDEdxEmb){
        hNHitsDEdxData->GetXaxis()->SetTitle("N_{Hits}^{dE/dx}");
        hNHitsDEdxEmb->GetXaxis()->SetTitle("N_{Hits}^{dE/dx}");
        hNHitsDEdxData->GetYaxis()->SetTitle("counts");
        hNHitsDEdxEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hNHitsDEdxData, hNHitsDEdxEmb, "hNHitsDEdxComparison");
    }else{
        cout << "Could not get nHitsDEdx histograms for comparison." << endl;
    }


    TH1* pTJPsiData = dynamic_cast<TH1*>(controlFile->Get("hPtJPsi"));
    tree->Draw("pt>>hist19(60,0,1.5)", condition);
    TH1* pTJPsiEmb = dynamic_cast<TH1*>(gPad->FindObject("hist19"));
    if(pTJPsiData && pTJPsiEmb){
        pTJPsiData->GetXaxis()->SetTitle("p_{T}^{J/#psi} [GeV/c]");
        pTJPsiEmb->GetXaxis()->SetTitle("p_{T}^{J/#psi} [GeV/c]");
        pTJPsiData->GetYaxis()->SetTitle("counts");
        pTJPsiEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(pTJPsiData, pTJPsiEmb, "hPtJPsiComparison");
    }else{
        cout << "Could not get pTJPsi histograms for comparison." << endl;
    }

    cout << "All plots from embedding have been saved to output file. " << endl;

    
    controlFile->Close();
}



bool PlotEmbeddingJPsi::plot2Dists(TH1 *hData, TH1* hEmb, TString outName){


    if(!hData || !hEmb){
        cout << "Could not get histograms. Returning." << endl;
        return false;
    }

    // rescale embedding histogram to match data histogram
    hEmb->Scale(hData->GetEntries()/hEmb->GetEntries());

    CreateCanvas(&canvas, outName, 1200, 800);
    SetGPad();
    SetHistStyle(hData, kBlack, markerStyleTypical);
    SetHistStyle(hEmb, kRed, markerStyleTypical+1);
    
    hEmb->GetYaxis()->SetRangeUser(0, hEmb->GetMaximum()*1.2);
    hEmb->Draw();
    hData->Draw("same");
    DrawSTARpp510JPsi(0.6,0.85,0.93,0.93, 0.01);

    CreateLegend(&legend, 0.2, 0.8, 0.35, 0.88);
    legend->SetTextSize(0.04);
    legend->AddEntry(hData, "Data", "LEP");
    legend->AddEntry(hEmb, "Embedding (rescaled)", "LEP");
    legend->Draw("same");

    outFile->cd();
    canvas->Write(outName);
    canvas->Close();
    return true;


}



double PlotEmbeddingJPsi::GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data){
    double chiSquare = frame->chiSquare();  // This calculates chi-square per degree of freedom

    int numBins = frame->GetNbinsX();  // Number of bins used in the histogram
    int nParams = model.getParameters(data)->getSize();  // Number of floating parameters in the model
    int ndf = numBins - nParams;

    return chiSquare/ndf;
}


TString PlotEmbeddingJPsi::convertToString(double val) {

    ostringstream streamA;
    streamA << fixed << setprecision(1) << val;
    TString formattedA = streamA.str();

    return formattedA;
}

int PlotEmbeddingJPsi::makeInt(double val) {
	int result = val;
	return result;
}


void PlotEmbeddingJPsi::nSigmaCorrPlot(int particles, bool justJPsi){
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



void PlotEmbeddingJPsi::roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange, bool noBcg, bool plotDataPeak){ 

    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("invMassJPsiCanvas_Misa"), "Fit Result", 1200, 800); 
    SetGPad();

    TH1D *hSignal = new TH1D("invMassJPsi_Misa", "invMassJPsi", numBins, minRange, maxRange);

    hSignal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    
    TString cmd, condition;
    tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
    hSignal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

    if(!hSignal || hSignal->GetEntries() == 0 ){
        cout << "Could not load signal or background to histograms" << endl;
        return;
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
    RooRealVar cbmean("mean","mean",3.1,2.8,3.2);
    RooRealVar cbsigma("sigma","sigma",0.05,0.01,0.1);
    //RooRealVar cbn("cb n","cb_n",6,0.1,10);
    RooRealVar cbn("cb n","cb_n",2,1,10);
    RooRealVar cbalpha("cb alpha","cb_alpha",2,1,10);
    RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;
    
    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.5,0.,1.);
    //RooAddPdf model("model","model",RooArgList(cb,bkg),bkgfrac);
    //RooAddPdf model("model","g+a",RooArgList(cb,bkg), bkgfrac) ;
    RooRealVar nsig("nsig","signal events",150,0,15000);
    RooRealVar nbkg("nbkg","signal background events",100,-10000,10000);
    
    RooRealVar dataCMMean("dataCMMean","dataCMMean",3.097,2.8,3.2);
    RooRealVar dataCMSigma("dataCMSigma","dataCMSigma",0.05,0.01,0.1);
    RooRealVar dataCBAlpha("dataCBAlpha","dataCBAlpha",4.662,1,10);
    RooRealVar dataCBN("dataCBN","dataCBN",2.876,1,10);
    RooCBShape dataCB("dataCB","dataCB",x,dataCMMean,dataCMSigma,dataCBAlpha,dataCBN);
    
    RooRealVar nsigdata("nsigdata","signal events",10823,0,15000);
    RooAddPdf *dataModel = new RooAddPdf("dataModel","dataModel",RooArgList(dataCB),RooArgList(nsigdata));
    
    RooAddPdf* model;
    if(noBcg){
        model = new RooAddPdf("model","model",RooArgList(cb),RooArgList(nsig));
    }else{
        model = new RooAddPdf("model","model",RooArgList(cb,bkg),RooArgList(nsig,nbkg));
    }
    //fit data with model
    RooFitResult* fitResult = model->fitTo(dh, RooFit::Save()); 
    
    if(plotDataPeak){
        dataModel->plotOn(frame,LineStyle(kDashed), LineColor(kGreen), Name("RescaledData"));
    }
    //plotting
    model->plotOn(frame, Name("model"));
    //bkg.plotOn(frame);
    if(!noBcg){
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
    
    if(!noBcg){
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
    cout << "nsig = " << nsig.getVal() << " ± " << nsig.getError() << "\n";
    
    cout << "yield = "<< makeInt(yield) << " ± " << makeInt(err_yield) << "\n";


    //lines around mass +3-4 sigma
    auto l1 = new TLine(2.8682, 0, 2.8682, 500);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->Draw("same hist E");

    auto l2 = new TLine(3.2462, 0, 3.2462, 500);
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
    //text->AddText(TString("With RP conditions"));
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
    if(!noBcg){
        leg1->AddEntry("background","Poly2", "LP");
    }

    if(plotDataPeak){
        leg1->AddEntry("RescaledData","Rescaled Data", "LP");
    }

    leg1->Draw("same hist E");


    outFile->cd();
    c->Write("invMassJPsi");
    c->SaveAs("invMassJPsi.pdf");
}



void PlotEmbeddingJPsi::reconstructionEfficiency(int SWITCH, TString nameOfStarlightFile, TString nameOfOutput){  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
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
        
        starlightTree->Draw("etad1>>hist2(40,-1,1)");
        starlightTree->Draw("etad2>>+hist2");
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
        tree->Draw("pt>>hist1(30,0,1.5)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw("ptVM>>hist2(30,0,1.5)");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        
        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
        
    }else if(SWITCH == 5){ // pT of daughters
        tree->Draw("pTInGev0>>hist1(40,0.5,2.5)", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("pTInGev1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw("ptd1>>hist2(40,0.5,2.5)");
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

        cout << "Reconstruction efficiency for " << nameOfOutput << ": " << h1->GetEntries() << " / " << h2->GetEntries() << endl;
        h1->SetBinContent(0, 0);
        h1->SetBinContent(h1->GetNbinsX()+1, 0);


        TEfficiency *eff = new TEfficiency(*h1, *h2);
        TGraphAsymmErrors *g = eff->CreateGraph();
        g->SetMarkerColor(kBlue);
        g->SetName(nameOfOutput);
        g->SetTitle(nameOfOutput);
        g->GetXaxis()->SetTitle(xAxisDescription);
        g->GetYaxis()->SetTitle("reconstruction efficiency");
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlue);
        g->SetLineColor(kBlue); 
        g->GetYaxis()->SetRangeUser(0,1);
        g->Draw("AEP");

        DrawSTARpp510JPsi(0.7, 0.8, 0.9, 0.9);


        // draw a TLine at overall efficiency
        double overallEfficiency = h1->GetEntries()/h2->GetEntries();
        TLine *line = new TLine(h1->GetXaxis()->GetXmin() , overallEfficiency, h1->GetXaxis()->GetXmax(), overallEfficiency);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        // create a legend
        TLegend *legend = new TLegend(0.2, 0.74, 0.35, 0.8);
        legend->SetTextSize(0.03);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->AddEntry(g, "Reconstruction efficiency", "LEP");
        legend->AddEntry(line, "Overall efficiency ", "L");
        legend->Draw("same");


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

