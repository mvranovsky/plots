#include "PlotEmbeddingJPsi.h"

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}


void PlotEmbeddingJPsi::Make(){



    TH1D* invMass = loadInvMassHist(recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim, getCondition("embedding"));  // automatically subtracts background if bcg tree is loaded
    if(!invMass){
        cout << "Could not load invariant mass histogram. Returning." << endl;
        return;
    }
    
    FitJPsi* fit = new FitJPsi(invMass, "Poly1");
    fit->fitPeak();
    TCanvas *fitCanvas = fit->getCanvas();
    lowLimInvMass = fit->getLowLimitFit();
    topLimInvMass = fit->getHighLimitFit();
    
    if(!fitCanvas){
        cout << "Could not get fit canvas. Returning." << endl;
        return;
    }

    fitCanvas->SetName("fitJPsiCanvas");
    fitCanvas->SetTitle("");
    mEfficiencyFinal = fit->getYield()/starlightTree->GetEntries();
    mEfficiencyErrFinal = fit->getErrYield()/starlightTree->GetEntries();

    DrawEmbeddingpp510JPsi();

    outFile->cd();
    outFile->cd(nameOfEmbeddingJPsiDir);
    fitCanvas->Write();
    fitCanvas->Close();

    cout << "Finished drawing J/psi peak for embedding." << endl;

    cout << "Now creating BEMC efficiency plots" << endl;
    //bemcEfficiency(1, "bemcEff_pairRap");
    bemcEfficiency(2, "bemcEff_eta");
    bemcEfficiency(3, "bemcEff_phi");
    //bemcEfficiency(4, "bemcEff_pTJPsi");
    bemcEfficiency(5, "bemcEff_pT");
    cerr << "Finished drawing BEMC efficiency plots." << endl;
    
    // plots for reconstruction efficiency
    reconstructionEfficiency(1,"recoEff_pairRap");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    reconstructionEfficiency(2,"recoEff_eta");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    reconstructionEfficiency(3,"recoEff_phi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    reconstructionEfficiency(4,"recoEff_pTJPsi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    reconstructionEfficiency(5,"recoEff_pT");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    
    cerr << "Finished drawing reconstruction efficiency plots." << endl;

    runSysStudy(); // run systematic study of the embedding

    cerr << "Finished running systematic study of the embedding." << endl;

    controlPlotsComparison(true); // true = JPsi, false = all
    
    cerr << "Finished comparing control plots " << endl;

    outFile->cd();
    // save all the histograms to canvases into outfile
    handleHistograms(nameOfEmbeddingJPsiDir, "EmbData");
    
    
}

void PlotEmbeddingJPsi::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();
    if(starlightFile) starlightFile->Close();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void PlotEmbeddingJPsi::Init(){
	//define the output file which will store all the canvases

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    mUtil = new Util();


    outFile->cd();
    outFile->mkdir("recoEffPlots");
    outFile->cd();
    outFile->mkdir("bemcEffPlots");
    outFile->cd();
    outFile->mkdir("BemcEfficiency");
    outFile->cd();
    
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, tree, bcgTree);
    
    // load MC tree
    
    TString MCFileName = "slight_EtaCut1_100k.root";
    starlightFile = shared_ptr<TFile>(new TFile(MCFileName, "READ"));
    if(!starlightFile || starlightFile->IsZombie() || !starlightFile->IsOpen()){
        cerr << "Could not get file with starlight data. Returning." << endl;
        return;
    }
    starlightTree = (TTree*)starlightFile->Get("T");
    if(!starlightTree || starlightTree->IsZombie() || starlightTree->GetEntries() <= 0){
        cerr << "Could not get starlight tree. Returning." << endl;
        return;
    }
    
    if(!tree){
        cerr << "Couldn't open tree with data. Returning." << endl;
    	return;
    }
    
    outFile->mkdir(nameOfEmbeddingJPsiDir);
    outFile->cd(nameOfEmbeddingJPsiDir);
    hSysStudyLoose = new TH1D("SystematicStudyLoose", "Systematic Study of embedding of J/psi photoproduction (Loose)",nVariables-3,NHITSFIT,nVariables-4 );
    hSysStudyLoose->SetTitle("Systematic Study of embedding of J/psi photoproduction (Loose)");
    hSysStudyLoose->GetXaxis()->SetTitle("");
    hSysStudyLoose->GetYaxis()->SetTitle("Ratio to Nominal");
    hSysStudyLoose->GetXaxis()->SetBinLabel(1, "N^{Fit}_{hits}");
    hSysStudyLoose->GetXaxis()->SetBinLabel(2, "N^{dEdx}_{hits}");
    hSysStudyLoose->GetXaxis()->SetBinLabel(3, "PID");
    hSysStudyLoose->GetXaxis()->SetBinLabel(4,"#eta");
    hSysStudyLoose->GetYaxis()->SetRangeUser(0.9, 1.1);
    hSysStudyLoose->SetLineColor(kRed);
    hSysStudyLoose->SetMarkerStyle(21);
    hSysStudyLoose->SetMarkerColor(kRed);
    hSysStudyLoose->SetMarkerSize(2);

    hSysStudyTight = (TH1*)hSysStudyLoose->Clone("SystematicStudyTight");
    hSysStudyTight->SetLineColor(kBlue);
    hSysStudyTight->SetMarkerStyle(22);
    hSysStudyTight->SetMarkerColor(kBlue);
    hSysStudyTight->SetMarkerSize(2);

}


void PlotEmbeddingJPsi::controlPlotsComparison(bool justJPsi){

    outFile->cd();
    outFile->mkdir("ComparisonPlots");
    outFile->cd("ComparisonPlots");

    if(!gPad){
        cerr << "gPad is not initialized. Returning." << endl;
        return;
    }

    TString condition = getCondition("embedding");
    if(justJPsi){
        condition += Form(" && invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
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
        hEtaData->GetXaxis()->SetTitle("#eta_{e} [-]");
        hEtaEmb->GetXaxis()->SetTitle("#eta_{e} [-]");
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
        hPhiData->GetXaxis()->SetTitle("#phi_{e} [rad]");
        hPhiEmb->GetXaxis()->SetTitle("#phi_{e} [rad]");
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
        hPtData->GetXaxis()->SetTitle("p^{e}_{T} [GeV/c]");
        hPtEmb->GetXaxis()->SetTitle("p^{e}_{T} [GeV/c]");
        hPtData->GetYaxis()->SetTitle("counts");
        hPtEmb->GetYaxis()->SetTitle("counts");
        plot2Dists(hPtData, hPtEmb, "hPtComparison");
    }else{
        cout << "Could not get pT histograms for comparison." << endl;
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
    outFile->cd();
}



bool PlotEmbeddingJPsi::plot2Dists(TH1 *hData, TH1* hEmb, TString outName){


    if(!hData || !hEmb){
        cout << "Could not get histograms. Returning." << endl;
        return false;
    }
    SetGPad();
    // rescale embedding histogram to match data histogram
    hEmb->Scale(hData->GetEntries()/hEmb->GetEntries());

    CreateCanvas(&canvas, outName, 1200, 800);
    SetGPad();
    SetHistStyle(hData, kBlack, markerStyleTypical);
    SetHistStyle(hEmb, kRed, markerStyleTypical+1);
    
    hEmb->GetYaxis()->SetRangeUser(0, hEmb->GetMaximum()*1.2);
    hEmb->Draw();
    hData->Draw("same");
    DrawEmbeddingpp510JPsi(0.6,0.85,0.9,0.93);

    CreateLegend(&legend, 0.2, 0.8, 0.35, 0.88);
    legend->SetTextSize(0.04);
    legend->AddEntry(hData, "Data", "LEP");
    legend->AddEntry(hEmb, "Embedding (rescaled)", "LEP");
    legend->Draw("same");

    outFile->cd();
    outFile->cd("ComparisonPlots");
    canvas->Write(outName);
    canvas->Close();
    outFile->cd();
    return true;


}

TGraphAsymmErrors* PlotEmbeddingJPsi::bemcEfficiency(int SWITCH, TString nameOfOutput, TString dir ){// 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters

    TH1 *h1, *h2;
    TString xAxisDescription;
    TString condition1, condition2;
    condition1 = Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass );
    condition2 = Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass );
    if(SWITCH ==1){
        tree->Draw(Form("pairRapidity>>hist1(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));

        tree->Draw(Form("pairRapidity>>hist2(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim));
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "y_{J/#psi} [-]";
    }else if(SWITCH == 2){
        tree->Draw(Form("etaHadron0>>hist1(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim) );
        tree->Draw("etaHadron1>>+hist1" );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));

        tree->Draw(Form("etaHadron0>>hist2(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim), "isBemcHit0 == 1" );
        tree->Draw("etaHadron1>>+hist2", "isBemcHit1 == 1" );
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#eta_{e} [-]";
    }else if(SWITCH == 3){
        tree->Draw(Form("phiHadron0>>hist1(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim));
        tree->Draw("phiHadron1>>+hist1");
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));

        tree->Draw(Form("phiHadron0>>hist2(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim), "isBemcHit0 == 1");
        tree->Draw("phiHadron1>>+hist2", "isBemcHit1 == 1");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#phi_{e} [-]";
    }else if(SWITCH == 4){
        tree->Draw(Form("pt>>hist1(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim));
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        tree->Draw(Form("pt>>hist2(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim), "isBemcHit0 == 1 && isBemcHit1 == 1");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
    }else if(SWITCH == 5){
        tree->Draw(Form("pTInGev0>>hist1(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim) );
        tree->Draw("pTInGev1>>+hist1");
        TH1 *hh1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));

        tree->Draw(Form("pTInGev0>>hist2(%d,%f,%f)", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim), "isBemcHit0 == 1");
        tree->Draw("pTInGev1>>+hist2", "isBemcHit1 == 1");
        TH1 *hh2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));

        //const int nEdges = 33;
        //Double_t edges[nEdges] = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.05, 2.15, 2.5};
        const int nEdges = 15;
        Double_t edges[nEdges] = {0.5,0.7,0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,  1.5,  1.6,  1.7,  1.8,  1.9, 2.5};
        
        h1 = (TH1*) hh1->Rebin(nEdges-1, "h1", edges);
        h2 = (TH1*) hh2->Rebin(nEdges-1, "h2", edges);
        xAxisDescription = "p_{T}^{e} [GeV/c]";
    }

    if( !(h1 && h2 && h1->GetEntries() > 0 && h2->GetEntries() > 0) ) {
        cerr << "No histograms for bemc efficiency " << nameOfOutput << endl; 
        return nullptr;
    }

    TGraphAsymmErrors *g = plotEfficiency(h2, h1, nameOfOutput, {0.75, 1.0}, xAxisDescription, "BEMC Efficiency", dir);
    return g;
}


TGraphAsymmErrors* PlotEmbeddingJPsi::reconstructionEfficiency(int SWITCH, TString nameOfOutput, bool sumLastBins){  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters

    TH1 *h1, *h2;
    TString xAxisDescription;
    if(SWITCH == 1){ // pair rapidity
        tree->Draw(Form("pairRapidity>>hist1(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim), Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw(Form("rapVM>>hist2(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim));
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "y_{J/#psi} [-]";
    }else if(SWITCH == 2){ // eta of daughters
        tree->Draw(Form("etaHadron0>>hist1(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim), Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("etaHadron1>>+hist1", Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw(Form("etad1>>hist2(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim));
        starlightTree->Draw("etad2>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#eta_{e} [-]";
    }else if(SWITCH == 3){ // phi of daughters
        tree->Draw(Form("phiHadron0>>hist1(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim), Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("phiHadron1>>+hist1", Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        
        starlightTree->Draw(Form("(phiatemchd1 - 3.14)>>hist2(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim));
        starlightTree->Draw("(phiatemchd2 - 3.14)>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "#phi [rad]";
        
    }else if(SWITCH == 4){ // pT of pair
        tree->Draw(Form("pt>>hist1(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim), Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw(Form("ptVM>>hist2(%d,%f,%f)",recoBins[SWITCH].nBins,recoBins[SWITCH].lowLim,recoBins[SWITCH].topLim));
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        
        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
        
    }else if(SWITCH == 5){ // pT of daughters
        tree->Draw("pTInGev0>>hist1(40,0.5,2.5)", Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("pTInGev1>>+hist1", Form("isBemcHit0 == 1 && isBemcHit1 == 1 && invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        h1 = dynamic_cast<TH1*>(gPad->FindObject("hist1"));
        starlightTree->Draw("ptd1>>hist2(40,0.5,2.5)");
        starlightTree->Draw("ptd2>>+hist2");
        h2 = dynamic_cast<TH1*>(gPad->FindObject("hist2"));
        xAxisDescription = "p_{T}^{e} [GeV/c]";
    }else{
        cout << "Unknown switch. Returning." << endl;
        return nullptr;
    }
    
    if( !(h1 && h2 && h1->GetEntries() > 0 && h2->GetEntries() > 0) ){
        cout << "Could not get histograms for efficiency plot " << nameOfOutput << endl;
        return nullptr;
    }

    TGraphAsymmErrors *g = plotEfficiency(h1, h2, nameOfOutput, {0.0, 0.3}, xAxisDescription, "reconstruction efficiency", "recoEffPlots", sumLastBins);

    return g;

}

TGraphAsymmErrors* PlotEmbeddingJPsi::plotEfficiency(TH1 *h1, TH1 *h2, TString nameOfOutput, vector<double> range, TString xAxisDescription,TString yAxisDescription,  TString dir, bool sumLastBins){

    
    CreateCanvas(&canvas, nameOfOutput, 1200, 800);
    SetGPad();
    h1->Sumw2();
    h2->Sumw2();

    // instead of dividing the graphs, se TGraphAssymErrors and calculate the efficiency of reconstruction

    //cout << "Efficiency for " << nameOfOutput << ": " << h1->GetEntries() << " / " << h2->GetEntries() << endl;
    h1->SetBinContent(0, 0);
    h1->SetBinContent(h1->GetNbinsX()+1, 0);


    TEfficiency *eff = new TEfficiency(*h1, *h2);
    TGraphAsymmErrors* g = eff->CreateGraph();
    g->SetMarkerColor(kBlue);
    g->SetName(nameOfOutput);
    g->SetTitle(nameOfOutput);
    g->GetXaxis()->SetTitle(xAxisDescription);
    g->GetYaxis()->SetTitle(yAxisDescription);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kBlue); 
    g->GetYaxis()->SetRangeUser(range[0], range[1]);
    g->Draw("AEP");

    DrawEmbeddingpp510JPsi(0.7, 0.85, 0.89, 0.9);

    // draw a TLine at overall efficiency
    double overallEfficiency = h1->GetEntries()/h2->GetEntries();
    TLine *line = new TLine(h1->GetXaxis()->GetXmin() , overallEfficiency, h1->GetXaxis()->GetXmax(), overallEfficiency);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    // create a legend
    TLegend *legend = new TLegend(0.2, 0.71, 0.35, 0.8);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(g, "Efficiency", "LEP");
    legend->AddEntry(line, "Overall efficiency ", "L");
    legend->Draw("same");


    outFile->cd();
    if(nameOfOutput != ""){
        outFile->cd(dir);
        canvas->Write(nameOfOutput);

        //save to a file
        TFile *effFile = new TFile("AnaJPsiRecoEff.root","UPDATE");
        if(effFile && !effFile->IsZombie()){
            effFile->cd();
            g->Write(nameOfOutput);
            effFile->Close();
            //cout << "Reconstruction efficiency saved to AnaJPsiRecoEff.root" << endl;
        }else{
            cout << "Could not create or open AnaJPsiRecoEff.root file." << endl;
        }
    }

    if(sumLastBins){
        
        TH1D *t1 = new TH1D("t1", "t1", 3, 0, 3);
        TH1D *t2 = new TH1D("t2", "t2", 3, 0, 3);
        for(int i = 0; i < 3; i++){
            double t1val = (h1->GetBinContent(h1->GetNbinsX() - i*3) + h1->GetBinContent(h1->GetNbinsX() - (i*3 + 1)) + h1->GetBinContent(h1->GetNbinsX() - (i*3 + 2)));
            double t1err = sqrt( pow(h1->GetBinError(h1->GetNbinsX() - i*3), 2) + pow(h1->GetBinError(h1->GetNbinsX() - (i*3 + 1)), 2) + pow(h1->GetBinError(h1->GetNbinsX() - (i*3 + 2)), 2) );
            double t2val = ( h2->GetBinContent(h2->GetNbinsX() - i*3) + h2->GetBinContent(h2->GetNbinsX() - (i*3 + 1)) + h2->GetBinContent(h2->GetNbinsX() - (i*3 + 2)) );
            double t2err = sqrt( pow(h2->GetBinError(h2->GetNbinsX() - i*3), 2) + pow(h2->GetBinError(h2->GetNbinsX() - (i*3 + 1)), 2) + pow(h2->GetBinError(h2->GetNbinsX() - (i*3 + 2)), 2) );
    
            //cout << "Now filling histograms for recoEff: " << t1val << "/" << t2val << endl;
            //cout << "with errors: " << t1err << ", " << t2err << endl;
            t1->SetBinContent(i+1, t1val);
            t1->SetBinError(i+1, t1err);
            t2->SetBinContent(i+1, t2val);
            t2->SetBinError(i+1, t2err);

        }

        
        TEfficiency *T = new TEfficiency(*t1,*t2);
        TGraphAsymmErrors *TG = T->CreateGraph();

        // remove the last 3 points in g
        for(int i = 0 ; i < 9; i++){
            g->RemovePoint(g->GetN() -1);
        }


        g->SetPoint(g->GetN(), 1.125, TG->GetY()[2]);
        g->SetPointError(g->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(2), TG->GetErrorYlow(2));
        //cout << "Just set point: " << TG->GetY()[2] << endl;

        g->SetPoint(g->GetN(), 1.275, TG->GetY()[1]);
        g->SetPointError(g->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(1), TG->GetErrorYlow(1));
        //cout << "Just set point: " << TG->GetY()[1] << endl;

        g->SetPoint(g->GetN(), 1.425, TG->GetY()[0]);
        g->SetPointError(g->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(0), TG->GetErrorYlow(0));
        //cout << "Just set point: " << TG->GetY()[2] << endl;

    }
    
    return g;

}
void PlotEmbeddingJPsi::changeBinning(int SWITCH, int nBins, double low, double top){
    if(SWITCH <= 5 && SWITCH >= 0){
        recoBins[SWITCH].nBins = nBins;
        recoBins[SWITCH].lowLim = low;
        recoBins[SWITCH].topLim = top;
    }else{
        cerr << "ERROR: To change binning, SWITCH argument has to be from 0 up to 5" << endl;
    }

    return;

}


bool PlotEmbeddingJPsi::runStudy(int VAR,TString condition){

    CreateCanvas(&canvas, TString("SysStudyEmb_") + mUtil->nameOfVariable(VAR), 1200, 800);

    vector<int> yields;
    outFile->cd();
    outFile->mkdir(mUtil->nameOfVariable(VAR));
    for(int i = 0; i < 3; i++){

        //cout << "Variation name: " << condition + TString(" ") + mUtil->variationName(i) << endl;
        //cout << "Condition for plotting: " << getCondition(condition + TString(" ") + mUtil->variationName(i)) << endl;
        // load signal
        tree->Draw(Form("invMass>>histMass(%d,%f,%f)", recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim), getCondition(condition + TString(" ") + mUtil->variationName(i)) );
        TH1 *invMassHist = (TH1*)gPad->GetPrimitive("histMass");
        if(!invMassHist || invMassHist->GetEntries() == 0){ 
            cout << "Empty or no hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        TH1D *hSignal = (TH1D*)invMassHist->Clone(Form("hSignal_%d", i));

        // load background
        bcgTree->Draw(Form("invMass>>massBcg(%d,%f,%f)", recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim), getCondition(condition + TString(" ") + mUtil->variationName(i)) );
        TH1D *hBackground = (TH1D*)gPad->GetPrimitive("massBcg");
        if(!hBackground){
            cout << "Empty or no background hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        hSignal->Add(hBackground, -1);

        // fit and obtain yield
        FitJPsi *fit = new FitJPsi(hSignal, "poly1");
        fit->fitPeak();
        yields.push_back(fit->getYield());
        DrawEmbeddingpp510JPsi();
        fit->saveCanvas(outFile, mUtil->nameOfVariable(VAR) + "_" + mUtil->variationName(i), mUtil->nameOfVariable(VAR));

    }

    yieldResults[mUtil->nameOfVariable(VAR)] = yields;

    hSysStudyTight->SetBinContent(VAR+1, 1.0*yields[1]/yields[0] );
    hSysStudyLoose->SetBinContent(VAR+1, 1.0*yields[2]/yields[0] );

    return true;

}



void PlotEmbeddingJPsi::runSysStudy(){

    TString opt = "embedding ";

    for(int i = 0; i < nVariables-3; i++){
        if(!runStudy(i, opt + mUtil->nameOfVariable(i))){
            cout << "Error when running systematic study of " << mUtil->nameOfVariable(i) << endl;
        }
    }

    saveSysStudyYieldsHists();

    cout << "Finished with systematic study for all variables" << endl;
    
}

void PlotEmbeddingJPsi::saveSysStudyYieldsHists(){
    
    
    CreateCanvas(&canvas, "SysStudyEmb", 1200, 800);

    SetGPad(false, 0.12, 0.06, 0.11, 0.06);

    canvas->Clear();
    canvas->SetName("SysStudyOverview");
    canvas->cd();

    
    //SetHistStyle(hSysStudyTight, kBlue, 20);
    //SetHistStyle(hSysStudyLoose, kRed, 21);
    hSysStudyLoose->SetTitle("");
    hSysStudyTight->SetTitle("");
    hSysStudyTight->SetMarkerSize(2);
    hSysStudyTight->SetMarkerSize(2);
    hSysStudyTight->Draw("P");
    hSysStudyLoose->Draw("same P");
    
    DrawSTARInternal();
    // draw the nominal line from the min of the histogram to the max
    
    TLine* line = new TLine(hSysStudyTight->GetXaxis()->GetBinCenter(NHITSFIT) + 0.5, 1, hSysStudyTight->GetXaxis()->GetBinCenter(ETA) + 1.5, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(1.2);
    line->Draw("same");
    
    CreateLegend(&legend, 0.2,0.8, 0.45, 0.92);
    legend->AddEntry(hSysStudyTight, "Tight Condition", "lp");
    legend->AddEntry(hSysStudyLoose, "Loose Condition", "lp");
    legend->AddEntry(line, "Nominal Condition", "lp");
    legend->Draw("same");
    
    
    canvas->Update();
    outFile->cd();
    hSysStudyTight->Write();
    hSysStudyLoose->Write();
    canvas->Write("SysStudyOverview");
    
}
