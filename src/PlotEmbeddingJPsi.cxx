#include "PlotEmbeddingJPsi.h"

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, unique_ptr<TFile> &file): Plot(mInputList, file) {}


void PlotEmbeddingJPsi::Make(){


	int nBins = 40;
    double lowerLim = 2.0;
    double upperLim = 4.0;

    TH1D* invMass = loadInvMassHist(nBins, lowerLim, upperLim, getCondition("embedding"));  // automatically subtracts background if bcg tree is loaded
    if(!invMass){
        cout << "Could not load invariant mass histogram. Returning." << endl;
        return;
    }
    
    FitJPsi* fit = new FitJPsi(invMass, "");
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

    // save all the histograms to canvases into outfile
    handleHistograms(nameOfEmbeddingJPsiDir);
    

    outFile->Close();
    histFile->Close();
    starlightFile->Close();
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
    outFile->mkdir(nameOfEmbeddingJPsiDir);
    outFile->cd();
    outFile->mkdir("recoEffPlots");

    histFile = unique_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, tree, bcgTree);

    // load MC tree

    TString MCFileName = "slight_EtaCut1_100k.root";
    starlightFile = unique_ptr<TFile>(new TFile(MCFileName, "READ"));
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



void PlotEmbeddingJPsi::reconstructionEfficiency(int SWITCH, TString nameOfOutput){  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
    // open tree from starlight


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
    //-----------------------------------------------------------------------------------
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

        DrawEmbeddingpp510JPsi();

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
        outFile->cd("recoEffPlots");
        canvas->Write(nameOfOutput);

        TFile *effFile = new TFile("AnaJPsiRecoEff.root","RECREATE");
        if(effFile && !effFile->IsZombie()){
            effFile->cd();
            g->Write(nameOfOutput);
            effFile->Close();
            cout << "Reconstruction efficiency saved to AnaJPsiRecoEff.root" << endl;
        }else{
            cout << "Could not create or open AnaJPsiRecoEff.root file." << endl;
        }


    }else{
        cout << "Could not get histograms." << endl;
        if(!h1 || h1->GetEntries() == 0){
            cout << "h1 is null or has no entries." << endl;
        }
        if(!h2 || h2->GetEntries() == 0){
            cout << "h2 is null or has no entries." << endl;
        }
    }
    // close the starlight file

    return;
}

bool PlotEmbeddingJPsi::runStudy(int VAR,TString condition){

    CreateCanvas(&canvas, "SystematicStudyOfCuts", 1200, 800);
    
    vector<int> yields;
    outFile->cd();
    outFile->mkdir(mUtil->nameOfVariable(VAR));
    for(int i = 0; i < 3; i++){

        cout << "Variation name: " << condition + TString(" ") + mUtil->variationName(i) << endl;
        cout << "Condition for plotting: " << getCondition(condition + TString(" ") + mUtil->variationName(i)) << endl;
        // load signal
        tree->Draw(Form("invMass>>histMass(%d,%f,%f)", nBins, lowerLim, upperLim), getCondition(condition + TString(" ") + mUtil->variationName(i)) );
        TH1 *invMassHist = (TH1*)gPad->GetPrimitive("histMass");
        if(!invMassHist || invMassHist->GetEntries() == 0){ 
            cout << "Empty or no hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        TH1D *hSignal = (TH1D*)invMassHist->Clone(Form("hSignal_%d", i));

        // load background
        bcgTree->Draw(Form("invMass>>massBcg(%d,%f,%f)", nBins, lowerLim, upperLim), getCondition(condition + TString(" ") + mUtil->variationName(i)) );
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
    
    
    CreateCanvas(&canvas, "SystematicStudyOfCuts", 1200, 800);

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
