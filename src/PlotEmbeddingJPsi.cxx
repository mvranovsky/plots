#include "PlotEmbeddingJPsi.h"

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}


void PlotEmbeddingJPsi::Make(){


    TH1D* invMass = loadInvMassHist(recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim, getCondition("embedding"), false);  // automatically subtracts background if bcg tree is loaded
    if(!invMass){
        cout << "Could not load invariant mass histogram. Returning." << endl;
        return;
    }

    cout << "InvMass number of entries: " << invMass->GetEntries() << endl;

    //check if JPsi embedding or not    
    if(isJPsiEmbedding(invMass)){ // only for JPsi embedding


        FitJPsi* fit = new FitJPsi(invMass, "cb Poly1");
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

        cout << "Now plotting reconstruction efficiency..." << endl;
        // plots for reconstruction efficiency
        reconstructionEfficiency(1,"recoEff_pairRap");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(2,"recoEff_eta");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        reconstructionEfficiency(3,"recoEff_phi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        TGraphAsymmErrors *gPtJPsi = reconstructionEfficiency(4,"recoEff_pTJPsi");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        TGraphAsymmErrors *gPtElectrons = reconstructionEfficiency(5,"recoEff_pT");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters

        outFile->cd();
        gPtElectrons->Write();
        gPtJPsi->Write();
        
    }else{ // for gamma gamma -> e+ e- embedding
        int nBins = 25;
        double low = 0.0;
        double top = 5.0;

        invMass = loadInvMassHist(nBins, low, top, getCondition("embedding"), false);  // automatically subtracts background if bcg tree is loaded

        if(!invMass || invMass->GetEntries() == 0){
            cout << "Could not load invariant mass histogram for gamma gamma -> e+ e- embedding. Returning." << endl;
            return;
        }

        FitJPsi* fit = new FitJPsi(invMass, "Continuum");
        fit->setLegendPosition(0.61,0.6,0.8,0.8);

        fit->setFitRangeLow(0.5);
        fit->setFitRangeHigh(4.0);
        fit->fitPeak();

        fit->addContinuumFunction(0.3772, -0.5518, 0.0487);

        TCanvas *c = fit->getCanvas();
        DrawEmbeddingpp510JPsi();
        fit->writeContinuumResult(0.65, 0.6, 0.8, 0.4);
        fit->writeFitResult();

        outFile->cd();
        outFile->cd(nameOfEmbeddingJPsiDir);
        c->SetName("ContinuumFit");
        c->SetTitle("");
        c->Write();

        // create text of continuum

        cout << "Finished drawing continuum for gamma gamma -> e+ e- embedding." << endl;
    }

    
    cout << "Now creating BEMC efficiency plots" << endl;
    // find plots in histFile
    histFile->cd();
    TDirectory* dir = histFile->GetDirectory(nameOfEmbeddingJPsiDir);

    vector<TH1D*> embed, McSim;

    TH1D* h1 = scanDir(dir, "hBemcPtHitEmbed");
    if(h1) {
        h1->SetTitle("Embedding Efficiency");
        embed.push_back(h1);
    }

    TH1D* h2 = scanDir(dir, "hBemcPtAllEmbed");
    if(h2) {
        h2->SetTitle("Embedding Efficiency");
        embed.push_back(h2);
    }

    TH1D* h3 = scanDir(dir, "hBemcPtHitMc");
    if(h3) {
        McSim.push_back(h3);
        h3->SetTitle("MC Simulation Efficiency");
    }
    TH1D* h4 = scanDir(dir, "hBemcPtAllMc");
    if(h4) McSim.push_back(h4);

    if(embed.size() == McSim.size() && embed.size() == 2 ){
        plotEfficiency(embed, McSim, "BemcEfficiency_pT", {0.5, 1.0}, "p_{T}^{e} [GeV/c]", "BEMC efficiency", bemcEffDir, true);
    }else{
        cerr << "Could not get histograms for BEMC efficiency. Skipping" << endl;
        cout << "embed size: " << embed.size() << ", McSim size: " << McSim.size() << endl;
        return;
    }
    

    cerr << "Finished drawing BEMC efficiency plots." << endl;
    

     // --- IGNORE ---
     if(isRunSysStudy()){
         runSysStudy(); // run systematic study of the embedding
     }

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


// Recursive lambda to scan directories
TH1D* PlotEmbeddingJPsi::scanDir(TDirectory *&dir, TString nameOfHist) {
    TIter nextkey(dir->GetListOfKeys());
    TKey *key;


    while ((key = (TKey*)nextkey())) {
        TObject *obj = key->ReadObj();

        // If it's a subdirectory, scan it recursively
        if (obj->InheritsFrom(TDirectory::Class())) {
            TDirectory *subdir = (TDirectory*)obj;
            TH1D* res = scanDir(subdir, nameOfHist);
            if (res) return res;
        }

        // If it's a histogram (TH1 or derived) and name contains substring
        else if (obj->InheritsFrom(TH1::Class())) {
            TString name = obj->GetName();
            if (name.Contains(nameOfHist)) {
                std::cout << "Found histogram: " << name
                            << " in " << dir->GetPath() << std::endl;
                return (TH1D*)obj;
            }
        }
    }

    return nullptr;
};

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
    outFile->mkdir(bemcEffDir);
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



TGraphAsymmErrors* PlotEmbeddingJPsi::reconstructionEfficiency(int SWITCH, TString nameOfOutput, bool plotFitFunc, bool sumLastBins){  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters

    TH1D *h1 = new TH1D("hist1", "hist1", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim);
    TH1D *h2 = new TH1D("hist2", "hist2", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim);
    TString xAxisDescription;

    if(!gPad) cout << "Warning: gPad is null before TTree::Draw()" << endl;
    if(SWITCH == 1){ // pair rapidity
        tree->Draw("pairRapidity>>hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );

        starlightTree->Draw("rapVM>>hist2");
        xAxisDescription = "y_{J/#psi} [-]";
    }else if(SWITCH == 2){ // eta of daughters
        tree->Draw("etaHadron0>>hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("etaHadron1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );


        starlightTree->Draw("etad1>>hist2");
        starlightTree->Draw("etad2>>+hist2");
        xAxisDescription = "#eta_{e} [-]";
    }else if(SWITCH == 3){ // phi of daughters
        tree->Draw("phiHadron0>>hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("phiHadron1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );


        starlightTree->Draw("(phiatemchd1 - 3.14)>>hist2");
        starlightTree->Draw("(phiatemchd2 - 3.14)>>+hist2");
        xAxisDescription = "#phi [rad]";
        
    }else if(SWITCH == 4){ // pT of pair
        tree->Draw("pt>>hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        starlightTree->Draw("ptVM>>hist2");

        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
            
    }else if(SWITCH == 5){ // pT of daughters
        tree->Draw("pTInGev0>>hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        tree->Draw("pTInGev1>>+hist1", Form("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) );
        
        starlightTree->Draw("ptd1>>hist2(40,0.5,2.5)");
        starlightTree->Draw("ptd2>>+hist2");

        xAxisDescription = "p_{T}^{e} [GeV/c]";
    }else{
        cout << "Unknown switch. Returning." << endl;
        return nullptr;
    }
    

    if( !(h1 && h2 && h1->GetEntries() > 0 && h2->GetEntries() > 0) ){
        cout << "Could not get histograms for efficiency plot " << nameOfOutput << endl;
        return nullptr;
    }

    vector<TH1D*> hists;
    h1->SetTitle("Reconstruction Efficiency");
    h2->SetTitle("Reconstruction Efficiency");
    hists.push_back(h1);
    hists.push_back(h2);

    cout << "Reconstruction efficiency histograms for " << nameOfOutput << " have been filled." << endl;

    TGraphAsymmErrors *g = plotEfficiency(hists,{nullptr, nullptr}, nameOfOutput, {0.0, 0.3}, xAxisDescription, "Reconstruction Efficiency", "recoEffPlots", plotFitFunc, sumLastBins);

    /*
    // delete hist1, hist2
    h1->Delete();
    h2->Delete();
    delete h1;
    delete h2;
    
    */

    return g;

}

TGraphAsymmErrors* PlotEmbeddingJPsi::plotEfficiency(vector<TH1D*> h1, vector<TH1D*> h2, TString nameOfOutput, vector<double> range, TString xAxisDescription,TString yAxisDescription,  TString dir, bool plotFitFunc, bool sumLastBins){


    CreateCanvas(&canvas, nameOfOutput, 1200, 800);
    SetGPad();

    
    if(!h1[0] || !h1[1]){
        cout << "Could not get histograms for efficiency plot " << nameOfOutput << endl;
        return nullptr;
    }

    bool plotSecondEff = false;
    if(h2[0] && h2[1]){
        plotSecondEff = true;
        h2[0]->SetBinContent(0, 0);
        h2[0]->SetBinContent(h2[0]->GetNbinsX()+1, 0);
        h2[1]->SetBinContent(0, 0);
        h2[1]->SetBinContent(h2[1]->GetNbinsX()+1, 0);
        //cout << "Efficiency plot: " << h2[0]->GetEntries() << " / " << h2[1]->GetEntries() << endl;

    }
    

    h1[0]->SetBinContent(0, 0);
    h1[0]->SetBinContent(h1[0]->GetNbinsX()+1, 0);
    h1[1]->SetBinContent(0, 0);
    h1[1]->SetBinContent(h1[1]->GetNbinsX()+1, 0);

    
    TEfficiency *eff1 = new TEfficiency(*h1[0], *h1[1]);
    TGraphAsymmErrors* g1 = eff1->CreateGraph();
    g1->SetMarkerColor(kBlue);
    g1->SetName(nameOfOutput + "_embedding");
    g1->SetTitle(nameOfOutput);
    g1->GetXaxis()->SetTitle(xAxisDescription);
    g1->GetYaxis()->SetTitle(yAxisDescription);
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue); 
    g1->SetMinimum(range[0]);
    g1->SetMaximum(range[1]);
    g1->Draw("AEP");

    
    TEfficiency *eff2;
    TGraphAsymmErrors* g2;

    if(plotSecondEff) {
        eff2 = new TEfficiency(*h2[0], *h2[1]);
        g2 = eff2->CreateGraph();
        g2->SetMarkerColor(kRed);
        g2->SetName(nameOfOutput + "_MCSimulator");
        g2->SetTitle(nameOfOutput);
        g2->GetXaxis()->SetTitle(xAxisDescription);
        g2->GetYaxis()->SetTitle(yAxisDescription);
        g2->SetMarkerStyle(20);
        g2->SetMarkerColor(kRed);
        g2->SetLineColor(kRed); 
        //g2->GetYaxis()->SetRangeUser(range[0], range[1]);
        g2->SetMinimum(range[0]);
        g2->SetMaximum(range[1]);
        g2->Draw("same EP");
    }
    cout << "Will plot second efficiency: " << plotSecondEff << endl;

    DrawEmbeddingpp510JPsi(0.7, 0.85, 0.89, 0.9);

    // create a legend
    TLegend *legend = new TLegend(0.2, 0.71, 0.35, 0.8);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(g1, h1[0]->GetTitle(), "LEP");
    if(plotSecondEff)   legend->AddEntry(g2, h2[0]->GetTitle(), "LEP");
    
    if(plotFitFunc){ // plot fit function for BEMC MC simulator

        TF1 *f = new TF1("f","[0] + [1]*(1 + TMath::Erf( (x-[2])/(sqrt(2.)*[3]) ) )",0.5, 2.5);
        // variables for BEMC MC simulator
        double eps0 = 0.00;
        double n = 0.423;
        double pTThr = 0.424;
        double sigma = 0.39;
        // Set initial parameters (important for convergence!)
        f->SetParameters(eps0, n, pTThr, sigma);
        f->SetParNames("eps0","n","pT^{thr}","sigma");
        f->SetLineColor(kGreen);
        f->Draw("same");
        legend->AddEntry(f,"Fit to data","L");
    }
    
    legend->Draw("same");
    
    outFile->cd();
    outFile->cd(dir);
    canvas->Write(nameOfOutput);


    
    if(sumLastBins){
        
    TH1D *t1 = new TH1D("t1", "t1", 3, 0, 3);
        TH1D *t2 = new TH1D("t2", "t2", 3, 0, 3);

        for(int i = 0; i < 3; i++){
            cout << "Now summing last 3 bins: " << i << endl;
            double t1val = (h1[0]->GetBinContent(h1[0]->GetNbinsX() - i*3) + h1[0]->GetBinContent(h1[0]->GetNbinsX() - (i*3 + 1)) + h1[0]->GetBinContent(h1[0]->GetNbinsX() - (i*3 + 2)));
            double t1err = sqrt( pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - i*3), 2) + pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - (i*3 + 1)), 2) + pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - (i*3 + 2)), 2) );
            double t2val = ( h1[1]->GetBinContent(h1[1]->GetNbinsX() - i*3) + h1[1]->GetBinContent(h1[1]->GetNbinsX() - (i*3 + 1)) + h1[1]->GetBinContent(h1[1]->GetNbinsX() - (i*3 + 2)) );
            double t2err = sqrt( pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - i*3), 2) + pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - (i*3 + 1)), 2) + pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - (i*3 + 2)), 2) );
            
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
            g1->RemovePoint(g1->GetN() -1);
        }


        g1->SetPoint(g1->GetN(), 1.125, TG->GetY()[2]);
        g1->SetPointError(g1->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(2), TG->GetErrorYlow(2));
        //cout << "Just set point: " << TG->GetY()[2] << endl;

        g1->SetPoint(g1->GetN(), 1.275, TG->GetY()[1]);
        g1->SetPointError(g1->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(1), TG->GetErrorYlow(1));
        //cout << "Just set point: " << TG->GetY()[1] << endl;

        g1->SetPoint(g1->GetN(), 1.425, TG->GetY()[0]);
        g1->SetPointError(g1->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(0), TG->GetErrorYlow(0));
        //cout << "Just set point: " << TG->GetY()[2] << endl;
        
    }

    
    return g1;
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




bool PlotEmbeddingJPsi::isJPsiEmbedding(TH1D* h){

    int bin = h->FindBin(3.097);
    return (h->GetBinContent(bin) > 500 ? true : false);  //if JPsi bin has more than 500 counts, we assume it's JPsi embedding
}

bool PlotEmbeddingJPsi::isRunSysStudy(){
    TH1D *h1 = loadInvMassHist(100, 2.0, 4.0, getCondition("embedding"), false);
    TH1D* h2 = loadInvMassHist(100, 2.0, 4.0, getCondition("embedding nhitsfit loose"), false);

    if(!h1 || !h2 || h1->GetEntries() == 0 || h2->GetEntries() == 0){
        cerr << "Could not load histograms to decide whether to run systematic study. Returning false." << endl;
        return false;
    }
    cout << "Entries in nominal embedding: " << h1->GetEntries() << ", entries in loose nhitsfit embedding: " << h2->GetEntries() << endl;

    return (h1->GetEntries() == h2->GetEntries() ? false : true);

}