#include "PlotEmbeddingJPsi.h"

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotEmbeddingJPsi::PlotEmbeddingJPsi(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}


void PlotEmbeddingJPsi::Make(){


    TH1D* invMass = loadInvMassHist(recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim, getCondition("embedding"), false);  // automatically subtracts background if bcg tree is loaded
    if(!invMass){
        cerr << "ERROR in PlotEmbeddingJPsi::Make(): Could not load invariant mass histogram. Returning." << endl;
        return;
    }

    //check if JPsi embedding or not    
    if(isJPsiEmbedding(invMass)){ // only for JPsi embedding

        if(DEBUGMODE) cout << "Running peak fitting for J/Psi embedding..." << endl;
        
        fitPeakJPsi(invMass);

        if(DEBUGMODE) cout << "Finished drawing J/psi peak for embedding." << endl;


        if(DEBUGMODE) cout << "Now plotting reconstruction efficiency..." << endl;
        // plots for reconstruction efficiency

        reconstructionEfficiency(1,"recoEff_pairRap");  // 1 == pair rapidity, 2 == daughter eta, 3 == daughter phi, 4 == pT of JPsi, 5 == pT of daughters
        
        reconstructionEfficiency(2,"recoEff_eta");  
        
        reconstructionEfficiency(3,"recoEff_phi");  
        
        reconstructionEfficiency(4,"recoEff_pTJPsi");  
        
        reconstructionEfficiency(5,"recoEff_pT");  

        if(DEBUGMODE) cout << "Running systematic study of the embedding..." << endl;
        
        runSysStudy(); // run systematic study of the embedding
        
        if(DEBUGMODE) cout << "Finished running systematic study of the embedding." << endl;
        
        controlPlotsComparison(true); // true = JPsi, false = all
        
        if(DEBUGMODE) cout << "Finished comparing control plots " << endl;
        
    }else{ // for gamma gamma -> e+ e- embedding

        if(DEBUGMODE) cout << "Running continuum fitting for gamma gamma -> e+ e- embedding..." << endl;

        fitContinuum();

        if(DEBUGMODE)  cout << "Finished drawing continuum for gamma gamma -> e+ e- embedding." << endl;

    }        

    if(DEBUGMODE) cout << "Creating BEMC efficiency plots..." << endl;

    createBemcEfficiencyPlots();

    if(DEBUGMODE) cout << "Finished creating BEMC efficiency plots." << endl;

    outFile->cd();

    if(DEBUGMODE) cout << "Handling histograms for embedding..." << endl;

    handleHistograms(nameOfEmbeddingJPsiDir, "EmbData");
    
}

void PlotEmbeddingJPsi::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();
    if(starlightFile) starlightFile->Close();
    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE) cout << "The output file is saved: " << outputPosition << endl;
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
                if(DEBUGMODE) std::cout << "Found histogram: " << name << " in " << dir->GetPath() << std::endl;
                return (TH1D*)obj;
            }
        }
    }

    return nullptr;
};

void PlotEmbeddingJPsi::Init(){
	//define the output file which will store all the canvases

    if(DEBUGMODE) cout << "Starting PlotEmbeddingJPsi::Init()..." << endl;

	if(!outFile || outFile->IsZombie() ){
		cerr << "ERROR in PlotEmbeddingJPsi::Init(): Couldn't open output file. Leaving..." << endl;
        return;
	}

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    mUtil = new Util();


    outFile->cd();
    if(!isDir("recoEffPlots")) outFile->mkdir("recoEffPlots");
    outFile->cd();
    if(!isDir(bemcEffDir)) outFile->mkdir(bemcEffDir);
    outFile->cd();
    if(!isDir("BemcEfficiency")) outFile->mkdir("BemcEfficiency");
    outFile->cd();
    
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "ERROR in PlotEmbeddingJPsi::Init(): Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
    if(DEBUGMODE) cout << "Connecting input tree for embedding..." << endl;
    ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, tree, bcgTree);

    if(!tree){
    cerr << "ERROR in PlotEmbeddingJPsi::Init(): Couldn't open tree with data. Returning." << endl;
    return;
}
    if(DEBUGMODE) cout << "Successfully connected input tree for embedding." << endl; 
    
    if(DEBUGMODE) cout << "Connecting starlight file..." << endl;
    loadStarlightTree();
    if(DEBUGMODE) cout << "Successfully connected starlight tree for embedding." << endl;


    
    outFile->mkdir(nameOfEmbeddingJPsiDir);
    outFile->cd(nameOfEmbeddingJPsiDir);

    if(DEBUGMODE) cout << "Finished with PlotEmbeddingJPsi::Init()" << endl;

}

bool PlotEmbeddingJPsi::isDir(TString dirPath) {
    TDirectory* targetDir = outFile.get(); // root file is a TDirectory
    if (!dirPath.IsNull()) {
        targetDir = dynamic_cast<TDirectory*>(outFile->Get(dirPath));
        if (!targetDir) {
            return false;
        }
    }
    return true;
}


Value PlotEmbeddingJPsi::fitPeakJPsi(TH1D* invMass){

    if(DEBUGMODE) cout << "Now fitting embedding peak in PlotEmbeddingJPsi::fitPeakJPsi()" << endl;

    if(!invMass || invMass == nullptr){
        invMass = loadInvMassHist(recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim, getCondition("embedding"), false);
    }

    if(DEBUGMODE) cout << "Invariant mass distribution ready" << endl;

    FitJPsi* fit = new FitJPsi(invMass, "cb");
    fit->setAlphaLoose(true);
    fit->setNLoose(true);
    fit->setIsEmbedding(true);
    fit->setEmbeddingDescription("Embedding (J/#psi #rightarrow e^{+}e^{-})");

    fit->showDataPeak();
    fit->fitPeak();
    fit->writeFitResult();


    TCanvas *fitCanvas = fit->getCanvas();
    lowLimInvMass = fit->getLowLimitFit();
    topLimInvMass = fit->getHighLimitFit();
    
    if(!fitCanvas){
        cerr << "ERROR in PlotEmbeddingJPsi::fitPeakJPsi(): Could not get fit canvas. Returning." << endl;
        return Value();
    }

    fitCanvas->SetName("fitJPsiCanvas");
    fitCanvas->SetTitle("");
    mEfficiencyFinal = fit->getYield()/starlightTree->GetEntries("abs(etad1) < 0.9 && abs(etad2) < 0.9 && abs(rapVM) < 0.9");
    mEfficiencyErrFinal = fit->getErrYield()/starlightTree->GetEntries("abs(etad1) < 0.9 && abs(etad2) < 0.9 && abs(rapVM) < 0.9");
    mGoodStarlightEntries = starlightTree->GetEntries("abs(etad1) < 0.9 && abs(etad2) < 0.9 && abs(rapVM) < 0.9");

    DrawEmbeddingpp510JPsi();
    
    outFile->cd();
    outFile->cd(nameOfEmbeddingJPsiDir);
    fitCanvas->Write();
    fitCanvas->Close();

    if(DEBUGMODE) cout << "    Finished fitting peak in PlotEmbeddingJPsi::fitPeakJPsi()" << endl;
    
    Value res;
    res.val = fit->getYield();
    res.errTop = fit->getErrYield();
    res.errLow = fit->getErrYield();

    return res;
}

void PlotEmbeddingJPsi::fitContinuum(){
    int nBins = 25;
    double low = 0.0;
    double top = 5.0;

    TH1D* invMass = loadInvMassHist(nBins, low, top, getCondition("embedding nonsigma"), false);  // automatically subtracts background if bcg tree is loaded

    if(!invMass || invMass->GetEntries() == 0){
        cerr << "ERROR in PlotEmbeddingJPsi::fitContinuum(): Could not load invariant mass histogram for gamma gamma -> e+ e- embedding. Returning." << endl;
        return;
    }

    FitJPsi* fit = new FitJPsi(invMass, "Continuum");
    fit->setLegendPosition(0.61,0.6,0.8,0.8);
    fit->setVisualizeError(true);

    fit->setIsEmbedding(true);
    fit->setEmbeddingDescription("Embedding (#gamma#gamma #rightarrow e^{+}e^{-})");
    fit->setFitRangeLow(0.5);
    fit->setFitRangeHigh(4.0);


    fit->showAdditionalContinuumFunction(true);
    fit->fitPeak();


    DrawEmbeddingpp510JPsi();
    fit->writeContinuumResult(0.65, 0.6, 0.8, 0.4);
    fit->writeFitResult();
    TCanvas *c = fit->getCanvas();

    outFile->cd();
    outFile->cd(nameOfEmbeddingJPsiDir);
    c->SetName("ContinuumFit");
    c->SetTitle("");
    c->Write();
}



void PlotEmbeddingJPsi::createBemcEfficiencyPlots(){
    

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
        plotEfficiency(embed, McSim, "BemcEfficiency_pT", {0.5, 1.0}, "p_{T}^{e} [GeV/c]", "BEMC efficiency", bemcEffDir, false);
    }else{
        cerr << "ERROR in PlotEmbeddingJPsi::createBemcEfficiencyPlots(): Could not get histograms for BEMC efficiency. Skipping..." << endl;
        return;
    }

    TEfficiency *eff1 = new TEfficiency(*embed[0], *embed[1]);
    TGraphAsymmErrors* g1 = eff1->CreateGraph();

    TFile *f = TFile::Open("BEMC_Efficiency_Graph.root", "READ");
    if(!f || f->IsZombie() || !f->IsOpen()){
        cerr << "ERROR in PlotEmbeddingJPsi::createBemcEfficiencyPlots(): Could not open BEMC efficiency graph file. Returning." << endl;
        return;
    }   
    TGraphAsymmErrors* g2 = (TGraphAsymmErrors*) f->Get("gBemcEfficiencyPt");
    if(!g2){
        cerr << "ERROR in PlotEmbeddingJPsi::createBemcEfficiencyPlots(): Could not get BEMC efficiency graph from file. Returning." << endl;
        return;
    }
    TGraphAsymmErrors* gClone = (TGraphAsymmErrors*) g2->Clone("gBemcEfficiencyPt_Data");
    plotEfficiency(g1, gClone, "BemcEfficiencyComparison_pT", {0.5, 1.0}, "p_{T}^{e} [GeV/c]", "BEMC efficiency", bemcEffDir, true);
    
    f->Close();
}


void PlotEmbeddingJPsi::controlPlotsComparison(bool justJPsi){

    outFile->cd();
    outFile->mkdir("ComparisonPlots");
    outFile->cd("ComparisonPlots");


    TString condition = getCondition("embedding");
    if(justJPsi){
        condition += Form(" && invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
    }

    TFile *controlFile;

    controlFile = new TFile("controlPlots.root", "read");
    if(!controlFile || controlFile->IsZombie() || !controlFile->IsOpen()){
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not open control plots file. Returning." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get eta histograms for comparison." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get phi histograms for comparison." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get pT histograms for comparison." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get nHitsFit histograms for comparison." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get nHitsDEdx histograms for comparison." << endl;
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
        cerr << "ERROR in PlotEmbeddingJPsi::controlPlotsComparison(): Could not get pTJPsi histograms for comparison." << endl;
    }

    
    controlFile->Close();
    outFile->cd();
}



bool PlotEmbeddingJPsi::plot2Dists(TH1 *hData, TH1* hEmb, TString outName){


    if(!hData || !hEmb){
        cerr << "ERROR in PlotEmbeddingJPsi::plot2Dists(): Could not get histograms. Returning." << endl;
        return false;
    }

    SetGPad();
    // rescale embedding histogram to match data histogram
    hEmb->Scale(hData->GetEntries()/hEmb->GetEntries());

    CreateCanvas(&canvas, outName, 1200, 800);
    SetGPad();
    SetHistStyle(hData, kBlue, markerStyleTypical);
    SetHistStyle(hEmb, kRed, markerStyleTypical+1);
    
    hEmb->GetYaxis()->SetRangeUser(0, hEmb->GetMaximum()*1.2);
    hEmb->Draw("EP");
    hData->Draw("same EP");
    TH1D* hEmbClone = (TH1D*)hEmb->Clone("hEmbClone");
    hEmbClone->SetLineColor(kRed);
    hEmbClone->SetLineWidth(2);
    hEmbClone->Draw("same hist");
    TH1D* hDataClone = (TH1D*)hData->Clone("hDataClone");
    hDataClone->SetLineColor(kBlue);
    hDataClone->SetLineWidth(2);
    hDataClone->Draw("same hist");
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

    TH1D *h1 = new TH1D("hist1" + TString(SWITCH), "hist1", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim);
    TH1D *h2 = new TH1D("hist2" + TString(SWITCH), "hist2", recoBins[SWITCH].nBins, recoBins[SWITCH].lowLim, recoBins[SWITCH].topLim);
    TString xAxisDescription;

    TString condMCFile = "abs(etad1) < 0.9 && abs(etad2) < 0.9 && abs(rapVM) < 0.9";

    if(DEBUGMODE) cout << "Condition:" << getCondition("embedding") << endl;

    if(SWITCH == 1){ // pair rapidity
        tree->Draw("pairRapidity>>hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );

        starlightTree->Draw("rapVM>>hist2", condMCFile);
        xAxisDescription = "y_{J/#psi} [-]";
    }else if(SWITCH == 2){ // eta of daughters
        tree->Draw("etaHadron0>>hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );
        tree->Draw("etaHadron1>>+hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );


        starlightTree->Draw("etad1>>hist2", condMCFile);
        starlightTree->Draw("etad2>>+hist2", condMCFile);
        xAxisDescription = "#eta_{e} [-]";
    }else if(SWITCH == 3){ // phi of daughters
        tree->Draw("phiHadron0>>hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );
        tree->Draw("phiHadron1>>+hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );

        starlightTree->Draw("(phiatemchd1 - 3.14)>>hist2", condMCFile);
        starlightTree->Draw("(phiatemchd2 - 3.14)>>+hist2", condMCFile);
        xAxisDescription = "#phi [rad]";
        
    }else if(SWITCH == 4){ // pT of pair
        tree->Draw("pt>>hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );
        starlightTree->Draw("ptVM>>hist2", condMCFile);

        xAxisDescription = "p_{T}^{J/#psi} [GeV/c]";
            
    }else if(SWITCH == 5){ // pT of daughters
        tree->Draw("pTInGev0>>hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );
        tree->Draw("pTInGev1>>+hist1", TString::Format("invMass > %f && invMass < %f",lowLimInvMass, topLimInvMass ) + " && " + getCondition("embedding") );
        
        starlightTree->Draw("ptd1>>hist2(40,0.5,2.5)", condMCFile);
        starlightTree->Draw("ptd2>>+hist2", condMCFile);

        xAxisDescription = "p_{T}^{e} [GeV/c]";
    }else{
        cerr << "ERROR in PlotEmbeddingJPsi::reconstructionEfficiency(): Unknown switch. Returning." << endl;
        return nullptr;
    }
    

    if( !(h1 && h2 && h1->GetEntries() > 0 && h2->GetEntries() > 0) ){
        cerr << "ERROR in PlotEmbeddingJPsi::reconstructionEfficiency(): Could not get histograms for efficiency plot " << nameOfOutput << endl;
        return nullptr;
    }

    if(DEBUGMODE) cout << "Successfully loaded histograms: contents: " << h1->GetEntries() << " / " << h2->GetEntries() << endl;
    outFile->cd();
    h1->SetName(nameOfOutput + "_reco");
    h2->SetName(nameOfOutput + "_gen");
    h1->Write();
    h2->Write();
    vector<TH1D*> hists;
    h1->SetTitle("Reconstruction Efficiency");
    h2->SetTitle("Reconstruction Efficiency");
    hists.push_back(h1);
    hists.push_back(h2);


    TGraphAsymmErrors *g = plotEfficiency(hists,{nullptr, nullptr}, nameOfOutput, {0.0, 0.3}, xAxisDescription, "Reconstruction Efficiency", "recoEffPlots", plotFitFunc, sumLastBins);

    return g;

}

TGraphAsymmErrors* PlotEmbeddingJPsi::plotEfficiency(vector<TH1D*> h1, vector<TH1D*> h2, TString nameOfOutput, vector<double> range, TString xAxisDescription,TString yAxisDescription,  TString dir, bool plotFitFunc, bool sumLastBins){


    CreateCanvas(&canvas, nameOfOutput, 1200, 800);
    SetGPad();

    
    if(!h1[0] || !h1[1]){
        cerr << "ERROR in PlotEmbeddingJPsi::plotEfficiency(): Could not get histograms for efficiency plot " << nameOfOutput << endl;
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
            double t1val = (h1[0]->GetBinContent(h1[0]->GetNbinsX() - i*3) + h1[0]->GetBinContent(h1[0]->GetNbinsX() - (i*3 + 1)) + h1[0]->GetBinContent(h1[0]->GetNbinsX() - (i*3 + 2)));
            double t1err = sqrt( pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - i*3), 2) + pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - (i*3 + 1)), 2) + pow(h1[0]->GetBinError(h1[0]->GetNbinsX() - (i*3 + 2)), 2) );
            double t2val = ( h1[1]->GetBinContent(h1[1]->GetNbinsX() - i*3) + h1[1]->GetBinContent(h1[1]->GetNbinsX() - (i*3 + 1)) + h1[1]->GetBinContent(h1[1]->GetNbinsX() - (i*3 + 2)) );
            double t2err = sqrt( pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - i*3), 2) + pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - (i*3 + 1)), 2) + pow(h1[1]->GetBinError(h1[1]->GetNbinsX() - (i*3 + 2)), 2) );
            

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

        g1->SetPoint(g1->GetN(), 1.275, TG->GetY()[1]);
        g1->SetPointError(g1->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(1), TG->GetErrorYlow(1));

        g1->SetPoint(g1->GetN(), 1.425, TG->GetY()[0]);
        g1->SetPointError(g1->GetN()-1, 0.075, 0.075 ,TG->GetErrorYhigh(0), TG->GetErrorYlow(0));
        
    }
    return g1;
}

void PlotEmbeddingJPsi::changeBinning(int SWITCH, int nBins, double low, double top){
    if(SWITCH <= 5 && SWITCH >= 0){
        recoBins[SWITCH].nBins = nBins;
        recoBins[SWITCH].lowLim = low;
        recoBins[SWITCH].topLim = top;
    }else{
        cerr << "ERROR in PlotEmbeddingJPsi::changeBinning(): To change binning, SWITCH argument has to be from 0 up to 5" << endl;
    }

    return;

}


bool PlotEmbeddingJPsi::runStudy(int VAR,TString condition){

    CreateCanvas(&canvas, TString("SysStudyEmb_") + mUtil->nameOfVariable(VAR), 1200, 800);

    vector<int> yields;
    outFile->cd();
    outFile->mkdir(mUtil->nameOfVariable(VAR));
    for(int i = 0; i < 3; i++){

        tree->Draw(Form("invMass>>histMass(%d,%f,%f)", recoBins[0].nBins, recoBins[0].lowLim, recoBins[0].topLim), getCondition(condition + TString(" ") + mUtil->variationName(i)) );
        TH1 *invMassHist = (TH1*)gPad->GetPrimitive("histMass");
        if(!invMassHist || invMassHist->GetEntries() == 0){ 
            cerr << "ERROR in PlotEmbeddingJPsi::runStudy(): Empty or no hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        TH1D *hSignal = (TH1D*)invMassHist->Clone(Form("hSignal_%d", i));

        // fit and obtain yield
        FitJPsi *fit = new FitJPsi(hSignal, "cb");
        fit->setAlphaLoose(true);
        fit->setNLoose(true);
        fit->setIsEmbedding(true);

        fit->fitPeak();
        yields.push_back(fit->getYield());
        fit->writeFitResult();
        fit->setEmbeddingDescription("Embedding (J/#psi #rightarrow e^{+}e^{-})");
        DrawEmbeddingpp510JPsi();
        fit->saveCanvas(outFile, mUtil->nameOfVariable(VAR) + "_" + mUtil->variationName(i), mUtil->nameOfVariable(VAR));

    }

    yieldResults[mUtil->nameOfVariable(VAR)] = yields;
    

    return true;

}

void PlotEmbeddingJPsi::runSysStudy(){

    TString opt = "embedding ";

    for(int i = 0; i < nVariables; i++){
        if(DEBUGMODE) cout << "Running systematic study for variable: " << mUtil->nameOfVariable(i) << endl;
        if(!runStudy(i, opt + mUtil->nameOfVariable(i))){
            cerr << "ERROR in PlotEmbeddingJPsi::runSysStudy(): when running systematic study of " << mUtil->nameOfVariable(i) << endl;
        }
    }
    
}


bool PlotEmbeddingJPsi::isJPsiEmbedding(TH1D* h){

    int bin = h->FindBin(3.097);
    return (h->GetBinContent(bin) > 500 ? true : false);  //if JPsi bin has more than 500 counts, we assume it's JPsi embedding
}

bool PlotEmbeddingJPsi::isRunSysStudy(){
    TH1D *h1 = loadInvMassHist(100, 2.0, 4.0, getCondition("embedding"), false);
    TH1D* h2 = loadInvMassHist(100, 2.0, 4.0, getCondition("embedding nhitsfit loose"), false);

    if(!h1 || !h2 || h1->GetEntries() == 0 || h2->GetEntries() == 0){
        cerr << "ERROR in PlotEmbeddingJPsi::isRunSysStudy(): Could not load histograms to decide whether to run systematic study. Returning false." << endl;
        return false;
    }
    if(DEBUGMODE) cout << "Entries in nominal embedding: " << h1->GetEntries() << ", entries in loose nhitsfit embedding: " << h2->GetEntries() << endl;

    return (h1->GetEntries() == h2->GetEntries() ? false : true);

}

TGraphAsymmErrors* PlotEmbeddingJPsi::plotEfficiency(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TString nameOfOutput, vector<double> range, TString xAxisDescription,TString yAxisDescription,  TString dir, bool plotFitFunc, bool sumLastBins){


    CreateCanvas(&canvas, nameOfOutput, 1200, 800);
    SetGPad();
    canvas->Clear();

    SetTGraphStyle(g1, kRed, 20);
    g1->SetMarkerColor(kRed);
    g1->SetMinimum(range[0]);
    g1->SetMaximum(range[1]);
    
    SetTGraphStyle(g2, kBlue, 21);
    g2->SetMinimum(range[0]);
    g2->SetMaximum(range[1]);
    

    g1->Draw("AEP");
    g2->Draw("same EP");

    DrawSTARInternal(0.7, 0.85, 0.89, 0.9);

    // create a legend
    TLegend *legend = new TLegend(0.2, 0.71, 0.35, 0.8);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(g1, "Efficiency (embedding)", "LEP");
    legend->AddEntry(g2, "Efficiency (data)", "LEP");
    
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
        f->SetLineColor(kBlue);
        f->Draw("same");
        legend->AddEntry(f,"Fit to data","L");
    }
    
    legend->Draw("same");
    
    outFile->cd();
    outFile->cd(dir);
    canvas->Write(nameOfOutput);
    
    return g1;
}
