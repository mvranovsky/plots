#include "PlotBemcEfficiency.h"


PlotBemcEfficiency::PlotBemcEfficiency(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotBemcEfficiency::PlotBemcEfficiency(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotBemcEfficiency::Make() {


    if(DEBUGMODE) cout << "Starting PlotBemcEfficiency::Make()..." << endl;

    int nBins = 40;
    double low = 0.5;
    double top = 2.5;

    bemcEfficiencyPt(nBins, low, top);

    saveEfficiencyGraph();

    deltaDipAnglePlot();

    handleHistograms(nameOfBemcEfficiencyDir, "BEData");

}


void PlotBemcEfficiency::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();

    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE) cout << "The output file is saved: " << outputPosition << endl;
}

void PlotBemcEfficiency::Init(){

    if(DEBUGMODE) cout << "Starting PlotBemcEfficiency::Init()..." << endl;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}

    outFile->cd();
    outFile->mkdir(nameOfBemcEfficiencyDir);
    outFile->cd();
    
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );
    
    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfBemcEfficiencyTree, tree, bcgTree);

    if(DEBUGMODE) cout << "Finished PlotBemcEfficiency::Init()..." << endl;
}

TGraphAsymmErrors* PlotBemcEfficiency::bemcEfficiencyPt(int nBins, double low, double top){

    if(DEBUGMODE) cout << "Starting PlotBemcEfficiency::bemcEfficiencyPt()..." << endl;

    tree->Draw(Form("pTInGev0>>hist(%d,%f,%f)", nBins, low, top) );
    tree->Draw("pTInGev1>>+hist");
    TH1* hh1 = (TH1*)gPad->FindObject("hist");

    tree->Draw(Form("pTInGev0>>hist2(%d,%f,%f)", nBins, low, top) , "isBemcHit0 == 1");
    tree->Draw("pTInGev1>>+hist2", "isBemcHit1 == 1");
    TH1* hh2 = (TH1*)gPad->FindObject("hist2");

    if(!hh1 || !hh2 ){
        cerr << "ERROR in PlotBemcEfficiency::bemcEfficiencyPt(): histograms for pT of BEMC efficiency are either null or empty." << endl;
        return nullptr;
    }

    //const int nEdges = 33;
    //Double_t edges[nEdges] = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.05, 2.15, 2.5};
    const int nEdges = 17;
    Double_t edges[nEdges] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,  1.5,  1.7,  1.9, 2.5};

    TH1 *h1 = (TH1*) hh1->Rebin(nEdges-1, "h1", edges);
    TH1 *h2 = (TH1*) hh2->Rebin(nEdges-1, "h2", edges);

    CreateCanvas(&canvas, "BEMC_Efficiency", 1200, 800);

    TEfficiency *eff = new TEfficiency(*h2, *h1);
    TGraphAsymmErrors *g = eff->CreateGraph();

    SetTGraphStyle(g, kBlue, 20);
    g->GetXaxis()->SetTitle("p_{T}^{e} [GeV/c]");
    g->GetYaxis()->SetTitle("BEMC Efficiency");

    g->Draw("AP");
    
    double overallEff = h2->GetEntries()/h1->GetEntries();
    cout << "Overall BEMC Efficiency from data " << overallEff << endl;


    DrawSTARInternal(0.82, 0.91, 0.94, 0.94);
    
    
    //fit error function to data
    TF1 *f = new TF1("f","[0]*(1 + TMath::Erf( (x-[1])/(sqrt(2.)*[2]) ) )",0.5, 2.5);
    
    // Set initial parameters (important for convergence!)
    f->SetParameters(0.4, 0.8, 0.2);
    f->SetParNames("n","pT^{thr}","sigma");
    f->SetLineColor(kBlue);
    TFitResultPtr r = g->Fit(f, "S R M");
    r->Print("V");
    CreateLegend(&legend, 0.2, 0.81, 0.45, 0.9);
    legend->AddEntry(g, "BEMC Efficiency", "lep");
    legend->AddEntry(f,"Fit Function", "l");
    legend->Draw("same");
    
    int NDF = g->GetN() - f->GetNpar();
    
    TPaveText *text = new TPaveText(0.65,0.17,0.85,0.42, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
    text->SetTextSize(0.03);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetFillColorAlpha(0, 0);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    text->AddText("f(p_{T}) = n#left[1 + erf#left( #frac{p_{T} - p_{T}^{thr}}{#sqrt{2} #sigma} #right) #right]");
    text->AddText("");
    //text->AddText(Form("#epsilon_{0} = %.2f #pm %.2f",f->GetParameter(0),f->GetParError(0)));
    text->AddText(Form("n = %.3f #pm %.3f", f->GetParameter(0),f->GetParError(0)));
    text->AddText(Form("p_{T}^{thr} = %.2f #pm %.2f", f->GetParameter(1),f->GetParError(1)));
    text->AddText(Form("#sigma = %.2f #pm %.2f", f->GetParameter(2),f->GetParError(2)));
    text->AddText(Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",f->GetChisquare(), NDF, f->GetChisquare()/NDF ) );
    
    text->Draw("same hist E");
    
    //setEpsilon(f->GetParameter(0));
    //setEpsilonError(f->GetParError(0));
    setN(f->GetParameter(0));
    setNError(f->GetParError(0));
    setPtThreshold(f->GetParameter(1));
    setPtThresholdError(f->GetParError(1));
    setSigma(f->GetParameter(2));
    setSigmaError(f->GetParError(2));
    setEfficiencyGraph(g);
    
    outFile->cd();
    outFile->cd(nameOfBemcEfficiencyDir);
    canvas->Write();
    h1->Write("hBemcEfficiencyAll");
    h2->Write("hBemcEfficiencyWithHit");
    canvas->Close();
    
    g->GetListOfFunctions()->Clear();

    setEfficiencyGraph(g);
    saveEfficiencyGraph();
    return g;


}

void PlotBemcEfficiency::deltaDipAnglePlot(){
    
    if(DEBUGMODE) cout << "Starting PlotBemcEfficiency::deltaDipAnglePlot()..." << endl;

    tree->Draw("deltaDipAngle>>hist(50,0.,0.05)");
    TH1* h = (TH1*)gPad->FindObject("hist");
    
    if(!h || h->GetEntries() == 0){
        cerr << "ERROR: histogram for deltaDipAngle is null." << endl;
        return;
    }

    CreateCanvas(&canvas, "DeltaDipAngle", 800, 600);
    SetHistStyle(h, kBlue, 20);
    h->GetXaxis()->SetTitle("#delta^{dip} [rad]");
    h->GetYaxis()->SetTitle("counts");
    h->Draw("EP");
    TH1* hClone = (TH1*) h->Clone("hDeltaDipAngle");
    hClone->Draw("same hist");

    DrawSTARInternal();

    canvas->Write();
    canvas->Close();



}

void PlotBemcEfficiency::saveEfficiencyGraph() {

    if(!mEfficiencyGraph){
        cerr << "ERROR in PlotBemcEfficiency::saveEfficiencyGraph(): Efficiency graph is null. Cannot save." << endl;
        return;
    }

    TFile *file = TFile::Open("BEMC_Efficiency_Graph.root", "RECREATE");
    if(!file || file->IsZombie() || !file->IsOpen()){
        cerr << "ERROR in PlotBemcEfficiency::saveEfficiencyGraph(): Could not open file to save BEMC efficiency graph. Returning." << endl;
        return; 
    }

    file->cd();
    mEfficiencyGraph->Write("gBemcEfficiencyPt");
    file->Close();
    outFile->cd();

}