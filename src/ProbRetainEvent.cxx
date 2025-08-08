#include "ProbRetainEvent.h"



ProbRetainEvent::ProbRetainEvent(const std::unique_ptr<TFile>& outFile, TTree* tree)
    : outFile(outFile.get()), tree(tree) {}


ProbRetainEvent::~ProbRetainEvent() {
    outFile = nullptr;
    tree = nullptr;
}

void ProbRetainEvent::Make( ) {  // argument is the name of the output file
    // This function reads a ROOT file containing a TTree with event data
    // and creates a histogram of the probability of retaining an event
    
    data = getData();
    if (data.empty()) {
        cout << "No data found in the input list." << endl;
        return;
    }
    cout << "Obtained data from input list. Size: " << data.size() << endl;
    
    TGraph *graph = fillGraph(data);
    if (!graph) {
        cout << "Error filling graph." << endl;
        return;
    }

    //saveTopologyEfficiency(data);
    
    cout << "Obtained graph. Size: " << graph->GetN() << endl;
    

    TCanvas *c = new TCanvas("ProbabilityOfRetainingAnEvent", "ProbRetainEvent", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    graph->SetTitle("");

    // define TGraph
    graph->SetName("ProbRetainEvent");
    graph->SetTitle("ProbRetainEvent");
    graph->GetXaxis()->SetTitle("L [#mu b^{-1}]");
    graph->GetYaxis()->SetTitle("probability of retaining an event");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    graph->SetLineStyle(1);

    // define an exponential function to fit the data
    TF1 *fit = new TF1("fit", "[0]*exp([1]*x)", 60, 160);
    fit->SetParameters(1, -0.01);
    fit->SetParNames("A", "B");
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    //fit->SetLineStyle(2);
    graph->Fit(fit, "R");
    graph->GetYaxis()->SetRangeUser(0, 1);

    graph->Draw("AEP");
    fit->Draw("same");

    a = fit->GetParameter(0);
    b = fit->GetParameter(1);
    aErr = fit->GetParError(0);
    bErr = fit->GetParError(1);

    
    TPaveText *textpp510 = new TPaveText(0.77, 0.77, 0.94, 0.96,"brNDC");
    textpp510 -> SetTextSize(0.04);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillStyle(4000);
    textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #rightarrow p #oplus J/#psi #oplus p");
    textpp510->AddText( "#sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");
    
    TPaveText *text = new TPaveText(0.77, 0.6, 0.94, 0.75,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetTextAlign(33);
    text -> SetFillStyle(4000);
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetTextFont(82);
    // print chi square and NDF
    double chi2 = fit->GetChisquare();
    double ndf = fit->GetNDF();
    text -> AddText(Form("#chi^{2}/ndf = %.2f/%.0f", chi2, ndf));
    text -> AddText("f(L) = A*exp(B*L)");
    text -> AddText(Form("A = %.2f #pm %.2f", a, aErr));
    text -> AddText(Form("B = %.4f #pm %.4f", b, bErr));
    text -> Draw("same");
    
    
    TH1D *h = new TH1D("h", "", 100, -0.5,0.5 );
    // loop over all points in the graph 
    for (int i = 0; i < graph->GetN(); i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        double funcVal = exponential(x, fit->GetParameter(0), fit->GetParameter(1)); 
        // difference between y and function value
        h->Fill(y-funcVal);   
    }
    
    vector<double> fitPars= fitGaussian(h,data);
    if(fitPars.size() == 0){
        cout << "Error fitting Gaussian." << endl;
        return;
    }
    
    goodRuns = filterProbRetainEvent(data, fitPars[0], fitPars[1], fit->GetParameter(0), fit->GetParameter(1));
        
    // draw lines with TF1 at +- 9 sigma
    c->cd();
    TF1 *linePlus = new TF1("linePlusP", "[0]*exp([1]*x) + [2]", 60,160);
    linePlus->SetParameters(fit->GetParameter(0), fit->GetParameter(1), 9 * fitPars[1]);
    linePlus->SetLineColor(kBlack);
    linePlus->SetLineStyle(2);
    linePlus->SetLineWidth(3);
    linePlus->Draw("same");
    
    TF1 *lineMinus = new TF1("lineMinusP", "[0]*exp([1]*x) - [2]", 60,160);
    lineMinus->SetParameters(fit->GetParameter(0), fit->GetParameter(1), 9 * fitPars[1]);
    lineMinus->SetLineColor(kBlack);
    lineMinus->SetLineStyle(2);
    lineMinus->SetLineWidth(3);
    lineMinus->Draw("same");
    
    TLegend *legend = new TLegend(0.12, 0.7, 0.3, 0.88);
    legend->AddEntry(graph, "Data", "lp");
    legend->AddEntry(fit, "Exponential fit", "lp");
    legend->AddEntry(linePlus, "Fit #pm 9#sigma", "l");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(0.04);
    legend->Draw("same");

    c->Modified();
    c->Update();


    outFile->cd();
    outFile->cd(dir);
    c->Write("ProbRetainEvent");

}

vector<double> ProbRetainEvent::fitGaussian(TH1D *h, vector<runInfo> data) {
    // Fit a Gaussian to the histogram
    TCanvas *c1 = new TCanvas("GaussianFit_ProbRetainEvent", "ProbRetainEvent", 1500, 1200);
    c1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    
    if(!h || h->GetEntries() == 0) {
        cout << "Error: histogram is null." << endl;
        vector<double> empty;
        return empty;
    }

    outFile->cd();
    dir = "ProbRetainEvent";
    outFile->mkdir(dir);
    outFile->cd(dir);

    
    TF1 *gaus = new TF1("gaus", "gaus", -0.5, 0.5);
    gaus->SetParameters(500, 0, 0.1); // initial parameters for the Gaussian
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    h->SetTitle("");
    h->Fit(gaus, "R");
    
    h->SetLineColor(kBlue);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1);
    h->SetMarkerColor(kBlue);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("data - fit");
    h->GetYaxis()->SetTitle("counts");
    c1->SetLogy();
    h->Draw("EP");

    gaus->Draw("same");
    
    double mean = gaus->GetParameter(1);
    double sigma = gaus->GetParameter(2);

    // create 2 TLine objects plus and minus 9 sigma
    TLine *linePlus = new TLine(mean + 9 * sigma, 0, mean + 9 * sigma, 25);
    linePlus->SetLineColor(kBlack);
    linePlus->SetLineStyle(2);
    linePlus->SetLineWidth(2);
    linePlus->Draw("same");
    TLine *lineMinus = new TLine(mean - 9 * sigma, 0, mean - 9 * sigma, 25);
    lineMinus->SetLineColor(kBlack);
    lineMinus->SetLineStyle(2);
    lineMinus->SetLineWidth(2);
    lineMinus->Draw("same");

    TLegend *legend = new TLegend(0.18 , 0.7, 0.35, 0.88);
    legend->AddEntry(h, "Data - Fit", "lp");
    legend->AddEntry(gaus, "Gaussian fit", "lp");
    legend->AddEntry(linePlus, "Mean #pm 9#sigma", "l");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(0.04);
    legend->Draw("same");

    TPaveText *textpp510 = new TPaveText(0.77, 0.77, 0.94, 0.96,"brNDC");
    textpp510 -> SetTextSize(0.04);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillStyle(4000);
    textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #rightarrow p #oplus J/#psi #oplus p");
    textpp510->AddText( "#sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");

    TPaveText *text = new TPaveText(0.77, 0.6, 0.94, 0.75,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetTextAlign(33);
    text -> SetFillStyle(4000);
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetTextFont(82);
    // print chi square and NDF
    double chi2 = gaus->GetChisquare();
    double ndf = gaus->GetNDF();
    text -> AddText(Form("#chi^{2}/ndf = %.2f/%.0f", chi2, ndf));
    text -> AddText(Form("#mu = %.3f #pm %.3f", gaus->GetParameter(1), gaus->GetParError(1)));
    text -> AddText(Form("#sigma = %.3f #pm %.3f", gaus->GetParameter(2), gaus->GetParError(2)));
    text -> Draw("same");

    outFile->cd();
    outFile->cd(dir);
    c1->Write("ProbRetainEventGaussianFit");
    c1->Close();

    vector<double> fitParams;
    fitParams.push_back(mean);
    fitParams.push_back(sigma);
    return fitParams;

}

vector<int> ProbRetainEvent::filterProbRetainEvent(const vector<runInfo>& data, double mean, double sigma, double A, double B) {

    vector<int> filteredData;
    cout << "Size of all data: " << data.size() << endl;
    int i = 0;
    for (const auto& entry : data) {
        
        double diff = entry.probRetainEvent - exponential(entry.instLumi, A, B);

        if(diff > mean + 9 * sigma || diff < mean - 9 * sigma){
            //cout << "Error: probRetainEvent out of range in file: " << filename << endl;
            continue;
        }
        filteredData.push_back(entry.runNumber);
        i++;
    }

    return filteredData;
}


double ProbRetainEvent::exponential(double x, double A, double B) {
    return A * exp(B * x);
}


vector<runInfo> ProbRetainEvent::getData(){

    vector<runInfo> data;

    int mRunNumber = 0;
    double luminosity = 0;
    double instLumi = 0;
    int nEventsZBVetoAll = 0;
    int nEventsZBVetoPassed = 0;
    int nEventsTEAll = 0;
    int nEventsTEPassed = 0;

    tree->SetBranchAddress("runNumber", &mRunNumber);
    tree->SetBranchAddress("luminosity", &luminosity);
    tree->SetBranchAddress("instLumi", &instLumi);
    tree->SetBranchAddress("nEventsZBVetoAll", &nEventsZBVetoAll);
    tree->SetBranchAddress("nEventsZBVetoPassed", &nEventsZBVetoPassed);
    tree->SetBranchAddress("nEventsTEAll", &nEventsTEAll);
    tree->SetBranchAddress("nEventsTEPassed", &nEventsTEPassed);

    if(!tree || tree->GetEntries() == 0){
        cout << "Error: tree is null." << endl;
        return data;
    }
    cout << "Number of entries in the tree: " << tree->GetEntries() << endl;

    for(int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        //cout << "Processing run number: " << mRunNumber << endl;
        if(mRunNumber == 0) continue;
        runInfo entry;
        entry.runNumber = mRunNumber;
        entry.luminosity = luminosity;
        entry.instLumi = instLumi;
        entry.probRetainEvent = 1.0*nEventsZBVetoPassed/nEventsZBVetoAll;
        entry.nEventsZBVetoAll = nEventsZBVetoAll;
        entry.nEventsZBVetoPassed = nEventsZBVetoPassed;
        entry.TEAll = nEventsTEAll;
        entry.TEPassed = nEventsTEPassed;
        
        data.push_back(entry);
        RunNumToLumi[mRunNumber] = luminosity;
        RunNumToInstLumi[mRunNumber] = instLumi;
    }

    cout << "Obtained data from TTree. Size: " << data.size() << endl;
    return data;

}

TGraph* ProbRetainEvent::fillGraph(vector<runInfo> data){
    
    
    TGraph *graph = new TGraph();
    graph->SetName("lumiPRE");
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
    graph->GetYaxis()->SetTitle("Probability of Retaining Event");


    for (const auto& entry : data) {
        //int runNumber = entry.runNumber;
        double luminosity = entry.instLumi;
        double probRetainEvent = entry.probRetainEvent;
        //cout << "luminosity: " << luminosity << ". ProbRetainEvent: " << probRetainEvent << endl;
        if(luminosity <= 0 || probRetainEvent <= 0){
            continue;
        }
        graph->SetPoint(graph->GetN(), luminosity, probRetainEvent);
    }
    return graph;
}


double ProbRetainEvent::calculateSysError(vector<int> runs) {
    // systematic error will stem from varying the parameters within their errors 

    if(runs.empty()) {
        // calculate systematic error from runs that passed the filter
        runs = goodRuns;
    }

    double lowA = a - aErr;
    double highA = a + aErr;
    double lowB = b - bErr;
    double highB = b + bErr;

    int nVariations = 100; // number of variations to perform

    TCanvas *c = new TCanvas("SysErrorProbRetainEvent", "Systematic Error", 1000, 800);
    TH1D *hLumi = new TH1D("hLumi", "Luminosity",30 , 20,80 );
    hLumi->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
    hLumi->GetYaxis()->SetTitle("counts");

    for(int i = 0; i < nVariations; i++) {
        double A = lowA + (highA - lowA) * (i / double(nVariations));
        for(int j = 0; j < nVariations; j++) {
            double B = lowB + (highB - lowB) * (j / double(nVariations));

            double lumi = getLuminosity(runs, A, B);
            hLumi->Fill(lumi);
        }
    }
    hLumi->SetLineColor(kBlue);
    hLumi->SetMarkerStyle(20);
    hLumi->SetMarkerSize(1);
    hLumi->SetMarkerColor(kBlue);
    hLumi->SetLineWidth(2);
    hLumi->Draw("E HIST");

    // create a gaussian fit to the histogram
    TF1 *gaus = new TF1("gaus", "gaus", hLumi->GetXaxis()->GetXmin(), hLumi->GetXaxis()->GetXmax());
    gaus->SetParameters(50, 50, 4);
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    hLumi->Fit(gaus, "R");
    gaus->Draw("same");

    // print the fit parameters to the canvas
    TPaveText *text = new TPaveText(0.75, 0.7, 0.89, 0.85,"brNDC");
    text -> SetTextSize(0.03);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(33);
    text -> AddText(Form("#mu = %.2f #pm %.2f", gaus->GetParameter(1), gaus->GetParError(1)));
    text -> AddText(Form("#sigma = %.2f #pm %.2f", gaus->GetParameter(2), gaus->GetParError(2)));
    text -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.75, 0.88, 0.89, 0.93,"brNDC");
    textSTAR -> SetTextSize(0.04);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(72);
    textSTAR -> SetTextAlign(33);
    textSTAR->AddText("STAR Internal");
    textSTAR -> Draw("same");

    TPaveText *textpp510;
    textpp510 = new TPaveText(0.75, 0.83, 0.89, 0.88,"brNDC");
    textpp510 -> SetTextSize(0.03);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillColor(0);
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");

    c->Update();
    outFile->cd();
    outFile->cd(dir);
    c->SetName("SysErrorProbRetainEvent");
    c->Write();
    c->Close();
    return gaus->GetParameter(2); // return the sigma of the fit as the systematic error    

}

double ProbRetainEvent::getLuminosity(const vector<int> &runs, double A, double B) {

    double luminosity = 0;
    for(const auto& run : runs) {
        luminosity += RunNumToLumi[run]*exponential(RunNumToLumi[run], A, B);
    }
    return luminosity;
}




void ProbRetainEvent::saveTopologyEfficiency(vector<runInfo> &dat){

    TCanvas *canvas = new TCanvas("TopologyEfficiency", "Topology Efficiency", 1000, 800);

    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.11);
    gPad->SetTickx();
    gPad->SetTicky(); 
    gStyle->SetOptStat(0);  
    gStyle->SetOptTitle(0);

    TGraph *hTrigEff = new TGraph();
    hTrigEff->SetName("hTrigEff");
    hTrigEff->SetTitle("Trigger efficiency");
    hTrigEff->GetXaxis()->SetTitle("Run number");
    hTrigEff->GetYaxis()->SetTitle("Trigger efficiency");

    int nEventsPassed = 0;
    int nEventsTotal = 0;
    int i = 0;
    for (const auto& entry : dat) {
        if(entry.runNumber == 0 ) continue; // skip invalid entries

        cout << "Trigger eff: " << entry.TEPassed << "/" << entry.TEAll << " for run number: " << entry.runNumber << endl;

        hTrigEff->SetPoint(hTrigEff->GetN(), i, 1.0*entry.TEPassed/entry.TEAll);
        nEventsPassed += entry.TEPassed;
        nEventsTotal += entry.TEAll;
        i++;
    }

    double overallEfficiency = 1.0*nEventsPassed / nEventsTotal;
    cout << "Overall topology efficiency: " << overallEfficiency << endl;

    if(!hTrigEff || hTrigEff->GetN() == 0){
        cerr << "Trigger efficiency graph is not initialized. Exiting..." << endl;
        return;
    }
    const int textFont = 42; 
    const Double_t textSize = 0.03;
    const Double_t labelSize = 0.05;


    hTrigEff->SetLineColor(kBlue);
    hTrigEff->SetLineWidth(2);
    hTrigEff->SetMarkerColor(kBlue);
    hTrigEff->SetMarkerStyle(20);
    hTrigEff->SetMarkerSize(1.5);
    hTrigEff->GetXaxis()->SetTitleFont(textFont);
    hTrigEff->GetYaxis()->SetTitleFont(textFont);
    hTrigEff->GetXaxis()->SetLabelFont(textFont);
    hTrigEff->GetYaxis()->SetLabelFont(textFont);
    hTrigEff->GetXaxis()->SetLabelSize(labelSize);
    hTrigEff->GetYaxis()->SetLabelSize(labelSize);
    hTrigEff->GetXaxis()->SetTitleSize(labelSize);
    hTrigEff->GetYaxis()->SetTitleSize(labelSize);
    hTrigEff->GetXaxis()->SetTitle("Run number");
    hTrigEff->GetYaxis()->SetTitle("Trigger efficiency");

    if(!hTrigEff || hTrigEff->GetN() == 0){
        cerr << "Topology efficiency graph is not initialized. Exiting..." << endl;
        return;
    }

    
    hTrigEff->SetMaximum(1.0);
    hTrigEff->SetMinimum(0.0);
    //hTrigEff->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
    hTrigEff->GetYaxis()->SetTitle("Trigger efficiency");
    hTrigEff->SetTitle("");
    hTrigEff->Draw("AEP");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.78, 0.89, 0.91, 0.93,"brNDC");
    textSTAR -> SetTextSize(textSize+0.01);
    textSTAR -> SetFillColorAlpha(0, 0);
    textSTAR -> SetTextFont(72);
    textSTAR -> SetTextAlign(33);
    textSTAR->AddText("STAR Internal");
    textSTAR -> Draw("same");

    TPaveText *textpp510;
    textpp510 = new TPaveText(0.78, 0.84, 0.91, 0.88,"brNDC");
    textpp510 -> SetTextSize(textSize);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillColorAlpha(0,0);
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");




    hTrigEff->GetYaxis()->SetRangeUser(0, 1);
    

    TLine *overallEff = new TLine(0, overallEfficiency, i, overallEfficiency);
    overallEff->SetLineColor(kGreen);
    overallEff->SetLineWidth(2);
    overallEff->Draw("same");

    TLegend *legend = new TLegend(0.15, 0.83, 0.4, 0.91);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextAlign(11);
    legend->SetTextSize(textSize);
    legend->SetTextFont(42);
    legend->SetMargin(0.1);       
    legend->AddEntry(hTrigEff, "Trigger Efficiency", "p");
    //legend->AddEntry(fit, "Poly0 Fit", "l");
    legend->AddEntry(overallEff, "Overall Efficiency", "l");
    legend->Draw("same");

    /*
    TPaveText *text = new TPaveText(0.23, 0.73, 0.45, 0.73,"brNDC");
    text -> SetTextSize(textSize);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(11);   
    //text -> AddText(Form("A = %.2f #pm %.2f", fit->GetParameter(0), fit->GetParError(0)));
    //text -> AddText(Form("B = %.4f #pm %.4f", fit->GetParameter(1), fit->GetParError(1)));
    text -> Draw("same");
    */
    
    outFile->cd();
    outFile->cd(dir);
    canvas->Write();
    canvas->Clear();
}



