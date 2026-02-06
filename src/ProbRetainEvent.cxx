#include "ProbRetainEvent.h"

ProbRetainEvent::ProbRetainEvent(const std::shared_ptr<TFile>& outFile, TTree* tree)
    : outFile(outFile), tree(tree) {}


ProbRetainEvent::~ProbRetainEvent() {
    outFile = nullptr;
    tree = nullptr;
}


void ProbRetainEvent::Make( TString fitOption ) {  // argument is the name of the output file
    // This function reads a ROOT file containing a TTree with event data
    // and creates a histogram of the probability of retaining an event
    usedFitOption = fitOption;
    if(DEBUGMODE) cout << "Starting ProbRetainEvent::Make()..." << endl;

    defineFitOption(fitOption);

    if(DEBUGMODE) cout << "Reading luminosity file..." << endl;
    readLumiFile(false);

    if(DEBUGMODE) cout << "Getting data from input list..." << endl;
    data = getData();

    if(!(tofMultVeto || bbcEVeto || bbcWVeto)){
        if(DEBUGMODE) cout << "Applying pile-up correction to data..." << endl;
        pileUpCorrection(data);
    }

    if (data.empty()) {
        cerr << "ERROR in ProbRetainEvent::Make(): No data found in the input list." << endl;
        return;
    }
    
    if(DEBUGMODE) cout << "Obtained data from input list. Size: " << data.size() << endl;
    
    if(DEBUGMODE) cout << "Filling graph and performing fit..." << endl;
    TGraphErrors *graph = fillGraph(data);
    if (!graph) {
        cerr << "ERROR in ProbRetainEvent::Make(): Error filling graph." << endl;
        return;
    }

    if(DEBUGMODE) cout << "Obtained graph. Size: " << graph->GetN() << endl;
    

    TCanvas *c = new TCanvas("ProbabilityOfRetainingAnEvent", "ProbRetainEvent", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    graph->SetTitle("");

    // define TGraph
    graph->SetName("ProbRetainEvent");
    graph->SetTitle("ProbRetainEvent");
    graph->GetXaxis()->SetTitle("#Lambda [#mu b^{-1}]");
    graph->GetYaxis()->SetTitle("probability of retaining an event");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    graph->SetLineStyle(1);


    // define an exponential function to fit the data
    TF1 *fit = new TF1("fit", fitFunc, 60, 160);
    
    fit->SetParameter(0, 1.);
    fit->FixParameter(0,1.);
    //fit->SetParLimits(0, 0.95, 1.05);
    fit->SetParameter(1, -0.01);
    
    fit->SetParNames("A", "B");
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);

    TFitResultPtr fitResult = graph->Fit(fit, "RS", "", 60, 160);
    
    if(tofMultVeto){
        graph->GetYaxis()->SetRangeUser(0.5, 1.5);
    }else{
        graph->GetYaxis()->SetRangeUser(0, 1);
    }
    graph->GetXaxis()->SetLimits(60, 160);  // Set x-axis limits
    graph->Draw("AP");
    fit->Draw("same");

    if (fitResult.Get()) {
        int status = fitResult->Status();
        if (status != 0) {
            cout << "Fit failed with status: " << status << endl;
            return;
        }
    }
    if(!(tofMultVeto || bbcEVeto || bbcWVeto)){
        a = fit->GetParameter(0);
        b = fit->GetParameter(1);
        aErr = fit->GetParError(0);
        bErr = fit->GetParError(1);
    }

    DrawSTARInternal();

    double yTop = 0.92, yBottom = 0.77;
    if(tofMultVeto || bbcEVeto || bbcWVeto){
        yTop = 0.5;
        yBottom = 0.32;
    }
    TPaveText *text = new TPaveText(0.77, yTop - 0.1, 0.94, yBottom - 0.1,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetTextAlign(33);
    text -> SetFillStyle(4000);
    text -> SetFillColorAlpha(kWhite, 0.0);  
    text -> SetFillColorAlpha(kWhite, 0.0);  
    text -> SetTextFont(82);
    // print chi square and NDF
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    if(DEBUGMODE) cout << "A = " << a << " +/- " << aErr << endl;
    if(DEBUGMODE) cout << "B = " << b << " +/- " << bErr << endl;
    if(DEBUGMODE) cout << "Chi square/NDF = " << chi2 << "/" << ndf << endl;
    text -> AddText(Form("#chi^{2}/ndf = %.0f/%d", chi2, ndf));
    text -> AddText(fitFuncLegend);
    //text -> AddText(Form("A = %.4f #pm %.4f", fit->GetParameter(0), fit->GetParError(0)));
    text -> AddText(Form("B = %.6f #pm %.6f", fit->GetParameter(1), fit->GetParError(1)));
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
        cerr << "ERROR in ProbRetainEvent::Make(): Error fitting Gaussian." << endl;
        return;
    }

    if(!(tofMultVeto || bbcEVeto || bbcWVeto)) goodRuns = filterProbRetainEvent(data, fitPars[0], fitPars[1], a, b);

    if(DEBUGMODE) cout << "Good runs after probRetainEvent: " << goodRuns.size() << endl;

    // draw lines with TF1 at +- 9 sigma
    c->cd();
    TF1* linePlus = nullptr;
    TF1* lineMinus = nullptr;
    if(!(tofMultVeto || bbcEVeto || bbcWVeto)){

        linePlus = new TF1("linePlusP", fitFunc + " + [2]", 60,160);
        linePlus->SetParameters(fit->GetParameter(0), fit->GetParameter(1), 9 * fitPars[1]);
        linePlus->SetLineColor(kBlack);
        linePlus->SetLineStyle(2);
        linePlus->SetLineWidth(3);
        linePlus->Draw("same");
        
        lineMinus = new TF1("lineMinusP", fitFunc + " - [2]", 60,160);
        lineMinus->SetParameters(fit->GetParameter(0), fit->GetParameter(1), 9 * fitPars[1]);
        lineMinus->SetLineColor(kBlack);
        lineMinus->SetLineStyle(2);
        lineMinus->SetLineWidth(3);
        lineMinus->Draw("same");
    }
    TLegend *legend = new TLegend(0.12, yBottom, 0.3, yTop);
    legend->AddEntry(graph, "Data", "p");
    legend->AddEntry(fit, fitDescription, "lp");
    legend->AddEntry(linePlus, "Fit #pm 9#sigma", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(4000);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(0.04);

    legend->Draw("same");

    c->Modified();
    c->Update();
    outFile->cd();
    outFile->cd(dir);
    if(tofMultVeto){
        c->SetName("ProbRetainEvent_TofMultVeto");
    }else if(bbcEVeto){
        c->SetName("ProbRetainEvent_BbcEVeto");
    }else if(bbcWVeto){
        c->SetName("ProbRetainEvent_BbcWVeto");
    }
    c->Write();

}

void ProbRetainEvent::plotTofMultVeto(TString option){
    if(DEBUGMODE) cout << "Starting ProbRetainEvent::plotTofMultVeto()..." << endl;
    tofMultVeto = true;
    Make(option);
    tofMultVeto = false;
    if(DEBUGMODE) cout << "Finished ProbRetainEvent::plotTofMultVeto()." << endl;
}

void ProbRetainEvent::plotBbcEVeto(TString option){
    if(DEBUGMODE) cout << "Starting ProbRetainEvent::plotBbcEVeto()..." << endl;
    bbcEVeto = true;
    Make(option);
    bbcEVeto = false;
    if(DEBUGMODE) cout << "Finished ProbRetainEvent::plotBbcEVeto()." << endl;
}

void ProbRetainEvent::plotBbcWVeto(TString option){
    if(DEBUGMODE) cout << "Starting ProbRetainEvent::plotBbcWVeto()..." << endl;
    bbcWVeto = true;
    Make(option);
    bbcWVeto = false;
    if(DEBUGMODE) cout << "Finished ProbRetainEvent::plotBbcWVeto()." << endl;
}

void ProbRetainEvent::runCustomList(vector<int> runList, TString option){
    if(DEBUGMODE) cout << "Starting ProbRetainEvent::runCustomList()..." << endl;
    customList = runList;
    runningCustomList = true;
    Make(option);
    runningCustomList = false;
    if(DEBUGMODE) cout << "Finished ProbRetainEvent::runCustomList()." << endl;
}


void ProbRetainEvent::DrawSTARInternal(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(0.04);
   textSTAR -> SetFillColorAlpha(0,0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetFillStyle(4000);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl - 0.05, xr, yr - 0.05,"brNDC");
   textpp510 -> SetTextSize(0.03);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillColor(0);
   textpp510 -> SetFillStyle(4000);
   textpp510 -> SetTextFont(62);
   textpp510->AddText("p+p #sqrt{s} = 510 GeV");
   textpp510 -> Draw("same");
}

vector<double> ProbRetainEvent::fitGaussian(TH1D *h, vector<runInfo> data) {
    // Fit a Gaussian to the histogram
    TCanvas *c1 = new TCanvas("GaussianFit_ProbRetainEvent", "ProbRetainEvent", 1500, 1200);
    c1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    
    if(!h || h->GetEntries() == 0) {
        cerr << "ERROR in ProbRetainEvent::fitGaussian(): histogram is null." << endl;
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

    DrawSTARInternal();

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
    if(tofMultVeto){
        c1->SetName("GaussianFit_ProbRetainEvent_TofMultVeto");
    }else if(bbcEVeto){
        c1->SetName("GaussianFit_ProbRetainEvent_BbcEVeto");
    }else if(bbcWVeto){
        c1->SetName("GaussianFit_ProbRetainEvent_BbcWVeto");
    }
    c1->Write();
    c1->Close();

    vector<double> fitParams;
    fitParams.push_back(mean);
    fitParams.push_back(sigma);
    return fitParams;

}

vector<int> ProbRetainEvent::filterProbRetainEvent(const vector<runInfo>& data, double mean, double sigma, double A, double B) {

    vector<int> filteredData;
    int i = 0;
    int nKickedOutLumi = 0, nKickedOutDiff = 0, nNotKickedOut = 0;
    for (const auto& entry : data) {
        
        double diff = entry.probRetainEvent - exponential(entry.instLumi, A, B);

        if(entry.probRetainEvent < 0.15 && entry.nEventsZBVetoAll > 200){
            cout << "Kicking out run " << entry.runNumber << " for having very low prob retain event: " << entry.nEventsZBVetoPassed << " / " << entry.nEventsZBVetoAll << " = " << entry.probRetainEvent << ". Kicking out " << entry.nEventsJPsi << " jpsi events." << endl;
            continue;
        }


        if(entry.instLumi < 10 ){
            if(DEBUGMODE) cout << "Kicking out run " << entry.runNumber << " for having incorrect inst lumi " << entry.instLumi << ". Kicking out " << entry.nEventsJPsi << " jpsi events." << "Run info: probRetainEvent = " << entry.nEventsZBVetoPassed << "/" << entry.nEventsZBVetoAll << " = " << entry.probRetainEvent << endl;    
            nKickedOutLumi += entry.nEventsJPsi;
            continue;
        } 

        if(diff > mean + 9 * sigma || diff < mean - 9 * sigma){ // tu mozno pridat podmienku, ze toto moze platit iba ak nEventsZBVetoPassed > 20 alebo 50
            if(entry.nEventsZBVetoPassed > 20){
                nKickedOutDiff += entry.nEventsJPsi;
                if(DEBUGMODE) cout << "HIGHLIGHT: kicking out run " << entry.runNumber << " for having pro retain event: " << entry.nEventsZBVetoPassed << " / " << entry.nEventsZBVetoAll << " = " << entry.probRetainEvent << "for luminosity: " << entry.instLumi << ". Because of this, kicking out " << entry.nEventsJPsi << "jpsi events." << endl;
                continue;
            }
        }
        nNotKickedOut += entry.nEventsJPsi;
        filteredData.push_back(entry.runNumber);
        i++;
    }


    if(DEBUGMODE) cout << "Number of JPsi events kicked out due to bad luminosity: " << nKickedOutLumi << " , fraction: " << 1.0*nKickedOutLumi/(nKickedOutLumi + nKickedOutDiff + nNotKickedOut) << endl;
    if(DEBUGMODE) cout << "Number of JPsi events kicked out due to bad probRetainEvent: " << nKickedOutDiff << " , fraction: " << 1.0*nKickedOutDiff/(nKickedOutLumi + nKickedOutDiff + nNotKickedOut) << endl;
    if(DEBUGMODE) cout << "Number of JPsi events NOT kicked out: " << nNotKickedOut << " , fraction: " << 1.0*nNotKickedOut/(nKickedOutLumi + nKickedOutDiff + nNotKickedOut) << endl;
    if(DEBUGMODE) cout << "Total number of JPsi events: " << nKickedOutLumi + nKickedOutDiff + nNotKickedOut << endl;

    return filteredData;
}

void ProbRetainEvent::defineFitOption(TString option){
    
    if(option.Contains("linear")){
        fitFunc = "[0] + [1]*x";
        fitDescription = "Linear Fit";
        fitFuncLegend = "#epsilon_{veto}(#Lambda) = 1 + B#Lambda";

    }else if(option.Contains("exponential")){
        fitFunc = "[0]*exp([1]*x)";
        fitDescription = "Exponential Fit";
        fitFuncLegend = "#epsilon_{veto}(#Lambda) = A exp(B#Lambda)";
    }else{
        cerr << "ERROR in ProbRetainEvent::defineFitOption(): Unknown fit option: " << option << ". Using default exponential fit." << endl;
        fitFunc = "[0]*exp([1]*x)";
        fitDescription = "Exponential Fit";
        fitFuncLegend = "#epsilon_{veto}(#Lambda) = A exp(B#Lambda)";
    }
}


double ProbRetainEvent::exponential(double x, double A, double B) {
    return (A * exp(B * x));
}

double ProbRetainEvent::linear(double x, double A, double B) {
    return (A + B * x);
}   

vector<runInfo> ProbRetainEvent::getData(){

    vector<runInfo> data;

    int mRunNumber = 0;
    double luminosity = 0;
    double instLumi = 0;
    int nEventsZBVetoAll = 0;
    int nEventsZBVetoPassed = 0;
    double nEventsJPsi = 0;
    int nEventsZBTofMultPassed = 0;
    int nEventsZBBbcEPassed = 0;
    int nEventsZBBbcWPassed = 0;
    int mNVerticesZB0, mNVerticesZB1, mNVerticesZB2More;



    tree->SetBranchAddress("runNumber", &mRunNumber);
    tree->SetBranchAddress("luminosity", &luminosity);
    tree->SetBranchAddress("instLumi", &instLumi);
    tree->SetBranchAddress("nEventsZBVetoAll", &nEventsZBVetoAll);
    tree->SetBranchAddress("nEventsZBVetoPassed", &nEventsZBVetoPassed);
    tree->SetBranchAddress("nEventsZBTofMultPassed", &nEventsZBTofMultPassed);
    tree->SetBranchAddress("nEventsZBBbcEPassed", &nEventsZBBbcEPassed);
    tree->SetBranchAddress("nEventsZBBbcWPassed", &nEventsZBBbcWPassed);
    tree->SetBranchAddress("nVerticesZB0", &mNVerticesZB0);
    tree->SetBranchAddress("nVerticesZB1", &mNVerticesZB1);
    tree->SetBranchAddress("nVerticesZB2More", &mNVerticesZB2More);


    tree->SetBranchAddress("nEventsJPsi", &nEventsJPsi);

    if(!tree || tree->GetEntries() == 0){
        cerr << "ERROR in ProbRetainEvent::getData(): tree is null." << endl;
        return data;
    }
    if(DEBUGMODE) cout << "Number of entries in the tree: " << tree->GetEntries() << endl;

    for(int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);

        if(mRunNumber == 0) continue;

        if(runningCustomList){
            if(find(customList.begin(), customList.end(), mRunNumber) == customList.end()){
                continue;
            }
        }

        runInfo entry;
        entry.runNumber = mRunNumber;
        
        if(mInstLumiPerRun.find(mRunNumber) == mInstLumiPerRun.end()){
            if(DEBUGMODE) cout << "WARNING in ProbRetainEvent::getData(): run number " << mRunNumber << " not found in instantaneous luminosity map, inst lumi " << mInstLumiPerRun[mRunNumber] << ", " << luminosity << " nEventsJpsi: " << nEventsJPsi << endl;
            continue;
        }
        if(mLumiPerRun.find(mRunNumber) == mLumiPerRun.end()){
            if(DEBUGMODE) cout << "WARNING in ProbRetainEvent::getData(): run number " << mRunNumber << " not found in luminosity map. Skipping this run." << endl;
            continue;
        }

        if(nEventsJPsi == 0 || mInstLumiPerRun[mRunNumber] == 0){
            if(DEBUGMODE) cout << "WARNING in ProbRetainEvent::getData(): Skipping run number " << mRunNumber << " for having zero nEventsJPsi or zero instLumi." << endl;
            continue;
        }

        entry.luminosity = mLumiPerRun[mRunNumber];
        entry.instLumi = mInstLumiPerRun[mRunNumber];
        entry.probRetainEvent = nEventsZBVetoAll == 0 ? 0 : 1.0*nEventsZBVetoPassed/nEventsZBVetoAll;
        entry.nEventsZBVetoAll = nEventsZBVetoAll;
        entry.nEventsZBVetoPassed = nEventsZBVetoPassed;
        entry.nEventsZBTofMultPassed = nEventsZBTofMultPassed;
        entry.nEventsZBBbcEPassed = nEventsZBBbcEPassed;
        entry.nEventsZBBbcWPassed = nEventsZBBbcWPassed;
        entry.nEventsJPsi = nEventsJPsi;
        entry.nVerticesZB0 = mNVerticesZB0;
        entry.nVerticesZB1 = mNVerticesZB1;
        entry.nVerticesZB2More = mNVerticesZB2More;
        
        data.push_back(entry);
    }

    return data;

}

TGraphErrors* ProbRetainEvent::fillGraph(vector<runInfo> &data){
    
    
    TGraphErrors *graph = new TGraphErrors();
    graph->SetName("lumiPRE");
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Luminosity [pb^{-1}]");
    graph->GetYaxis()->SetTitle("Probability of Retaining Event");

    int i = 0;
    for (const auto& entry : data) {
        //int runNumber = entry.runNumber;
        double luminosity = entry.instLumi;
        double probRetainEvent = entry.probRetainEvent;
        int ZBAll = entry.nEventsZBVetoAll;
        int ZBVeto = entry.nEventsZBVetoPassed;
        //cout << "luminosity: " << luminosity << ". ProbRetainEvent: " << probRetainEvent << endl;
        
        /*
        if(luminosity <= 60 || luminosity > 160){
            i++;
            continue;
        }
        */
        if(probRetainEvent < 0 || probRetainEvent > 1) {
            i++;
            continue;
        }
        if(ZBAll < 20){
            i++;
            continue;
        }

        double err = probRetainEvent * sqrt( 1.0/ZBAll + 1.0/ZBVeto); 
        
        if(tofMultVeto){
            int ZBTofMultPassed = entry.nEventsZBTofMultPassed;
            probRetainEvent =  1.0*ZBTofMultPassed/ZBAll;
            err = probRetainEvent * sqrt( 1.0/ZBAll + 1.0/ZBTofMultPassed);
        }else if(bbcEVeto){
            int ZBBbcEPassed = entry.nEventsZBBbcEPassed;
            probRetainEvent =  1.0*ZBBbcEPassed/ZBAll;
            err = probRetainEvent * sqrt( 1.0/ZBAll + 1.0/ZBBbcEPassed);
        }else if(bbcWVeto){
            int ZBBbcWPassed = entry.nEventsZBBbcWPassed;
            probRetainEvent =  1.0*ZBBbcWPassed/ZBAll;
            err = probRetainEvent * sqrt( 1.0/ZBAll + 1.0/ZBBbcWPassed);
        }
        graph->SetPoint(graph->GetN(), luminosity, probRetainEvent);
        graph->SetPointError(graph->GetN()-1, 0, err);
        i++;
    }
    return graph;
}


double ProbRetainEvent::getLuminosity(const vector<int> &runs, double A, double B, double pileUpA, double pileUpB) {

    double luminosity = 0;
    for(const auto& run : runs) {
        if(usedFitOption.Contains("linear")){
            luminosity += mLumiPerRun[run]*linear(mInstLumiPerRun[run], A, B)*linear(mInstLumiPerRun[run], pileUpA, pileUpB);
        }else{
            luminosity += mLumiPerRun[run]*exponential(mInstLumiPerRun[run], A, B)*exponential(mInstLumiPerRun[run], pileUpA, pileUpB);
        }
    }
    return luminosity;
}



void ProbRetainEvent::readLumiFile( bool isZB ){
   ifstream lumiFile;

   if(isZB){
      lumiFile.open( "/star/u/mvranovsk/star-upcDst/work/lists/luminosityForZB.list" );
   }else{
      lumiFile.open( "/star/u/mvranovsk/star-upcDst/work/lists/lum_perrun_JPsi_HTTP.txt" );
   }

   if (!lumiFile.is_open() ){
      cerr << "\nERROR in AnaGoodRun::readLumiFile(): Problems with opening a file: " << endl;
      return;
   }


   //18149030(runnumber) 12908182(time0) 12909102(time1) 20931(fillnumber) 4.48836088583608e-08(lumi) 1689030.12500000(prescale) 0.706462022811793(livetime) ZDC-trgonly(zdctrg) 0.996090861818113(bs[0]) 1.58758396379823e-08(bs[1]) 5.47782395356853e-09(bs[2]) 3594(nEvents) 438.630043124484(bs[3]) 1271.2384122(bs[4]) 1.83613529620265e-08(bs[5]) 

   int runNumber, time0, time1, fillNumber;
   double lumi, prescale, livetime;
   unsigned int timeDiff;
   string line, zdctrg;
   Double_t bs[6];
   int nEventsLumiFile;

   while ( getline(lumiFile, line) )
   {
        stringstream ss(line);

        ss >> runNumber >> time0 >> time1 >> fillNumber >> lumi >> prescale >> livetime >> zdctrg >> bs[0] >> bs[1] >> bs[2] >> nEventsLumiFile >> bs[3] >> bs[4] >> bs[5];   
        
        timeDiff = time1-time0;
        double instantinousLumi = (livetime*timeDiff > 0 ) ? prescale*lumi*1000000/double(livetime*timeDiff) : 0;
        if( isZB ){
        }else{
            mInstLumiPerRun[runNumber] = instantinousLumi;
            mLumiPerRun[runNumber] = lumi;
            mNEventsLumiFile[runNumber] = nEventsLumiFile;
        }
        //cout << "Read from lumi file, isZB: " << isZB << " - runNumber: " << runNumber << ", inst lumi: " << mInstLumiPerRun[runNumber] << ", just lumi: " << mLumiPerRun[runNumber] << endl;

   }

   lumiFile.close();

   return;
}  



void ProbRetainEvent::pileUpCorrection(vector<runInfo> &data) {
    // apply pile-up correction to the probRetainEvent values in data vector

    TCanvas *cPileUp = new TCanvas("PileUpCorrection", "PileUpCorrection", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t

    TGraphErrors *gPileUp = new TGraphErrors();
    gPileUp->SetName("PileUpCorrection");
    for (auto& entry : data) {
        double instLumi = entry.instLumi;
        double correctionFactor = 1.0 * entry.nVerticesZB1/(entry.nVerticesZB1 + entry.nVerticesZB2More);
        double R = correctionFactor * sqrt(1.0/entry.nVerticesZB1 + 1.0/(entry.nVerticesZB1 + entry.nVerticesZB2More));

        double errLow, errHigh;
        (correctionFactor + R) > 1.0 ? errHigh = 1.0 - correctionFactor : errHigh = R;
        (correctionFactor - R) < 0.0 ? errLow = correctionFactor : errLow = R;
        errHigh > errLow ? R = errLow : R = errHigh;

        gPileUp->SetPoint(gPileUp->GetN(), instLumi, correctionFactor);
        gPileUp->SetPointError(gPileUp->GetN()-1, 0, R);
    }
    gPileUp->SetMarkerColor(kBlue);
    gPileUp->SetMarkerStyle(20);
    gPileUp->SetMarkerSize(1);
    gPileUp->SetLineColor(kBlue);
    gPileUp->SetLineWidth(2);
    gPileUp->GetXaxis()->SetTitle("#Lambda [#mu b^{-1}]");
    gPileUp->GetYaxis()->SetTitle("Pile-up Correction");
    gPileUp->GetXaxis()->SetLimits(50, 170);
    gPileUp->GetYaxis()->SetRangeUser(-0.05, 1.2);
    gPileUp->Draw("AP");

    DrawSTARInternal();

    //double overallCorrection = 1.0 * overallN1V / overallNAll;
    TF1 *fitPileUp = new TF1("fitPileUp", fitFunc, 0, 200);
    fitPileUp->SetParameter(0, 0.77);
    fitPileUp->SetParameter(1, -0.002);
    fitPileUp->SetLineColor(kRed);
    fitPileUp->SetLineWidth(2);
    TFitResultPtr fitResult = gPileUp->Fit(fitPileUp, "RS", "", 60, 150);
    fitPileUp->Draw("same");

    if(DEBUGMODE){
        cout << "Pile-up correction fit parameters: " << endl;
        cout << "A = " << fitPileUp->GetParameter(0) << " +/- " << fitPileUp->GetParError(0) << endl;
        cout << "B = " << fitPileUp->GetParameter(1) << " +/- " << fitPileUp->GetParError(1) << endl;
    }
    
    setPileUpCorrectionA(fitPileUp->GetParameter(0));
    setPileUpCorrectionB(fitPileUp->GetParameter(1));
    setPileUpCorrectionAError(fitPileUp->GetParError(0));
    setPileUpCorrectionBError(fitPileUp->GetParError(1));

    TPaveText *text = new TPaveText(0.75, 0.78, 0.89, 0.67,"brNDC");
    text -> SetTextSize(0.03);
    text -> SetFillColor(0);
    text -> SetTextFont(42);
    text -> SetTextAlign(11);
    text -> SetFillStyle(4000);
    text -> AddText(Form("#chi^{2}/ndf = %.2f/%.d",fitResult->Chi2(), fitResult->Ndf() ));
    text -> AddText(fitFuncLegend);
    text -> AddText(Form("A = %.2f #pm %.2f", fitPileUp->GetParameter(0), 0.02));
    text -> AddText(Form("B = %.4f #pm %.4f", fitPileUp->GetParameter(1), 0.0002));
    text -> Draw("same");

    //create tlegend
    TLegend *legend = new TLegend(0.12 , 0.87, 0.28, 0.93);
    legend->SetTextSize(0.04);
    legend->AddEntry(gPileUp, "Pile-up correction", "p");
    legend->AddEntry(fitPileUp, "#epsilon_{veto}(#Lambda) = A + B#Lambda", "lp");
    legend->SetBorderSize(0);
    legend->SetFillStyle(4000);
    legend->SetFillColorAlpha(0,0);
    legend->Draw("same");

    outFile->cd();
    outFile->cd(dir);
    cPileUp->SetName("PileUpCorrection");
    cPileUp->Write();


}