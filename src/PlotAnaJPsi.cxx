#include "PlotAnaJPsi.h"

PlotAnaJPsi::PlotAnaJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotAnaJPsi::PlotAnaJPsi(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotAnaJPsi::Make(){
    

    if(noRomanPots){
        handleHistograms(nameOfAnaJPsiDir, "JPsiDataNoRP");
    }else{
        handleHistograms(nameOfAnaJPsiDir, "JPsiDataWithRP");
    }
    
    // main fit of peak for all triggers combined
    double yieldAll = fitPeakJPsi("", "JPsiPeakFit", 0);
    cout << "Yield of J/Psi for all triggers combined: " << yieldAll << endl;

    
    if(inputPosition.find("sysStudy") != string::npos)  {
        
        if(DEBUGMODE) cout << "About to run systematic study..." << endl;
        
        runSysStudy();
    
        if(DEBUGMODE) cout << "About to run vertexZStudy()" << endl;
        
        vertexZStudy();
    }
    
    if(noRomanPots){
        if(DEBUGMODE) cout << "Creating control plots comparison with embedding..." << endl;
        controlPlotsComparison(true); // true = JPsi, false = all
    }

    if(!noRomanPots){
        if(DEBUGMODE) cout << "Creating pT missing plot..." << endl;
        pTMissingPlot();
        pTMissingCorrelationPlot();
    }

    if(DEBUGMODE) cout << "Loading Starlight file..." << endl;

    loadStarlightTree();

    if(DEBUGMODE) cout << "Creating rapidity dependence plot..." << endl;

    createDependencePlot("rapidity", "rapidityDependence");

    if(DEBUGMODE) cout << "Creating pT dependence plot..." << endl;

    createDependencePlot("pT", "pTDependence");

    if(DEBUGMODE) cout << "Finished creating dependence plots." << endl;

    
    if(DEBUGMODE) cout << "About to run peak fitting study..." << endl;
    
    //peakFittingStudy();   // takes quite long, comment out if running fast
    
    if(DEBUGMODE) cout << "Finished PlotAnaJPsi::Make()" << endl;
}

void PlotAnaJPsi::Finish(){
    
    if(outFile) outFile->Close();
    if(histFile) histFile->Close();
    if(starlightFile) starlightFile->Close();
    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE) cout << "The output file is saved: " << outputPosition << endl;

}

void PlotAnaJPsi::Init(){

	if(!outFile || outFile->IsZombie() ){
		cerr << "ERROR in PlotAnaJPsi::Init():Couldn't open output file. Leaving..." << endl;
        return;
	}

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);


    mUtil = new Util();
    outFile->cd();
    outFile->mkdir(nameOfAnaJPsiDir);

    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "ERROR in PlotAnaJPsi::Init(): Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaJPsiTree, tree, bcgTree);


    if(!tree || !bcgTree){
    	cerr << "ERROR in PlotAnaJPsi::Init(): Couldn't open tree or background tree with data. Returning." << endl;
    	return;
    }

    if(inputPosition.find("noRP") != string::npos || inputPosition.find("NoRP") != string::npos ){
        noRomanPots = true;
    }else{
        noRomanPots = false;
    }
    if(DEBUGMODE) cout << "Running analysis with Roman Pots: " << !noRomanPots << endl;

    // print level for RooFit
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);


    VzStudyMeanGraph = new TGraph();
    VzStudyMeanGraph->SetName("meanGraph");
    VzStudyMeanGraph->SetTitle("Mean V_{Z} vs Fill Number");

    VzStudySigmaGraph = new TGraph();
    VzStudySigmaGraph->SetName("sigmaGraph");
    VzStudySigmaGraph->SetTitle("Sigma V_{Z} vs Fill Number");

    VzStudyEfficiencyFitGraph = new TGraph();
    VzStudyEfficiencyFitGraph->SetName("efficiencyFitGraph");
    VzStudyEfficiencyFitGraph->SetTitle("Efficiency (Fit) vs Fill Number");

    VzStudyEfficiencyHistGraph = new TGraph();
    VzStudyEfficiencyHistGraph->SetName("efficiencyHistGraph");
    VzStudyEfficiencyHistGraph->SetTitle("Efficiency (Hist) vs Fill Number");


}


double PlotAnaJPsi::fitPeakJPsi(TString c, TString outName, int triggerCondition){
    
    // save all the histograms to canvases into outfile
    int nBins = 100;
    double lowerLim = 0.0;
    double upperLim = 5.0;

    TString fitOption = "Continuum cb";
    if(!noRomanPots){
        fitOption = "cb";
        lowerLim = 2.5;
        upperLim = 3.5;
        nBins = 40;
    }

    //cout << "Condition for JPsi analysis: " << getCondition("nsigmaproton nsigmakaon") << endl;
    c == "" ? c = getCondition() : c = c + " && " + getCondition();
    if(DEBUGMODE) cout << "Condition for JPsi fitting: " << c << endl;
    TH1D* invMass = loadInvMassHist(nBins, lowerLim, upperLim, c, true);  // automatically subtracts background if bcg tree is loaded


    if(!invMass){
        cout << "ERROR in PlotAnaJPsi::fitPeakJPsi(): Could not load invariant mass histogram. Returning." << endl;
        return -1;
    }
    invMass->SetName(outName + "_invMass");
    //cout << "InvMass number of entries: " << invMass->GetEntries() << endl;
    
    FitJPsi* fit = new FitJPsi(invMass, fitOption);
    fit->setAlphaLoose(false);
    fit->setNLoose(false);
    fit->setFitRangeLow(0.4);
    fit->setFitRangeHigh(5.0);
    fit->fitPeak();
    fit->writeFitResult();
    fit->saveRooFitResult("fitResultRoofit.root");
    
    if(fitOption.Contains("continuum") || fitOption.Contains("Continuum")) fit->writeContinuumResult();
    setInvMassLimitLow(fit->getLowLimitFit());
    setInvMassLimitHigh(fit->getHighLimitFit());
    
    //fit->printCovarianceMatrix();
    fit->saveRooFitResult("fitResult.root");
    
    TCanvas *fitCanvas = fit->getCanvas();
    if(!fitCanvas){
        cout << "ERROR in PlotAnaJPsi::fitPeakJPsi(): Could not get fit canvas. Returning." << endl;
        return -1;
    }
    

    fitCanvas->SetName(outName);
    fitCanvas->SetTitle("");
    if(triggerCondition == 0 ){  // only in the case of all triggers combined
        mYieldFinal = fit->getYield();
        mYieldErrFinal = fit->getErrYield();
    }

    DrawSTARpp510JPsi();

    if(noRomanPots){
        CreateText("No Roman Pots", 62,0.7,0.75,0.85,0.77);
    }else{
        CreateText("With Roman Pots", 62,0.7,0.75,0.85,0.77);
    }

    if(triggerCondition != 0){
        CreateText(Form("Trigger: HTTP*JPsi (%d)", triggerCondition), 62,0.7,0.65,0.85,0.72);
    }

    outFile->cd(); // Explicitly set our file as current directory
    outFile->cd(nameOfAnaJPsiDir);
    invMass->Write();
    fitCanvas->Write();
    fitCanvas->Close();
    outFile->cd();

    return fit->getYield();
}


void PlotAnaJPsi::controlPlotsComparison(bool justJPsi){

    TString condition = "";
    if(justJPsi){
        condition = Form("invMass > %f && invMass < %f", lowLimInvMass, topLimInvMass);
    }
    if(DEBUGMODE) cout << "Condition for control plots: " << condition << endl;

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
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get etaHadron histogram." << endl;
    }

    if(phiHadron && phiHadron->GetEntries() > 0){
        phiHadron->Write("hPhi");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get phiHadron histogram." << endl;
    }

    if(pTInGev && pTInGev->GetEntries() > 0){
        pTInGev->Write("hPt");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get pTInGev histogram." << endl;
    }

    
    if(vertexZ && vertexZ->GetEntries() > 0){
        vertexZ->Write("hVertexZ");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get vertexZ histogram." << endl;
    }
    
    if(dcaZ && dcaZ->GetEntries() > 0){
        dcaZ->Write("hDCAZ");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get dcaZ histogram." << endl;
    }

    if(dcaXY && dcaXY->GetEntries() > 0){
        dcaXY->Write("hDCAXY");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get dcaXY histogram." << endl;
    }
    

    if(nHitsFit && nHitsFit->GetEntries() > 0){
        nHitsFit->Write("hNHitsFit");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get nHitsFit histogram." << endl;
    }

    if(nHitsDEdx && nHitsDEdx->GetEntries() > 0){
        nHitsDEdx->Write("hNHitsDEdx");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get nHitsDEdx histogram." << endl;
    }

    if(pt && pt->GetEntries() > 0){
        pt->Write("hPtJPsi");
    }else{
        cerr << "ERROR in PlotAnaJPsi::controlPlotsComparison(): Could not get pt histogram." << endl;
    }

    if(DEBUGMODE) cout << "All data plots for comparison have been saved to controlPlots.root file" << endl;
    
    controlFile->Close();
}



void PlotAnaJPsi::vertexZStudy(){
    // function which studies the vertex Z position systematic error


    // loop over all TH1D histograms in histFile
    vector<pair<TH1*, TString>> histograms = GetAllTH1("hVtxZFillNum");
    vector<pair<int,TH1*>> vertexZHistograms;

    int minFillNum = 0;
    int maxFillNum = 0;
    for(auto &hist : histograms){
        
        TString nameStr(hist.first->GetName());
        if( nameStr.Contains("hVtxZFillNum") ){
            // name of hist is in the form hVtxZFillNum_<some_number>
            nameStr.ReplaceAll("hVtxZFillNum_","");
            int fillNum = nameStr.Atoi();
            vertexZHistograms.push_back(make_pair(fillNum, hist.first));
            if(fillNum < minFillNum || minFillNum == 0) minFillNum = fillNum;
            if(fillNum > maxFillNum) maxFillNum = fillNum;
        }
    }

    if(vertexZHistograms.empty()){
        cout << "No vertex Z histograms found. Returning." << endl;
        return;
    }

    cout <<" Found " << vertexZHistograms.size() << " vertex Z histograms." << endl;


    outFile->mkdir("vtxZSysStudy");
    outFile->cd("vtxZSysStudy");

    // create TGraphs where x is the fill number and y is the mean value, then sigma and then the efficiency calculated from the fit and then histogram
    // Define TGraphs for vertex Z study (x: fill number, y: value)
    TGraph *meanGraph = new TGraph();
    meanGraph->SetName("meanGraph");
    meanGraph->SetTitle("Mean V_{Z} vs Fill Number");

    TGraph *sigmaGraph = new TGraph();
    sigmaGraph->SetName("sigmaGraph");
    sigmaGraph->SetTitle("Sigma V_{Z} vs Fill Number");

    TGraph *efficiencyFitGraph = new TGraph();
    efficiencyFitGraph->SetName("efficiencyFitGraph");
    efficiencyFitGraph->SetTitle("Efficiency (Fit) vs Fill Number");

    TGraph *efficiencyHistGraph = new TGraph();
    efficiencyHistGraph->SetName("efficiencyHistGraph");
    efficiencyHistGraph->SetTitle("Efficiency (Hist) vs Fill Number");

    TGraph *efficiencyFitVzGraphLoose = new TGraph();
    efficiencyFitVzGraphLoose->SetName("effFitVzLoose");

    TGraph *efficiencyFitVzGraphTight = new TGraph();
    efficiencyFitVzGraphTight->SetName("effFitVzTight");

    TGraph *efficiencyHistVzGraphLoose = new TGraph();
    efficiencyHistVzGraphLoose->SetName("effHistVzLoose");

    TGraph *efficiencyHistVzGraphTight = new TGraph();
    efficiencyHistVzGraphTight->SetName("effHistVzTight");

    cout << "Starting to fill graphs for vertex Z study..." << endl;
    double averageFit = 0.0, averageFitTight = 0.0, averageFitLoose = 0.0;
    double averageHist = 0.0, averageHistTight = 0.0, averageHistLoose = 0.0;
    for(auto &histPair : vertexZHistograms){
        int fillNum = histPair.first;
        TH1 *hist = histPair.second;

        if(fillNum == 0){
            cout << "Fill number is 0. Skipping." << endl;
            continue;
        }

        if(!hist || hist->GetEntries() == 0){
            cout << "Empty histogram for fill number " << fillNum << ". Skipping." << endl;
            continue;
        }

        CreateCanvas(&canvas, Form("vtxZSysStudy_%d", fillNum), 1200, 800);
        SetHistStyle(hist, kBlue, markerStyleTypical);
        hist->GetXaxis()->SetTitle("V_{Z} [cm]");
        hist->GetYaxis()->SetTitle("Counts");
        hist->Draw("EP");
        TH1* hClone = (TH1*)hist->Clone("hClone");
        hClone->SetLineColor(kBlue);
        hClone->Draw("same hist");

        
        // fit the distribution with a Gaussian
        TF1 *gauss = new TF1("gauss", "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        gauss->SetLineColor(kRed);
        gauss->SetLineWidth(2);
        // fit only in range -120 120
        hist->Fit(gauss, "R", "", -120, 120);
        gauss->Draw("same");
        // integrate the fit function inside 
        double integral = gauss->Integral(-vertexZInCm[0],vertexZInCm[0])/ gauss->Integral(gauss->GetXmin(), gauss->GetXmax());
        double integralTight = gauss->Integral(-vertexZInCm[1],vertexZInCm[1])/ gauss->Integral(gauss->GetXmin(), gauss->GetXmax());
        double integralLoose = gauss->Integral(-vertexZInCm[2],vertexZInCm[2])/ gauss->Integral(gauss->GetXmin(), gauss->GetXmax());

        averageFit += integral;
        averageFitTight += integralTight;
        averageFitLoose += integralLoose;
        
        // integrate the histogram in the same range
        int sum = 0;
        int sumTight = 0;
        int sumLoose = 0;
        for(int i = 1; i <= hist->GetNbinsX(); ++i) {
            double binCenter = hist->GetBinCenter(i);
            if(binCenter >= -vertexZInCm[0] && binCenter <= vertexZInCm[0]) {
                sum += hist->GetBinContent(i);
            }
            if(binCenter <= vertexZInCm[1]) {
                sumTight += hist->GetBinContent(i);
            }
            if(binCenter <= vertexZInCm[2]) {
                sumLoose += hist->GetBinContent(i);
            }
        }
        
        //print the fit results onto the canvas
        DrawSTARInternal(0.77, 0.89, 0.9, 0.93);
        CreateText(Form("Fill number: %d", fillNum), 62, 0.7, 0.8, 0.9, 0.77);
        CreateText("Gaussian fit:", 62, 0.7, 0.74, 0.9, 0.77);
        CreateText(Form("Mean = %.2f +/- %.2f cm", gauss->GetParameter(1), gauss->GetParError(1)),textFont, 0.7, 0.74, 0.9, 0.71);
        CreateText(Form("Sigma = %.2f +/- %.2f cm", gauss->GetParameter(2), gauss->GetParError(2)), textFont, 0.7, 0.71, 0.9, 0.68);
        CreateText(Form("#epsilon_{fit} = %.2f", integral), textFont, 0.72, 0.68, 0.9, 0.65);
        CreateText("Bin integration:", 62, 0.7, 0.65, 0.9, 0.62);
        CreateText(Form("#epsilon_{hist} = %.2f", 1.0*sum/hist->GetEntries()),textFont, 0.7, 0.62, 0.9, 0.59);


        averageHist += 1.0*sum/hist->GetEntries();
        averageHistLoose += 1.0*sumLoose/hist->GetEntries();
        averageHistTight += 1.0*sumTight/hist->GetEntries();

        canvas->Update();
        
        canvas->Write(Form("vtxZSysStudy_%d", fillNum));
        canvas->Close();
        
        // add the fill number, mean, sigma, efficiency to the graphs
        efficiencyFitGraph->SetPoint(efficiencyFitGraph->GetN(), fillNum, integral);
        efficiencyHistGraph->SetPoint(efficiencyHistGraph->GetN(), fillNum, 1.0*sum/hist->GetEntries());

        efficiencyFitVzGraphTight->SetPoint(efficiencyFitVzGraphTight->GetN(), fillNum, integralTight);
        efficiencyFitVzGraphLoose->SetPoint(efficiencyFitVzGraphLoose->GetN(), fillNum, integralLoose);

        efficiencyHistVzGraphTight->SetPoint(efficiencyHistVzGraphTight->GetN(), fillNum, 1.0*sumTight/hist->GetEntries());
        efficiencyHistVzGraphLoose->SetPoint(efficiencyHistVzGraphLoose->GetN(), fillNum, 1.0*sumLoose/hist->GetEntries());

        meanGraph->SetPoint(meanGraph->GetN(), fillNum, gauss->GetParameter(1));
        sigmaGraph->SetPoint(sigmaGraph->GetN(), fillNum, gauss->GetParameter(2));
    }
    if(DEBUGMODE) cout << "Finished filling graphs for vertex Z study..." << endl;

    averageFit /= meanGraph->GetN();
    averageFitLoose /= meanGraph->GetN();
    averageFitTight /= meanGraph->GetN();

    averageHist /= meanGraph->GetN();
    averageHistLoose /= meanGraph->GetN();
    averageHistTight /= meanGraph->GetN();

    VzStudyFitEff = averageFit;
    VzStudyHistEff = averageHist;

    // for systematic error calculation
    VzStudyEffTight = (averageFitTight + averageHistTight) / 2.0;
	VzStudyEffLoose = (averageFitLoose + averageHistLoose) / 2.0;

    if(DEBUGMODE) {
        cout << "Average efficiency (fit): " << averageFit << endl;
        cout << "Average efficiency (hist): " << averageHist << endl;
        cout << "Average efficiency tight (fit): " << averageFitTight << endl;
        cout << "Average efficiency tight (hist): " << averageHistTight << endl;
        cout << "Average efficiency loose (fit): " << averageFitLoose << endl;
        cout << "Average efficiency loose (hist): " << averageHistLoose << endl;
    }

    CreateCanvas(&canvas, "vtxZSysStudyGraphs", 1200, 800);
    SetGPad();
    canvas->SetName("meanVtxZOfFills");
    meanGraph->GetXaxis()->SetTitle("Fill Number");
    meanGraph->GetYaxis()->SetTitle("<V_{Z}> [cm]");
    meanGraph->GetXaxis()->SetLimits(minFillNum, maxFillNum);
    meanGraph->GetYaxis()->SetRangeUser(-10, 10);
    SetTGraphStyle(meanGraph, kBlue, markerStyleTypical);
    meanGraph->Draw("AP");   
    DrawSTARInternal(0.77, 0.89, 0.9, 0.93); 
    canvas->Update();
    canvas->Write("meanGraphVtxZ");
    meanGraph->Write("meanGraphVtxZ_graph");
    //----------------------------------------------------------------------------

    canvas->Clear();
    SetGPad();
    canvas->SetName("sigmaVtxZOfFills");
    sigmaGraph->GetXaxis()->SetTitle("Fill Number");
    sigmaGraph->GetYaxis()->SetTitle("#sigma(V_{Z}) [cm]");
    sigmaGraph->GetXaxis()->SetLimits(minFillNum, maxFillNum);
    sigmaGraph->GetYaxis()->SetRangeUser(40, 80);
    SetTGraphStyle(sigmaGraph, kBlue, markerStyleTypical);
    sigmaGraph->Draw("AP");
    DrawSTARInternal(0.77, 0.89, 0.9, 0.93);
    canvas->Update();
    canvas->Write("sigmaGraphVtxZ");
    sigmaGraph->Write("sigmaGraphVtxZ_graph");
    //----------------------------------------------------------------------------


    canvas->Clear();
    SetGPad();
    canvas->SetName("efficiencyVtxZOfFills");

    efficiencyFitGraph->GetXaxis()->SetTitle("Fill Number");
    efficiencyFitGraph->GetYaxis()->SetTitle("#epsilon");
    efficiencyFitGraph->GetXaxis()->SetLimits(minFillNum, maxFillNum);
    efficiencyFitGraph->GetYaxis()->SetRangeUser(0.8, 1.0);
    SetTGraphStyle(efficiencyFitGraph, kBlue, markerStyleTypical);
    SetTGraphStyle(efficiencyHistGraph, kRed, markerStyleTypical+1);
    efficiencyFitGraph->Draw("APE");
    efficiencyHistGraph->Draw("same PE");

    TF1* fit1 = new TF1("fit", "[0]", minFillNum, maxFillNum);
    fit1->SetParameter(0, averageFit);
    fit1->SetLineColor(kBlue);
    efficiencyFitGraph->Fit(fit1, "R");
    fit1->Draw("same");

    TF1* fit2 = new TF1("fit", "[0]", minFillNum, maxFillNum);
    fit2->SetParameter(0, averageHist);
    fit2->SetLineColor(kRed);
    efficiencyHistGraph->Fit(fit2, "R");
    fit2->Draw("same");


    CreateLegend(&legend, 0.2, 0.2, 0.45, 0.4);
    legend->AddEntry(efficiencyFitGraph, "Fit Integration", "p");
    legend->AddEntry(fit1, Form("Average: %.2f +/- %.2f", round(1000*fit1->GetParameter(0))/10 , 100*fit1->GetParError(0)), "l");
    legend->AddEntry(efficiencyHistGraph, "Histogram Integration", "p");
    legend->AddEntry(fit2, Form("Average: %.2f +/- %.2f", round(1000*fit2->GetParameter(0))/10 , 100*fit2->GetParError(0)), "l");

    legend->Draw("same");

    DrawSTARInternal(0.77, 0.89, 0.9, 0.93);
    canvas->Update();
    canvas->Write("efficiencyVtxZ");
    efficiencyFitGraph->Write("efficiencyFitVtxZ_graph");
    efficiencyHistGraph->Write("efficiencyHistVtxZ_graph");
    canvas->Close();
    outFile->cd();


    return;
}


bool PlotAnaJPsi::runStudy(int VAR,TString condition){

    CreateCanvas(&canvas, TString("SysStudyData_") + mUtil->nameOfVariable(VAR), 1200, 800);

    vector<int> yields;
    outFile->cd();
    outFile->mkdir(mUtil->nameOfVariable(VAR));
    for(int i = 0; i < 3; i++){

        TH1D *hSignal = loadInvMassHist(100, 0.0, 5.0, getCondition(condition, i) , true);
        if(!hSignal){
            cout << "Could not load invariant mass histogram. Returning." << endl;
            continue;
        }

        // fit and obtain yield
        FitJPsi *fit = new FitJPsi(hSignal, "cb continuum");
        fit->fitPeak();
        fit->setFitRangeLow(0.4);
        fit->setFitRangeHigh(5.0);
        DrawSTARpp510JPsi();
        fit->writeFitResult();
        CreateText(mUtil->variationName(i) + TString(" ") + mUtil->variableLatex(VAR), 62,0.7,0.75,0.85,0.77);
        yields.push_back(fit->getYield());
        fit->saveCanvas(outFile, mUtil->nameOfVariable(VAR) + "_" + mUtil->variationName(i), mUtil->nameOfVariable(VAR));
    }

    cout << "Finished systematic study for variable " << mUtil->nameOfVariable(VAR) << ". Yields: ";
    for(auto &y : yields){
        cout << y << " ";
    }
    cout << endl;
    
    yieldResults[mUtil->nameOfVariable(VAR)] = yields;

    return true;

}



void PlotAnaJPsi::runSysStudy(){

    TString opt = "";

    for(int i = 0; i < nVariables; i++){
        if(!runStudy(i,mUtil->nameOfVariable(i))){
            cout << "Error when running systematic study of " << mUtil->nameOfVariable(i) << endl;
        }
    }

    cout << "Finished with systematic study for all variables" << endl;
    
}



void PlotAnaJPsi::peakFittingStudy(){

    TString dir = "PeakFittingStudy";
    outFile->mkdir(dir);

    double minYield, maxYield;
    int nBins;
    if(noRomanPots){
        minYield = 1400;
        maxYield = 1800;
        nBins = 80;
    }else{
        minYield = 50;
        maxYield = 150;
        nBins = 50;
    }
    double minN = 0.5;
    double maxN = 5.0;
    double minAlpha = 0.5;
    double maxAlpha = 5.0;

    TH1D *hYield = new TH1D("hYield", "hMeans; Yield ; counts", nBins, minYield, maxYield);
    TH1D *hN = new TH1D("hN", "hMeans; N ; counts", 90, minN, maxN);
    TH1D *hAlpha = new TH1D("hAlpha", "hMeans; Alpha ; counts",  90, minAlpha, maxAlpha);

    vector<double> lowLim = {0.0,0.2,0.4,0.6,0.8,1.0};
    vector<double> topLim = {4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

    for(auto &low : lowLim){
        for(auto &top : topLim){
            for(int nBins = 80; nBins <= 120; nBins += 5){
                

                TH1D* hSignal = loadInvMassHist(nBins, low, top, getCondition(""), true );  // automatically subtracts background if bcg tree is loaded
            
                FitJPsi *fit = new FitJPsi(hSignal, "cb continuum");
                fit->setAlphaLoose(false);
                fit->setNLoose(false);
                fit->fitPeak();
                int yield = fit->getYield();
                fit->saveCanvas(outFile, Form("PeakFittingStudy_%d_%.1f_%.1f", nBins, low, top), dir);
                if(yield < 0){
                    cout << "Something went wrong with peak fitting. Leaving." << endl;
                    return;
                }

                if(fit->getChiSquare() > 300){
                    cout << "Chi square is too high: " << fit->getChiSquare() << endl;
                    continue;
                }

                hYield->Fill(yield);
                hN->Fill(fit->getN());
                hAlpha->Fill(fit->getAlpha());
            }
        }
    }


    //drawAndFitHist(hN, "FittingStudy_n", minN, maxN );
    //drawAndFitHist(hAlpha, "FittingStudy_alpha", minAlpha, maxAlpha);    
    drawAndFitHist(hYield, "FittingStudy_rawYield", minYield, maxYield);
    
}


void PlotAnaJPsi::drawAndFitHist(TH1D *h, TString outName, double min, double max){


    CreateCanvas(&canvas,outName, 1200, 800);
    SetGPad();

    // Fit the histogram with a Gaussian function
    TF1 *gauss = new TF1("gauss", "gaus", min, max);
    gauss->SetLineColor(kRed);
    gauss->SetLineWidth(2);
    h->Fit(gauss, "R");

    // Draw the results
    canvas->cd();
    h->Draw("E hist");
    gauss->Draw("same");

    DrawSTARInternal(0.77, 0.89, 0.9, 0.93);

    // create text similarly under which display the mean and sigma of the fit
    TPaveText *res = new TPaveText(0.7, 0.63, 0.9, 0.72, "brNDC");
    res->SetTextColor(kBlack);
    res->SetFillStyle(0);
    res->SetFillColorAlpha(0, 0);
    res->SetFillColor(0);
    res->SetBorderSize(0);
    res->SetTextSize(textSize);
    res->SetTextFont(textFont);
    res->SetTextAlign(11);
    res->AddText(Form("Mean: %.0f +/- %.0f", gauss->GetParameter(1), gauss->GetParError(1)));
    res->AddText(Form("Sigma: %.1f +/- %.1f", gauss->GetParameter(2), gauss->GetParError(2)));
    res->Draw("same");

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    canvas->Update();
    outFile->cd("FittingStudy");
    canvas->Write(outName);

    if(outName.Contains("FittingStudy_rawYield") ){
        setPeakFittingMean(gauss->GetParameter(1));
        setPeakFittingSigma(gauss->GetParameter(2));
        setPeakFittingMeanError(gauss->GetParError(1));
        setPeakFittingSigmaError(gauss->GetParError(2));
    }


    return;

}


void PlotAnaJPsi::loadRapidityDependence(int nRapBins){

    double lowRap = 0.0;
    double topRap = 1.0;


    TGraphErrors *rapidityGraph = new TGraphErrors();
    rapidityGraph->SetName("rapidityGraph; y_{J/#psi}; counts");
    for(int bin = 1; bin <= nRapBins; bin++){

        double binLow = ((topRap - lowRap)/nRapBins) * (bin - 1) ;
        double binTop = ((topRap - lowRap)/nRapBins) * bin;

        TString condition = Form("abs(pairRapidity) > %f && abs(pairRapidity) < %f && ", binLow, binTop);

        TH1D* invMass = loadInvMassHist(40, 2.5, 3.5, condition + getCondition(""), true );  // automatically subtracts background if bcg tree is loaded
        
        if(!invMass){
            cerr << "ERROR in PlotAnaJPsi::plotRapidityDependence: Could not load invariant mass histogram for rapidity bin " << bin << ". Returning." << endl;
            return;
        }
        
        FitJPsi* fit = new FitJPsi(invMass, "cb + poly2");
        fit->fitPeak();

        double yield = fit->getYield();
        double yieldErr = fit->getErrYield();

        TCanvas *fitCanvas = fit->getCanvas();
        if(!fitCanvas){
            cerr << "ERROR in PlotAnaJPsi::plotRapidityDependence: Could not get fit canvas for rapidity bin " << bin << ". Returning." << endl;
            return;
        }
        fitCanvas->SetName(Form("fitJPsiCanvas_rapBin_%d", bin));
        fitCanvas->SetTitle("");
        DrawSTARpp510JPsi();
        CreateText(Form("Rapidity J/#psi: [%.1f, %.1f]", binLow, binTop), 62,0.7,0.75,0.85,0.77);

        outFile->cd(); // Explicitly set our file as current directory
        outFile->cd(nameOfAnaJPsiDir);
        invMass->Write(Form("invMass_rapBin_%d", bin));
        fitCanvas->Write();
        fitCanvas->Close();
        outFile->cd();

        mRapidityBins.push_back(make_pair(binLow, binTop));
        mRapidityYields.push_back(yield);
        mRapidityYieldsErrors.push_back(yieldErr);


    }

    return;
}


void PlotAnaJPsi::createDependencePlot(TString option, TString outName){

    if(DEBUGMODE) cout << "Creating dependence plot for option " << option << " ..." << endl;
    

    TH1* histSignal = nullptr;
    TH1* histStarlight = nullptr;

    TString invMassCond = TString::Format("invMass > %f && invMass < %f && ",getInvMassLimitLow(), getInvMassLimitHigh());  // JPsi mass window
    TString starlightCond = "abs(etad1) < 0.9 && abs(etad2) < 0.9";  // starlight condition for both daughters within TPC acceptance
    TString title = "";
    if(option.Contains("rapidity") || option.Contains("y")){
        tree->Draw("pairRapidity>>hist1(20, -1, 1)",invMassCond + getCondition());
        histSignal = (TH1*)gPad->GetPrimitive("hist1");
        bcgTree->Draw("pairRapidity>>hist2(20, -1, 1)",invMassCond + getCondition());
        TH1* histBackground = (TH1*)gPad->GetPrimitive("hist2");

        histSignal->Add(histBackground, -1);  // subtract background

        starlightTree->Draw("rapVM>>hist3(20, -1, 1)", starlightCond);  // no condition for starlight
        histStarlight = (TH1*)gPad->GetPrimitive("hist3");

        title = ";y_{J/#psi} [-]; counts";

    }else if(option.Contains("pT") || option.Contains("pt")){

        tree->Draw("pt>>hist1(15, 0, 1.5)",invMassCond + getCondition());
        histSignal = (TH1*)gPad->GetPrimitive("hist1");
        bcgTree->Draw("pt>>hist2(15, 0, 1.5)",invMassCond + getCondition());
        TH1* histBackground = (TH1*)gPad->GetPrimitive("hist2");

        histSignal->Add(histBackground, -1);  // subtract background

        starlightTree->Draw("ptVM>>hist3(15, 0, 1.5)", starlightCond);  // no condition for starlight
        histStarlight = (TH1*)gPad->GetPrimitive("hist3");

        title = ";p_{T}^{J/#psi} [GeV/c]; counts";

    }else{
        cerr << "ERROR in PlotAnaJPsi::createDependencePlot(): Unknown option " << option << ". Returning." << endl;
        return; 
    }

    
    histStarlight->Scale(1.0 * histSignal->GetEntries() / histStarlight->GetEntries() );  // scale starlight to data entries

    CreateCanvas(&canvas, outName, 1200, 800);
    SetGPad();

    SetHistStyle(histSignal, kBlue, markerStyleTypical);
    SetHistStyle(histStarlight, kRed, markerStyleTypical+1);
    histSignal->SetTitle(title);
    histSignal->SetName(outName);
    histSignal->GetYaxis()->SetRangeUser(0.9, histSignal->GetMaximum() > histStarlight->GetMaximum() ? histSignal->GetMaximum()*1.2 : histStarlight->GetMaximum()*1.2 );
    histSignal->Draw("EP");
    histStarlight->Draw("same EP");
    TH1* hS2 = (TH1*)histSignal->Clone("hS2");
    TH1* hSL2 = (TH1*)histStarlight->Clone("hSL2");
    SetHistStyle(hS2, kBlue, markerStyleTypical);
    SetHistStyle(hSL2, kRed, markerStyleTypical+1);

    hS2->Draw("same hist");
    hSL2->Draw("same hist");

    if(noRomanPots){
        CreateText("No Roman Pots", 62,0.7,0.75,0.85,0.77);
    }else{
        CreateText("With Roman Pots", 62,0.7,0.75,0.85,0.77);
    }


    CreateLegend(&legend,0.7, 0.75, 0.85, 0.67);
    legend->AddEntry(histSignal, "Data", "lep");
    legend->AddEntry(histStarlight, "Starlight MC (scaled)", "lep");
    legend->Draw("same");

    DrawSTARpp510JPsi(0.77, 0.87, 0.9, 0.93);   

    canvas->Update();
    outFile->cd(nameOfAnaJPsiDir);
    canvas->Write(outName);
    canvas->Close();

}


void PlotAnaJPsi::pTMissingPlot(){

    CreateCanvas(&canvas, "pTMissingPlot", 1200, 800);
    SetGPad();

    TString invMassCond = TString::Format("invMass > %f && invMass < %f && ",getInvMassLimitLow(), getInvMassLimitHigh());  // JPsi mass window

    tree->Draw("pTMissing>>hist1(20, 0, 1.5)", invMassCond + getCondition());
    TH1* histSignal = (TH1*)gPad->GetPrimitive("hist1");

    bcgTree->Draw("pTMissing>>hist2(20, 0, 1.5)", invMassCond + getCondition());
    TH1* histBackground = (TH1*)gPad->GetPrimitive("hist2");

    histSignal->Add(histBackground, -1);  // subtract background

    SetHistStyle(histSignal, kBlue, markerStyleTypical);
    TH1* hS2 = (TH1*)histSignal->Clone("hS2");
    SetHistStyle(hS2, kBlue, markerStyleTypical);
    histSignal->SetTitle(";p_{T}^{missing} [GeV/c]; counts");
    histSignal->SetName("pTMissingPlot");
    histSignal->Draw("EP");
    hS2->Draw("same hist");

    DrawSTARpp510JPsi(0.77, 0.87, 0.9, 0.93);

    canvas->Update();
    outFile->cd(nameOfAnaJPsiDir);
    canvas->Write("pTMissingPlot");
    canvas->Close(); 

}


void PlotAnaJPsi::pTMissingCorrelationPlot(){

    CreateCanvas(&canvas, "pTMissingCorrelationPlot", 1200, 800);
    SetGPad();

    TString invMassCond = TString::Format("invMass > %f && invMass < %f && ",getInvMassLimitLow(), getInvMassLimitHigh());  // JPsi mass window

    tree->Draw("pt:pTMissing>>hist1(15, 0, 1.5, 15, 0, 1.5)", invMassCond + getCondition());
    TH2* histSignal = (TH2*)gPad->GetPrimitive("hist1");

    //bcgTree->Draw("pt:pTMissing>>hist2(15, 0, 1.5, 15, 0, 1.5)", invMassCond + getCondition());
    //TH2* histBackground = (TH2*)gPad->GetPrimitive("hist2");

    if(!histSignal || histSignal->GetEntries() == 0 ){
        cerr << "ERROR in PlotAnaJPsi::pTMissingCorrelationPlot(): Could not get histograms from trees. Returning." << endl;
        return;
    }

    //histSignal->Add(histBackground, -1);  // subtract background

    histSignal->GetZaxis()->SetRangeUser(0.9, histSignal->GetMaximum()*1.1);
    SetTH2Style(histSignal);
    histSignal->SetTitle(";p^{J/#psi}_{T} [GeV/c]; p_{T}^{missing} [GeV/c]");
    histSignal->SetName("pTMissingCorrelationPlot");
    histSignal->Draw("COLZ");

    DrawSTARpp510JPsi(0.77, 0.87, 0.9, 0.93);

    canvas->Update();
    outFile->cd(nameOfAnaJPsiDir);
    canvas->Write("pTMissingCorrelationPlot");
    canvas->Close();
}