#include "SysStudy.h"

SysStudy::SysStudy(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

SysStudy::SysStudy(const string mInputList, unique_ptr<TFile> &file): Plot(mInputList, file) {}

void SysStudy::Make(){

    cout << "Starting systematic study..." << endl;
    
    peakFittingStudy( );  

    runSysStudy();

    vertexZStudy();

    TTree *GRTree, *GRTreeBcg;
    double sysErrorLumi = -1;;
    if(ConnectInputTree(inputPosition, nameOfAnaGoodRunTree, GRTree, GRTreeBcg)) {
        if(GRTree){
            ProbRetainEvent* probRetainEvent = new ProbRetainEvent(outFile, GRTree);
            probRetainEvent->Make();
            sysErrorLumi = probRetainEvent->calculateSysError();
        }
    }
    cout << "Systematic error of luminosity is " << sysErrorLumi << endl;   



    outFile->Close();
    histFile->Close();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void SysStudy::Init(){
	//define the output file which will store all the canvases

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    histFile = unique_ptr<TFile>( new TFile("histFile.root", "read") );

    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }

	//load the tree chain from the input file
    ConnectInputTree(inputPosition, nameOfSysStudyTree, tree, bcgTree);

    
    if(!tree || !bcgTree){
        ConnectInputTree(inputPosition, nameOfEmbeddingJPsiTree, tree,bcgTree);
        
        if(!tree){
            cerr << "Couldn't open tree or background tree with data. Returning." << endl;
            return;
        }else{
            isEmbedding = true;
            cout << "Will be running systematic study of embedding" << endl;
        }
    }

    if(inputPosition.find("noRP") != string::npos || inputPosition.find("NoRP") != string::npos){
        noRomanPots = true;
    }else{
        noRomanPots = false;
    }

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


    hSysStudyLoose = new TH1D("SystematicStudyLoose", "Systematic Study of J/psi photoproduction (Loose)",nVariables-1,NHITSFIT,nVariables-1 );
    hSysStudyLoose->SetTitle("Systematic Study of J/psi photoproduction (Loose)");
    hSysStudyLoose->GetXaxis()->SetTitle("");
    hSysStudyLoose->GetYaxis()->SetTitle("Ratio to Nominal");
    hSysStudyLoose->GetXaxis()->SetBinLabel(NHITSFIT, nameOfHist[NHITSFIT-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(NHITSDEDX, nameOfHist[NHITSDEDX-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(DCAZINCM, nameOfHist[DCAZINCM-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(DCAXYINCM, nameOfHist[DCAXYINCM-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(PID, nameOfHist[PID-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(VERTEXZ, nameOfHist[VERTEXZ-1]);
    hSysStudyLoose->GetXaxis()->SetBinLabel(ETA, nameOfHist[ETA-1]);
    hSysStudyLoose->GetYaxis()->SetRangeUser(0.9, 1.1);
    hSysStudyLoose->SetLineColor(kRed);
    hSysStudyLoose->SetMarkerStyle(21);
    hSysStudyLoose->SetMarkerColor(kRed);
    hSysStudyLoose->SetMarkerSize(2);

    hSysStudyTight = (TH1D*)hSysStudyLoose->Clone("SystematicStudyTight");
    hSysStudyTight->SetLineColor(kBlue);
    hSysStudyTight->SetMarkerStyle(22);
    hSysStudyTight->SetMarkerColor(kBlue);
    hSysStudyTight->SetMarkerSize(2);




    cout << "Successfully initialized SysStudy..." << endl;

}


double SysStudy::GetAverageYValue(TH1 *hist){
    if(!hist || hist->GetEntries() == 0){
        cout << "Empty histogram. Cannot calculate average." << endl;
        return -999;
    }

    double sum = 0.0;
    int count = 0;

    for(int i = 1; i <= hist->GetNbinsX(); ++i) {
        double binContent = hist->GetBinContent(i);
        if(binContent > 0) {
            sum += binContent;
            count++;
        }
    }

    if(count == 0) {
        cout << "No valid bins found. Cannot calculate average." << endl;
        return -999;
    }

    return sum / count;
}


int SysStudy::studyDiscriminator(TString canvasName){

    if(canvasName.Contains("nHitsFit")){
        return 0;
    }else if(canvasName.Contains("nHitsDedx")){
        return 1;
    }else if(canvasName.Contains("dcaZInCm")){
        return 2;
    }else if(canvasName.Contains("dcaXYInCm")){
        return 3;
    }

    return -1;

}

void SysStudy::vertexZStudy(){
    // function which studies the vertex Z position systematic error


    // loop over all TH1D histograms in histFile
    vector<pair<TH1*, TString>> histograms = GetAllTH1();
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


    cout << "Starting to fill graphs for vertex Z study..." << endl;
    double averageFit = 0.0;
    double averageHist = 0.0;
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
        hist->GetXaxis()->SetTitle("V_{Z} (cm)");
        hist->GetYaxis()->SetTitle("Counts");
        hist->Draw("hist");

        CreateText(Form("Fill number: %d", fillNum), 62);

        // fit the distribution with a Gaussian
        TF1 *gauss = new TF1("gauss", "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        gauss->SetLineColor(kRed);
        gauss->SetLineWidth(2);
        // fit only in range -120 120
        hist->Fit(gauss, "R", "", -120, 120);
        gauss->Draw("same");

        //print the fit results onto the canvas
        CreateText(Form("Mean = %.2f +/- %.2f cm", gauss->GetParameter(1), gauss->GetParError(1)),textFont, 0.72, 0.75, 0.9, 0.72);
        CreateText(Form("Sigma = %.2f +/- %.2f cm", gauss->GetParameter(2), gauss->GetParError(2)), textFont, 0.72, 0.72, 0.9, 0.69);

        DrawSTARInternal();
        
        // integrate the fit function inside +- 100 cm
        double minRange = -100;
        double maxRange = 100;

        double integral = gauss->Integral(minRange, maxRange)/ gauss->Integral(gauss->GetXmin(), gauss->GetXmax());
        CreateText(Form("#epsilon_{fit} = %.2f", integral), textFont, 0.15, 0.85, 0.4, 0.8);
        averageFit += integral;
        
        // integrate the histogram in the same range
        int sum = 0;
        for(int i = 1; i <= hist->GetNbinsX(); ++i) {
            double binCenter = hist->GetBinCenter(i);
            if(binCenter >= minRange && binCenter <= maxRange) {
                sum += hist->GetBinContent(i);
            }
        }

        CreateText(Form("#epsilon_{hist} = %.2f", 1.0*sum/hist->GetEntries()),textFont, 0.15, 0.8, 0.4, 0.75);
        averageHist += 1.0*sum/hist->GetEntries();
        canvas->Update();
        
        canvas->Write(Form("vtxZSysStudy_%d", fillNum));
        canvas->Close();
        
        // add the fill number, mean, sigma, efficiency to the graphs
        efficiencyFitGraph->SetPoint(efficiencyFitGraph->GetN(), fillNum, integral);
        efficiencyHistGraph->SetPoint(efficiencyHistGraph->GetN(), fillNum, 1.0*sum/hist->GetEntries());
        meanGraph->SetPoint(meanGraph->GetN(), fillNum, gauss->GetParameter(1));
        sigmaGraph->SetPoint(sigmaGraph->GetN(), fillNum, gauss->GetParameter(2));
    }
    cout << "Finished filling graphs for vertex Z study..." << endl;
    averageFit /= meanGraph->GetN();
    averageHist /= meanGraph->GetN();

    VzStudyFitEff = averageFit;
    VzStudyHistEff = averageHist;


    CreateCanvas(&canvas, "vtxZSysStudyGraphs", 1200, 800);
    SetGPad();
    canvas->SetName("meanVtxZOfFills");
    meanGraph->GetXaxis()->SetTitle("Fill Number");
    meanGraph->GetYaxis()->SetTitle("<V_{Z}> [cm]");
    meanGraph->GetXaxis()->SetLimits(minFillNum, maxFillNum);
    meanGraph->GetYaxis()->SetRangeUser(-10, 10);
    SetTGraphStyle(meanGraph, kBlue, markerStyleTypical);
    meanGraph->Draw("AP");   
    DrawSTARInternal(); 
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
    DrawSTARInternal();
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

    TLine *averageFitLine = new TLine(minFillNum, averageFit, maxFillNum, averageFit);
    averageFitLine->SetLineColor(kBlue);
    averageFitLine->SetLineStyle(2);
    averageFitLine->SetLineWidth(2);
    TLine *averageHistLine = new TLine(minFillNum, averageHist, maxFillNum, averageHist);
    averageHistLine->SetLineColor(kRed);
    averageHistLine->SetLineStyle(2);
    averageHistLine->SetLineWidth(2);
    averageFitLine->Draw("same");
    averageHistLine->Draw("same");

    CreateLegend(&legend, 0.2, 0.2, 0.45, 0.4);
    legend->AddEntry(efficiencyFitGraph, "Fit Integration", "p");
    legend->AddEntry(averageFitLine, Form("Average Fit: %.2f", averageFit), "l");
    legend->AddEntry(efficiencyHistGraph, "Histogram Integration", "p");
    legend->AddEntry(averageHistLine, Form("Average Hist: %.2f", averageHist), "l");

    legend->Draw("same");

    DrawSTARInternal();
    canvas->Update();
    canvas->Write("efficiencyVtxZ");
    efficiencyFitGraph->Write("efficiencyFitVtxZ_graph");
    efficiencyHistGraph->Write("efficiencyHistVtxZ_graph");
    canvas->Close();
    outFile->cd();


    return;
}


bool SysStudy::runStudy(int VAR,TString condition){

    CreateCanvas(&canvas, "SystematicStudyOfCuts", 1200, 800);
    
    vector<int> yields;
    outFile->cd();
    outFile->mkdir(varNames[VAR]);
    for(int i = 0; i < 3; i++){

        // load signal
        tree->Draw(Form("invMass>>histMass(%d,%f,%f)", nBins, lowerLim, upperLim), condition);
        TH1 *invMassHist = (TH1*)gPad->GetPrimitive("histMass");
        if(!invMassHist || invMassHist->GetEntries() == 0){ 
            cout << "Empty or no hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        TH1D *hSignal = (TH1D*)invMassHist->Clone(Form("hSignal_%d", i));

        // load background
        bcgTree->Draw(Form("invMass>>massBcg(%d,%f,%f)", nBins, lowerLim, upperLim), condition);
        TH1D *hBackground = (TH1D*)gPad->GetPrimitive("massBcg");
        if(!hBackground || hBackground->GetEntries() == 0){
            cout << "Empty or no background hist when trying to run systematic study. Leaving." << endl;
            return false;
        }
        hSignal->Add(hBackground, -1);

        // fit and obtain yield
        FitJPsi *fit = new FitJPsi(hSignal, "poly1");
        fit->fitPeak();
        yields.push_back(fit->getYield());
        fit->saveCanvas(outFile, varNames[VAR] + "_" + variations[i], varNames[VAR]);

    }

    yieldResults[varNames[VAR]] = yields;

    hSysStudyTight->SetBinContent(VAR+1, 1.0*yields[1]/yields[0] );
    hSysStudyLoose->SetBinContent(VAR+1, 1.0*yields[2]/yields[0] );

    return true;

}



void SysStudy::runSysStudy(){

    TString opt = "";
    if(isEmbedding) opt = "embedding ";


    for(int i = 0; i < nVariables; i++){
        if(!runStudy(i, getCondition(opt + varNames[i]))){
            cout << "Error when running systematic study of " << varNames[i] << endl;
        }
    }


    saveSysStudyYieldsHists();

    cout << "Finished with systematic study for all variables" << endl;
    
}

void SysStudy::saveSysStudyYieldsHists(){
    
    
    CreateCanvas(&canvas, "SystematicStudyOfCuts", 1200, 800);
    
    
    canvas->Clear();
    canvas->SetName("SysStudyOverview");
    canvas->cd();
    SetGPad();


    hSysStudyTight->SetMarkerSize(2);
    hSysStudyTight->SetMarkerSize(2);
    hSysStudyTight->Draw("P");
    hSysStudyLoose->Draw("same P");
    
    DrawSTARInternal();
    // draw the nominal line from the min of the histogram to the max
    
    TLine* line = new TLine(hSysStudyTight->GetXaxis()->GetBinCenter(NHITSFIT) - 0.5, 1, hSysStudyTight->GetXaxis()->GetBinCenter(ETA) + 0.5, 1);
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




void SysStudy::peakFittingStudy(){

    TString dir = "SysStudy";
    outFile->mkdir(dir);

    double minYield, maxYield;
    if(noRomanPots){
        minYield = 1200;
        maxYield = 1600;
    }else{
        minYield = 50;
        maxYield = 150;
    }
    
    TH1D *hYield = new TH1D("hYield", "hMeans; Yield ; counts", 50, minYield, maxYield);

    vector<double> lowLim = { 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    vector<double> topLim = { 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};

    for(auto &low : lowLim){
        for(auto &top : topLim){
            for(int nBins = 20; nBins <= 80; ++nBins){
                
                // load signal and background
                tree->Draw(Form("invMass>>hist(%d,%f,%f)", nBins, low, top) , getCondition( isEmbedding ? "embedding" : "" ));
                TH1D *hist = (TH1D*)gPad->GetPrimitive("hist");
                if(!hist || hist->GetEntries() == 0){
                    cout << "Empty or no hist when trying to run systematic study. Leaving." << endl;
                    return;
                }
                TH1D *hSignal = (TH1D*)hist->Clone("hSignal");
                
                if(!isEmbedding){  //embedding has no background
                    bcgTree->Draw(Form("invMass>>massBcg(%d,%f,%f)", nBins, low, top), getCondition());
                    TH1D *hBackground = (TH1D*)gPad->GetPrimitive("massBcg");
                    if(!hBackground || hBackground->GetEntries() == 0){
                        cout << "Empty or no background hist when trying to run systematic study. Leaving." << endl;
                        return;
                    }
                    hSignal->Add(hBackground, -1);
                    delete hBackground;
                }
                FitJPsi *fit = new FitJPsi(hSignal, "poly1");
                fit->fitPeak();
                int yield = fit->getYield();
                fit->saveCanvas(outFile, Form("PeakFittingStudy_%d_%.1f_%.1f", nBins, low, top), dir);

                if(yield < 0){
                    cout << "Something went wrong with peak fitting. Leaving." << endl;
                    return;
                }
                hYield->Fill(yield);
            }
        }
    }


    CreateCanvas(&canvas,"PeakFittingStudy", 1200, 800);
    SetGPad();


    // Fit the histogram with a Gaussian function
    TF1 *gauss = new TF1("gauss", "gaus", minYield, maxYield);
    gauss->SetParameters(30, 1400, 10);
    gauss->SetLineColor(kRed);
    gauss->SetLineWidth(2);
    hYield->Fit(gauss, "R");

    // Draw the results
    canvas->cd();
    hYield->Draw("E hist");
    gauss->Draw("same");

    DrawSTARInternal();

    // create text similarly under which display the mean and sigma of the fit
    TPaveText *res = new TPaveText(0.65, 0.63, 0.88, 0.72, "brNDC");
    res->SetTextColor(kBlack);
    res->SetFillColor(0);
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
    canvas->Write("FittingStudy");

    peakFittingMean = gauss->GetParameter(1);
    peakFittingSigma = gauss->GetParameter(2);
    peakFittingMeanError = gauss->GetParError(1);
    peakFittingSigmaError = gauss->GetParError(2);


}


/*

    auto drawSysPlot = [&](TString canvasName, vector<TH1*> &hists,const vector<TString> &histLegends,TString xAxisDescription, TString outFileDir, bool plotAverage = true)
    {
        CreateCanvas(&canvas,canvasName);
        CreateLegend(&legend, 0.2,0.8,0.45,0.95);

        if(hists.size() != 3 || !hists[0] || !hists[1] || !hists[2]){
            cout << "empty hists in SysStudy::drawSysPlot(), leaving" << endl;
            return false;
        }
        
        
        for (unsigned int i = 0; i < hists.size(); ++i) {


            TH1 *hist = (TH1*) hists[i]->Clone( canvasName + Form("tmp%d_",i));
            hist->Divide(hists[0]);
            SetHistStyle(hist, colorSet[i], markerSet[i]);
            
            if( i == 0){
                hist->GetXaxis()->SetTitle(xAxisDescription);
                hist->GetYaxis()->SetTitle("Ratio to Nominal");
                hist->GetYaxis()->SetRangeUser(0.9, 1.1);
                hist->Draw("hist");
    
            }else{
                hist->Draw("same hist");
            }
    
            legend->AddEntry(hist, histLegends[i],"l");
            if(plotAverage && i!=0)
            { 
                CreateLine(hist->GetXaxis()->GetXmin(),GetAverageYValue(hist),hist->GetXaxis()->GetXmax(),GetAverageYValue(hist));
                line->SetLineWidth( lineWidth );
                line->SetLineStyle(10);
                line->SetLineColor(colorSet[i]);
                line->Draw("same"); 
            }
        }
        
        legend->Draw("same");
        outFile->cd(outFileDir);
        canvas->Write(canvasName);
        canvas->Close();
        
        return true;
    };


    auto drawInvMassPlot = [&](TString canvasName, vector<TH1*> &hists, const vector<TString> &histLegends,TString xAxisDescription, TString outFileDir, bool plotAverage = true){

        CreateCanvas(&canvas,canvasName);
        CreateLegend(&legend,0.2,0.8,0.45,0.95);
        hists[0]->GetYaxis()->SetRangeUser(0.9,1.1);
        for (unsigned int i = 0; i < hists.size(); ++i) 
        {
            TH1 *hist = (TH1*)hists[i]->Clone( canvasName + Form("tmp%i_",i));
    
            hist->Divide(hists[0]);
            SetHistStyle(hist, colorSet[i], markerSet[i]);
            
            if( i == 0){
                hist->GetXaxis()->SetTitle(xAxisDescription);
                hist->GetYaxis()->SetTitle("Ratio to Nominal");

                hist->GetYaxis()->SetRangeUser(yRange[studyDiscriminator(canvasName)].first, yRange[studyDiscriminator(canvasName)].second);
                hist->Draw("hist");
            }else{
                hist->Draw("same hist");
            }
            
            legend->AddEntry(hist, histLegends[i],"l");
            if(plotAverage && i!=0)
            { 
                CreateLine(hist->GetXaxis()->GetXmin(),GetAverageYValue(hist),hist->GetXaxis()->GetXmax(),GetAverageYValue(hist));
                line->SetLineWidth( lineWidth );
                line->SetLineStyle(10);
                line->SetLineColor(colorSet[i]);
                line->Draw("same"); 
            }
        }
        
        legend->Draw("same");
        outFile->cd();
        outFile->cd(outFileDir);
        canvas->Write(canvasName);
        canvas->Close();
        
        return true;
    };
*/

