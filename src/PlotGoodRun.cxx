#include "PlotGoodRun.h"


PlotGoodRun::PlotGoodRun(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotGoodRun::Make() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    cout << "Starting analysis of BEMC tracks from .txt file..." << endl;

    TH1D* hEtaTpc = new TH1D("hEtaTpc", "TPC eta distribution", 100, -1, 1);
    TH1D* hEtaBemc = new TH1D("hEtaBemc", "BEMC eta distribution", 100, -1, 1);
    TH1D* hPhiTpc = new TH1D("hPhiTpc", "TPC phi distribution", 100, -3.14, 3.14);
    TH1D* hPhiBemc = new TH1D("hPhiBemc", "BEMC phi distribution", 100, -3.14, 3.14);

    tree->Draw("tpcEtaAverage>>hist1(100, -1,1)","");
    hEtaTpc->Add( (TH1D*)gPad->GetPrimitive("hist1") );

    tree->Draw("bemcEtaAverage>>hist2(100, -1,1)","");
    hEtaBemc->Add( (TH1D*)gPad->GetPrimitive("hist2") );

    tree->Draw("tpcPhiAverage>>hist3(100, -3.14,3.14)","");
    hPhiTpc->Add( (TH1D*)gPad->GetPrimitive("hist3") );

    tree->Draw("bemcPhiAverage>>hist4(100, -3.14,3.14)","");
    hPhiBemc->Add( (TH1D*)gPad->GetPrimitive("hist4") );

    pair<double,double> bemcFitPars = make_pair(0,0);
    if(hEtaBemc->GetEntries() > 0 ){
        bemcFitPars = FitEtaDistributions(hEtaBemc, 2);
    }else{
        cout << "Empty hEtaBemc histogram" << endl;
    }

    pair<double,double> tpcFitPars = make_pair(0,0);
    if(hEtaTpc->GetEntries() > 0){
        tpcFitPars = FitEtaDistributions(hEtaTpc, 1);
    }else{
        cout << "Empty hEtaTpc histogram" << endl;
    }

    pair<double,double> bemcPhiFitPars = make_pair(0,0);
    if(hPhiBemc->GetEntries() >0){
        bemcPhiFitPars = FitEtaDistributions(hPhiBemc, 3);
    }else{
        cout << "Empty hPhiBemc histogram" << endl;
    }
    
    
    // run over all runs and check if the run is good
    vector<int> goodRunList;
    vector<int> probRetainEventList = loadProbRetainEventList();
    vector<int> badEtaPhi = loadBadEtaPhiRuns("badEtaPhiFile.txt");


    ofstream badEtaPhiFile("badEtaPhiFile_pickedBad.txt");
    //check if the file is opened
    if(!badEtaPhiFile.is_open()){
        cerr << "Could not open file with bad eta phi events. Leaving..." << endl;
        return;
    }

    for(int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);

        hGoodRunListFlow->Fill(1);
        
        
        if(nEventsAll == 0 || atLeast1JPsiTrigger == 0){ // no events were processed or no JPSi trigger
            continue;
        }
        
        hGoodRunListFlow->Fill(2);
        
        if(bemcPhiAverage < bemcPhiFitPars.first - 9*bemcPhiFitPars.second || bemcPhiAverage > bemcPhiFitPars.first + 9*bemcPhiFitPars.second){
            if(find(badEtaPhi.begin(),badEtaPhi.end(),mRunNumber) != badEtaPhi.end()){
                badEtaPhiFile << mRunNumber << endl;
                continue;
            }
        }

        hGoodRunListFlow->Fill(3); 
        
        if(bemcEtaAverage < bemcFitPars.first - 9*bemcFitPars.second || bemcEtaAverage > bemcFitPars.first + 9*bemcFitPars.second){ // outside fit range bemc
            if(find(badEtaPhi.begin(),badEtaPhi.end(),mRunNumber) != badEtaPhi.end()){
                badEtaPhiFile << mRunNumber << endl;
                continue;
            }
        }
        
        hGoodRunListFlow->Fill(4);
        if(find(probRetainEventList.begin(), probRetainEventList.end(), mRunNumber) == probRetainEventList.end()){
            continue;
        }
        
        hGoodRunListFlow->Fill(5);

        //if(RPsClose == 0 ){  // RP not good
        //    continue;
        //}

        hGoodRunListFlow->Fill(6);
        
        goodRunList.push_back(mRunNumber);        
    }
    badEtaPhiFile.close();

    //createAverageTrackHists(runNumbersMapFiltered);
    
    createGoodRunList(goodRunList);
    
    cout << "Finished analysis of BEMC tracks..." << endl;
    
    outFile->cd();
    //hGoodRunListFlow->Write("hGoodRunListFlow");
    TH1General("hGoodRunListFlow", hGoodRunListFlow);
    outFile->Close();
    histFile->Close();

    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    cout << "All good runs saved to goodRuns.list file..." << endl;
    
}

void PlotGoodRun::Init(){
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
	ConnectInputTree(inputPosition, "RunInfo", false);

    
    hGoodRunListFlow = new TH1D("hGoodRunListFlow", "Good run list flow", 6, 1, 7);
    hGoodRunListFlow->GetXaxis()->SetBinLabel(1, "All");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(2, "J/#psi trigger");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(3, "<#eta_{TPC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(4, "<#eta_{BEMC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(5, "p_{retain event}");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(6, "RP in range");
    

    tree->SetBranchAddress("RunNumber", &mRunNumber);
    tree->SetBranchAddress("AtLeast1JPsiTrigger", &atLeast1JPsiTrigger);
    tree->SetBranchAddress("RPsClose", &RPsClose);
    tree->SetBranchAddress("nEventsAll", &nEventsAll);
    tree->SetBranchAddress("nEventsPassed", &nEventsPassed);
    tree->SetBranchAddress("luminosity", &luminosity);
    tree->SetBranchAddress("luminosityError", &luminosityError);
    tree->SetBranchAddress("nTracksBEMC", &nTracksBEMC);
    tree->SetBranchAddress("nClustersBEMC", &nClustersBEMC);
    tree->SetBranchAddress("nTracksTPC", &nTracksTPC);
    tree->SetBranchAddress("nTracksTOF", &nTracksTOF);
    tree->SetBranchAddress("nVertices", &nVertices);
    tree->SetBranchAddress("tpcEtaAverage", &tpcEtaAverage);
    tree->SetBranchAddress("bemcEtaAverage", &bemcEtaAverage);
    tree->SetBranchAddress("tpcPhiAverage", &tpcPhiAverage);
    tree->SetBranchAddress("bemcPhiAverage", &bemcPhiAverage);


    cout << "Finished PlotGoodRun::Init()..." << endl;

}

vector<int> PlotGoodRun::loadBadEtaPhiRuns(TString nameOfFile){

    ifstream file(nameOfFile);
    if(!file.is_open()){
        cout << "No file with name " << nameOfFile << endl;
        return {};
    }

    string line;
    vector<int> badRuns;
    while(getline(file, line)){
        if(line.empty())  continue;

        int run, isGood;
        istringstream iss(line);

        iss >> run >> isGood;

        if(isGood == 0){
            badRuns.push_back(run);
        }
    }
    return badRuns;
}



vector<int> PlotGoodRun::loadProbRetainEventList(){
    // load the prob retain event list from a file
    vector<int> probRetainEventList;
    ifstream input("/star/u/mvranovsk/plots/goodRun17.list");
    if(!input.is_open()){
        cerr << "Could not open file with prob retain event. Leaving..." << endl;
        return {};
    }
    
    string line;
    while(getline(input, line)){
        istringstream iss(line);
        int runNumber;
        double luminosity, probRetainEvent;
        iss >> runNumber >> luminosity >> probRetainEvent;
        probRetainEventList.push_back(runNumber);
    }
    input.close();
    return probRetainEventList;
}



void PlotGoodRun::createAverageTrackHists(vector<pair<array<int, 2>, array<double, 10>>> runNumberMap){

    double minimum = - 0.5;
    double maximum = runNumberMap.size() - 0.5;
    int nBins = runNumberMap.size();

    TH1* hBemcTracksAverage = new TH1D("hBemcTracksAverage", "Number of BEMC tracks per event", nBins, minimum, maximum);
    hBemcTracksAverage->GetXaxis()->SetTitle("run");
    hBemcTracksAverage->GetYaxis()->SetTitle("average number of BEMC tracks per event");

    TH1* hBemcClustersAverage = new TH1D("hBemcClustersAverage", "Number of BEMC clusters per event",nBins, minimum, maximum);
    hBemcClustersAverage->GetXaxis()->SetTitle("run");
    hBemcClustersAverage->GetYaxis()->SetTitle("average number of BEMC clusters per event");

    TH1* hTofTracksAverage = new TH1D("hTofTracksAverage", "Number of Tof tracks per event", nBins, minimum, maximum);
    hTofTracksAverage->GetXaxis()->SetTitle("run");
    hTofTracksAverage->GetYaxis()->SetTitle("average number of TOF tracks per event");

    TH1* hTpcTracksAverage = new TH1D("hTpcTracksAverage", "Number of TPC tracks per event", nBins, minimum, maximum);
    hTpcTracksAverage->GetXaxis()->SetTitle("run");
    hTpcTracksAverage->GetYaxis()->SetTitle("average number of TPC tracks per event");

    TH1* hVerticesAverage = new TH1D("hVerticesAverage", "Number of Vertices per event", nBins, minimum, maximum);
    hVerticesAverage->GetXaxis()->SetTitle("run");
    hVerticesAverage->GetYaxis()->SetTitle("average number of Vertices per event");

    // loop over content in runNumberMap
    int j = 1;
    for(auto runNumber : runNumberMap){
        // get the number of BEMC tracks per event
        double nTracksTpc = runNumber.second[1];
        double nTracksBemc = runNumber.second[2];
        double nClustersBemc = runNumber.second[3];
        double nTracksTof = runNumber.second[4];
        double nVertices = runNumber.second[5];
        //int idx = runNumber.first[0];

        if(nTracksTpc < 0 || nTracksBemc < 0 || nClustersBemc < 0 || nTracksTof < 0 || nVertices < 0){
            cout << "Negative values in line: " << runNumber.first[1] << endl;
            continue;
        }
        hBemcTracksAverage->SetBinContent(j, nTracksBemc);
        hTpcTracksAverage->SetBinContent(j, nTracksTpc);
        hBemcClustersAverage->SetBinContent(j, nClustersBemc);
        hTofTracksAverage->SetBinContent(j, nTracksTof);
        hVerticesAverage->SetBinContent(j, nVertices);
        j++;   

    }


    outFile->cd();

    TH1General("hBemcTracksAverage", hBemcTracksAverage);
    TH1General("hBemcClustersAverage", hBemcClustersAverage);
    TH1General("hTofTracksAverage", hTofTracksAverage);
    TH1General("hTpcTracksAverage", hTpcTracksAverage);
    TH1General("hVerticesAverage", hVerticesAverage);


}

void PlotGoodRun::createGoodRunList(vector<int> goodRunList){

    // get all run numbers from a separate file
    TString allRuns = "/star/u/mvranovsk/star-upcDst/work/lists/fullRun17.list";
    ifstream runNumbersFile(allRuns);
    if(!runNumbersFile.is_open()){
        cerr << "Could not open file with run numbers. Leaving..." << endl;
        return;
    }

    //while loop over all lines 
    string line;
    ofstream goodRunFile("goodRun17Fast.list");
    while(getline(runNumbersFile, line)){

        //get the run number
        int RUNNUMBER = stoi(line.substr(line.find_last_of("/")+1, line.find_last_of(".") - line.find_last_of("/") - 1));

        for(auto run : goodRunList){
            if(run == RUNNUMBER){
                goodRunFile << line << endl;
                break; 
            }
        }
    }
    goodRunFile.close();
    runNumbersFile.close();
    
    cout << "Finished creating fast good run list..." << endl;

    
}


pair<double, double> PlotGoodRun::FitEtaDistributions(TH1D* hist, int SWITCH){
    
    //fit the histogram with a gaussian function
    TF1* gaus = new TF1("gaus", "gaus", -1, 1);
    
    gaus->SetParameters(100, 0, 0.1);
    hist->Fit(gaus, "R");
    
    
    //create a canvas and draw the histogram
    TCanvas* c = new TCanvas(Form("DistFit%d",SWITCH), "c",1500, 1200);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    c->SetLogy();
    hist->SetMarkerStyle(20);
    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
    hist->SetMarkerSize(1);
    if(SWITCH == 1){
        hist->GetXaxis()->SetTitle("<#eta_{TPC}>");
    }else if(SWITCH == 2){
        hist->GetXaxis()->SetTitle("<#eta_{BEMC}>");
    }else if(SWITCH == 3){
        hist->GetXaxis()->SetTitle("<#phi_{BEMC}> [rad]");
    }
    hist->GetYaxis()->SetTitle("counts");
    hist->SetTitle("");
    
    hist->Draw("E HIST");
    
    //draw the fit function
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    gaus->SetLineStyle(1);
    gaus->Draw("SAME");
    
    // print the fit parameters to the canvas
    TPaveText* text = new TPaveText(0.6, 0.55, 0.9, 0.7, "brNDC");
    text->SetFillColorAlpha(kWhite, 0);
    text->SetBorderSize(0);
    text->SetTextSize(0.035);
    text->SetTextFont(82);
    text->SetTextAlign(33);
    text->AddText(Form("Mean: %.4f #pm %.4f", gaus->GetParameter(1), gaus->GetParError(1)));
    text->AddText(Form("Sigma: %.4f #pm %.4f", gaus->GetParameter(2), gaus->GetParError(2)));
    text->Draw("SAME");
    
    gStyle->SetOptFit(0);
    DrawSTARpp510JPsi(0.6,0.82,0.9,0.9, 0.01);
    
    // draw lines at values +- 9 sigma of the fit
    TLine* line1 = new TLine(gaus->GetParameter(1) - 9*gaus->GetParameter(2), 0, gaus->GetParameter(1) - 9*gaus->GetParameter(2), hist->GetMaximum());
    TLine* line2 = new TLine(gaus->GetParameter(1) + 9*gaus->GetParameter(2), 0, gaus->GetParameter(1) + 9*gaus->GetParameter(2), hist->GetMaximum());
    SetLineStyle(line1);
    SetLineStyle(line2);
    line1->Draw("same");
    line2->Draw("same");
    
    
    //save the canvas to the output file
    outFile->cd();
    c->Write(Form("DistributionFit%d", SWITCH));
    
    return make_pair(gaus->GetParameter(1), gaus->GetParameter(2));

    
}


void PlotGoodRun::AnalysisBEMC(vector<int> runNumberList){
    // this function will analyze the number of BEMC tracks per event w.r.t. run number 
    cout << "Starting analysis of BEMC tracks..." << endl;
    // how to access the minimum and maximum of the run number list
    double minimum = - 0.5;
    double maximum = runNumberList.size() - 0.5;
    int nBins = runNumberList.size();

    TH1D* hBemcTracksAverage = new TH1D("hBemcTracksAverage", "Number of BEMC tracks per event", nBins, minimum, maximum);
    hBemcTracksAverage->GetXaxis()->SetTitle("run");
    hBemcTracksAverage->GetYaxis()->SetTitle("average number of BEMC tracks per event");
    
    TH1D* hBemcClustersAverage = new TH1D("hBemcClustersAverage", "Number of BEMC clusters per event",nBins, minimum, maximum);
    hBemcClustersAverage->GetXaxis()->SetTitle("run");
    hBemcClustersAverage->GetYaxis()->SetTitle("average number of BEMC clusters per event");
    
    TH1D* hTofTracksAverage = new TH1D("hTofTracksAverage", "Number of Tof tracks per event", nBins, minimum, maximum);
    hTofTracksAverage->GetXaxis()->SetTitle("run");
    hTofTracksAverage->GetYaxis()->SetTitle("average number of TOF tracks per event");
    // loop over run numbers
    int j = 1;
    ofstream runNumberFile("runNumberToBin.txt");
    cout << "Before loop" << endl;
    for(auto runNumber : runNumberList){
        // get the number of BEMC tracks per event
        tree->Draw("nTracksBemc>>hist1(20,-0.5,19.5)", Form("runNumber == %d", runNumber));
        TH1* hBemcTracks = (TH1*)gDirectory->Get("hist1");
        tree->Draw("nClustersBemc>>hist2(20,-0.5,19.5)", Form("runNumber == %d", runNumber));
        TH1* hBemcClusters = (TH1*)gDirectory->Get("hist2");
        tree->Draw("nTracksTof>>hist3(20,-0.5,19.5)", Form("runNumber == %d", runNumber));
        TH1* hTofTracks = (TH1*)gDirectory->Get("hist3");
        cout << "Obtained the histograms for run number: " << runNumber << endl;

        if(!hBemcTracks || !hBemcClusters || !hTofTracks){
            cout << "Could not get BEMC tracks histogram for run number: " << runNumber << endl;
            j++;
            continue;
        }

        Double_t average = GetMean(hBemcTracks);
        hBemcTracksAverage->SetBinContent(j, average);
        average = GetMean(hBemcClusters);
        hBemcClustersAverage->SetBinContent(j, average);
        average = GetMean(hTofTracks);
        hTofTracksAverage->SetBinContent(j, average);
        
        // save binToRunNumber to a separate .txt file
        runNumberFile << j << " " << runNumber << endl;

        j++;
    }

    runNumberFile.close();


    outFile->cd();
    hBemcTracksAverage->Write("hBemcTracksAverage");
    hBemcClustersAverage->Write("hBemcClustersAverage");
    hTofTracksAverage->Write("hTofTracksAverage");

}

Double_t PlotGoodRun::GetMean(TH1* hist){
    Double_t sum = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++){
        sum += hist->GetBinContent(i) * (i-1);
    }
    return sum/hist->GetEntries();
}




/*
void PlotGoodRun::Make() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    cout << "Starting analysis of BEMC tracks from .txt file..." << endl;
    // how to access the minimum and maximum of the run number list
    ifstream input("average_tracks.txt");
    if(!input.is_open()){
        cerr << "Could not open file with average tracks. Leaving..." << endl;
        return;
    }

    TH1D* hEtaTpc = new TH1D("hEtaTpc", "TPC eta distribution", 100, -1, 1);
    TH1D* hEtaBemc = new TH1D("hEtaBemc", "BEMC eta distribution", 100, -1, 1);
    TH1D* hPhiTpc = new TH1D("hPhiTpc", "TPC phi distribution", 100, -3.14, 3.14);
    TH1D* hPhiBemc = new TH1D("hPhiBemc", "BEMC phi distribution", 100, -3.14, 3.14);

    int nGoodLines = 0;
    string line;        
    vector<pair<array<int, 2>, array<double, 10>>> runNumberMap;
    // load the data from the .txt file
    while(getline(input, line)){
        
        //skip empty lines
        if(line.empty()){
            continue;
        }
        //each line contains either 2 numbers, or 6 numbers. Load them into variables
        istringstream iss(line);
        double nTracksTpc, nTracksBemc, nClustersBemc, nTracksTof, nVertices, etaTpcAverage, etaBemcAverage, phiTpcAverage, phiBemcAverage;
        int j, runnumber, key;
        if(count(line.begin(),line.end(),' ') <=1){
            iss >> j;
            runNumberMap.push_back({{j, 0},{-1.0,0 ,0 ,0, 0,0, 0, 0, 0,0}});
            continue;
        }else if (count(line.begin(), line.end(), ' ') < 3 ) { // means either 0 number of events, or RP is bad
            iss >> j >> runnumber >> key;
            // save it to the map
            runNumberMap.push_back({{j, runnumber},{-1.0*key,0 ,0 ,0, 0,0, 0, 0, 0,0}});
            continue;
        }else if(line.find("-nan") != string::npos){ // check if line contains -nan -> means no good tracks were found in the run
            iss >> j >> runnumber;
            runNumberMap.push_back({{j, runnumber},{-1.0,0 ,0 ,0, 0,0, 0, 0, 0,0}});
            continue;
        }

        iss >> j >> runnumber >> key >> nTracksTpc >> nTracksBemc >> nClustersBemc >> nTracksTof >> nVertices >> etaTpcAverage >> etaBemcAverage >> phiTpcAverage >> phiBemcAverage;
        runNumberMap.push_back({{j, runnumber},{key*1.0,nTracksTpc, nTracksBemc, nClustersBemc, nTracksTof, nVertices, etaTpcAverage, etaBemcAverage, phiTpcAverage, phiBemcAverage}});
        nGoodLines++;
    
        hEtaTpc->Fill(etaTpcAverage);
        hEtaBemc->Fill(etaBemcAverage);
        hPhiTpc->Fill(phiTpcAverage);
        hPhiBemc->Fill(phiBemcAverage);        
    }
    input.close();
    

    pair<double,double> tpcFitPars = FitEtaDistributions(hEtaTpc, 1);
    pair<double,double> bemcFitPars = FitEtaDistributions(hEtaBemc, 2);
    
    
    // run over all runs and check if the run is good
    vector<int> goodRunList;
    vector<pair<array<int, 2>, array<double, 10>>> runNumbersMapFiltered;
    ofstream runNumberToBin("runNumberToBin.txt");
    int idx = 1;

    vector<int> probRetainEventList = loadProbRetainEventList();
    for(auto run : runNumberMap){
        hGoodRunListFlow->Fill(1);
        
        
        hGoodRunListFlow->Fill(2);
        
        if(run.second[0] <= 0 ){ // no events were processed or no JPSi trigger
            continue;
        }
        
        hGoodRunListFlow->Fill(3); 
        
        if(run.second[6] < tpcFitPars.first - 9*tpcFitPars.second || run.second[6] > tpcFitPars.first + 9*tpcFitPars.second){  // outside fit range tpc
            continue;
        }
        hGoodRunListFlow->Fill(4);
        
        if(run.second[7] < bemcFitPars.first - 9*bemcFitPars.second || run.second[7] > bemcFitPars.first + 9*bemcFitPars.second){ // outside fit range bemc
            continue;
        }
        
        hGoodRunListFlow->Fill(5);
        if(find(probRetainEventList.begin(), probRetainEventList.end(), run.first[2]) == probRetainEventList.end()){
            continue;
        }
        
        //if(run.second[0] < 0.5){  // RP not good
        //    continue;
        //}

        hGoodRunListFlow->Fill(6);
        
        goodRunList.push_back(run.first[1]);
        runNumbersMapFiltered.push_back(run);
        runNumberToBin << idx << " " << run.first[1] << " " << run.second[0] << " " <<  run.second[1] << " " <<  run.second[2] << " " << run.second[3] << " " << run.second[4] << " " << run.second[5] << endl;
        idx +=1;
        
    }
    runNumberToBin.close();

    //createAverageTrackHists(runNumbersMapFiltered);
    
    createGoodRunList(goodRunList);
    
    cout << "Finished analysis of BEMC tracks..." << endl;
    
    outFile->cd();
    //hGoodRunListFlow->Write("hGoodRunListFlow");
    TH1General("hGoodRunListFlow", hGoodRunListFlow);
    outFile->Close();
    histFile->Close();

    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    cout << "All good runs saved to goodRuns.list file..." << endl;
    
}
*/