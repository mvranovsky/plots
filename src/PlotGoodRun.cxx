#include "PlotGoodRun.h"


PlotGoodRun::PlotGoodRun(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotGoodRun::PlotGoodRun(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotGoodRun::Make() {

    createAverageTracksPlots();

    vector<int> probRetainEventList = getPREList();
    //cout << "Size of Prob Retain event list in PlotGoodRun: " << probRetainEventList.size() << endl;


    vector<int> badEtaPhi = loadBadEtaPhiRuns("badEtaPhiFile.txt");

    goodRunList.clear();
    rpGoodRunList.clear();

    cout << "Starting with run loop in PlotGoodRun" << endl;
    setBranchAddresses();
    for(int i = 0; i < tree->GetEntries(); i++){

        tree->GetEntry(i);

        if(mRunNumber == 0){
            continue;
        }

        if(nEventsJPsi < 0 || nEventsJPsi > 1e+9){
            cout << "bad nEventsJPsi " << nEventsJPsi << " for run: " << mRunNumber << endl;
        }

        hGoodRunListFlow->Fill(1);

        if( atLeast1JPsiTrigger == 0)  continue;

        hGoodRunListFlow->Fill(2);
        
        if(!inRangePhiBemc(bemcPhiAverage)) continue;
        
        hGoodRunListFlow->Fill(3);
        
        if(!inRangeEtaBemc(bemcEtaAverage)) continue;

        //check if it is in badeta phi runs
        if(find(badEtaPhi.begin(), badEtaPhi.end(), mRunNumber) != badEtaPhi.end()){
            continue;
        }

        hGoodRunListFlow->Fill(4);

        if(find(probRetainEventList.begin(), probRetainEventList.end(), mRunNumber) == probRetainEventList.end()){
            continue;
        }
        
        hGoodRunListFlow->Fill(5);
        
        goodRunList.push_back(mRunNumber);    
    
        if(RPsClose == 0 ){  
            continue;   
        } 
        rpGoodRunList.push_back(mRunNumber);

        hGoodRunListFlow->Fill(6);
    }

    cout << "Good run list without RPs size: " << goodRunList.size() << endl;

    cout << "Good run list with RPs size: " << rpGoodRunList.size() << endl;

    
    createGoodRunList(goodRunList, "goodRunNoRP.list");
    createGoodRunList(rpGoodRunList, "goodRunWithRP.list");



    cout << "Corrected luminosity for analysis with RP: " << mProbRetainEvent->getLuminosity(rpGoodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB()) << " pb^-1" << endl;
    cout << "Corrected luminosity for analysis without RP: " << mProbRetainEvent->getLuminosity(goodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB()) << " pb^-1" << endl;


    outFile->cd();

    handleHistograms(nameOfAnaGoodRunDir, "GoodRunData");

    //hGoodRunListFlow->Write("hGoodRunListFlow");
    TH1General("hGoodRunListFlow", hGoodRunListFlow);
}

void PlotGoodRun::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();

    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
}

void PlotGoodRun::Init(){
	//define the output file which will store all the canvases
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}

    outFile->cd();
    outFile->mkdir(nameOfAnaGoodRunDir);
    outFile->cd();
    
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );
    
    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaGoodRunTree, tree, bcgTree);

    //ConnectInputTree(inputPosition,nameOfZBTree, ZBTree, bcgTree);

    hGoodRunListFlow = new TH1D("hGoodRunListFlow", "Good run list flow", 6, 1, 7);
    hGoodRunListFlow->GetXaxis()->SetBinLabel(1, "All");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(2, "J/#psi trigger");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(3, "<#phi_{BEMC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(4, "<#eta_{BEMC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(5, "p_{retain event}");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(6, "RP in range");

    
    cout << "Finished PlotGoodRun::Init()..." << endl;
}
void PlotGoodRun::setBranchAddresses(){

    if(!tree){
        cerr << "Tree not loaded" << endl;
        return;
    }
    tree->SetBranchAddress("runNumber", &mRunNumber);
    tree->SetBranchAddress("AtLeast1JPsiTrigger", &atLeast1JPsiTrigger);
    tree->SetBranchAddress("RPsClose", &RPsClose);
    tree->SetBranchAddress("nEventsLumiFile", &nEventsLumiFile);
    tree->SetBranchAddress("nEventsJPsi", &nEventsJPsi);
    tree->SetBranchAddress("nEventsZBVetoAll", &nEventsZBVetoAll);
    tree->SetBranchAddress("nEventsZBVetoPassed", &nEventsZBVetoPassed);
    tree->SetBranchAddress("nEventsTEAll", &nEventsTEAll);
    tree->SetBranchAddress("nEventsTEPassed", &nEventsTEPassed);
    tree->SetBranchAddress("instLumi", &instLumi);

    tree->SetBranchAddress("nTracksBEMC", &nTracksBEMC);
    tree->SetBranchAddress("nClustersBEMC", &nClustersBEMC);
    tree->SetBranchAddress("nTracksTPC", &nTracksTPC);
    tree->SetBranchAddress("nTracksTOF", &nTracksTOF);
    tree->SetBranchAddress("nVertices", &nVertices);
    tree->SetBranchAddress("tpcEtaAverage", &tpcEtaAverage);
    tree->SetBranchAddress("bemcEtaAverage", &bemcEtaAverage);
    tree->SetBranchAddress("tpcPhiAverage", &tpcPhiAverage);
    tree->SetBranchAddress("bemcPhiAverage", &bemcPhiAverage);
}

vector<int> PlotGoodRun::getPREList(){ // this function has to exist because i need class ProbRetainEvent closed, so getting entries from tree works as i want it to
    // run over all runs and check if the probability of retaining event is in range
    mProbRetainEvent = new ProbRetainEvent(outFile, tree);

    cout << "Running ProbRetainEvent" << endl;
    mProbRetainEvent->Make();
    vector<int> probRetainEventList = mProbRetainEvent->getGoodRunList();

    if(probRetainEventList.empty()){
        cout << "No good runs found. Exiting..." << endl;
        return {};
    }
    return probRetainEventList;
    
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


/*

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
    
    TH1General("hBemcTracksAverage", hBemcTracksAverage, nameOfAnaGoodRunDir);
    TH1General("hBemcClustersAverage", hBemcClustersAverage, nameOfAnaGoodRunDir);
    TH1General("hTofTracksAverage", hTofTracksAverage, nameOfAnaGoodRunDir);
    TH1General("hTpcTracksAverage", hTpcTracksAverage, nameOfAnaGoodRunDir);
    TH1General("hVerticesAverage", hVerticesAverage, nameOfAnaGoodRunDir);
    
    
}
*/

void PlotGoodRun::createGoodRunList(vector<int> goodRunList, TString outName){

    // get all run numbers from a separate file
    TString allRuns = "/star/u/mvranovsk/star-upcDst/work/lists/fullRun17.list";
    ifstream runNumbersFile(allRuns);
    if(!runNumbersFile.is_open()){
        cerr << "Could not open file with run numbers. Leaving..." << endl;
        return;
    }

    //while loop over all lines 
    string line;
    ofstream goodRunFile(outName);
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
    
    TString name[3] = {"averageEtaTPC", "averageEtaBEMC", "averagePhiBEMC"};
    Double_t topLine[3] = {100,100,20};
    //create a canvas and draw the histogram
    TCanvas* c = new TCanvas(name[SWITCH-1], name[SWITCH-1],1500, 1200);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
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
    
    c->SetLogy();
    hist->Draw("E HIST");
    
    //draw the fit function
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    gaus->SetLineStyle(1);
    gaus->Draw("SAME");
    
    // print the fit parameters to the canvas
    TPaveText* text = new TPaveText(0.15, 0.7, 0.4, 0.85, "brNDC");
    text->SetFillColorAlpha(kWhite, 0);
    text->SetBorderSize(0);
    text->SetTextSize(0.035);
    text->SetTextFont(82);
    text->SetTextAlign(13);
    text->AddText(Form("Mean: %.4f #pm %.4f", gaus->GetParameter(1), gaus->GetParError(1)));
    text->AddText(Form("Sigma: %.4f #pm %.4f", gaus->GetParameter(2), gaus->GetParError(2)));
    text->Draw("SAME");
    
    gStyle->SetOptFit(0);
    DrawSTARInternal(0.65,0.87,0.95,0.95);
    
    // draw lines at values +- 9 sigma of the fit
    TLine* line1 = new TLine(gaus->GetParameter(1) - 9*gaus->GetParameter(2), 0, gaus->GetParameter(1) - 9*gaus->GetParameter(2), topLine[SWITCH-1]);
    TLine* line2 = new TLine(gaus->GetParameter(1) + 9*gaus->GetParameter(2), 0, gaus->GetParameter(1) + 9*gaus->GetParameter(2), topLine[SWITCH-1]);
    SetLineStyle(line1);
    SetLineStyle(line2);

    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    
    line1->Draw("same");
    line2->Draw("same");
    
    c->SetLogy();
    c->Modified();
    c->Update();
    
    //save the canvas to the output file
    outFile->cd();
    outFile->cd(nameOfAnaGoodRunDir);
    c->Write(name[SWITCH-1]);
    
    return make_pair(gaus->GetParameter(1), gaus->GetParameter(2));

    
}



Double_t PlotGoodRun::GetMean(TH1* hist){
    Double_t sum = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++){
        sum += hist->GetBinContent(i) * (i-1);
    }
    return sum/hist->GetEntries();
}



void PlotGoodRun::createAverageTracksPlots(){

    
    TH1D* hEtaTpc = new TH1D("hEtaTpc", "TPC eta distribution", 100, -0.5, 0.5);
    TH1D* hEtaBemc = new TH1D("hEtaBemc", "BEMC eta distribution", 100, -0.5, 0.5);
    TH1D* hPhiTpc = new TH1D("hPhiTpc", "TPC phi distribution", 100, -0.5, 0.5);
    TH1D* hPhiBemc = new TH1D("hPhiBemc", "BEMC phi distribution", 100, -0.5, 0.5);

    tree->Draw("tpcEtaAverage>>hist1(100, -0.5,0.5)","");
    hEtaTpc->Add( (TH1D*)gPad->GetPrimitive("hist1") );

    tree->Draw("bemcEtaAverage>>hist2(100, -0.5,0.5)","");
    hEtaBemc->Add( (TH1D*)gPad->GetPrimitive("hist2") );
    
    tree->Draw("tpcPhiAverage>>hist3(100, -0.5,0.5)","");
    hPhiTpc->Add( (TH1D*)gPad->GetPrimitive("hist3") );

    tree->Draw("bemcPhiAverage>>hist4(100, -0.5,0.5)","");
    hPhiBemc->Add( (TH1D*)gPad->GetPrimitive("hist4") );
    
    bemcFitPars = make_pair(0,0);
    if(hEtaBemc->GetEntries() > 0 ){
        bemcFitPars = FitEtaDistributions(hEtaBemc, 2);
    }else{
        cout << "Empty hEtaBemc histogram" << endl;
    }

    tpcFitPars = make_pair(0,0);
    if(hEtaTpc->GetEntries() > 0){
        tpcFitPars = FitEtaDistributions(hEtaTpc, 1);
    }else{
        cout << "Empty hEtaTpc histogram" << endl;
    }

    bemcPhiFitPars = make_pair(0,0);
    if(hPhiBemc->GetEntries() >0){
        bemcPhiFitPars = FitEtaDistributions(hPhiBemc, 3);
    }else{
        cout << "Empty hPhiBemc histogram" << endl;
    }
    
}

bool PlotGoodRun::inRangePhiBemc(double average){

    if(average < bemcPhiFitPars.first - 9*bemcPhiFitPars.second || average > bemcPhiFitPars.first + 9*bemcPhiFitPars.second){
        return false;
    }else{
        return true;
    }
}

bool PlotGoodRun::inRangeEtaBemc(double average){
    
    if(average < bemcFitPars.first - 9*bemcFitPars.second || average > bemcFitPars.first + 9*bemcFitPars.second){ // outside fit range bemc
        return false;
    }else{
        return true;
    }
}