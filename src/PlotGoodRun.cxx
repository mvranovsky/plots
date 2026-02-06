#include "PlotGoodRun.h"


PlotGoodRun::PlotGoodRun(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

PlotGoodRun::PlotGoodRun(const string mInputList, shared_ptr<TFile> file): Plot(mInputList, file) {}

void PlotGoodRun::Make() {

    createAverageTracksPlots();

    vector<int> probRetainEventList = getPREList();
    if(DEBUGMODE) cout << "Size of Prob Retain event list in PlotGoodRun: " << probRetainEventList.size() << endl;

    vector<int> badEtaPhi = loadBadEtaPhiRuns("badEtaPhiFile.txt");

    vector<int> badRunsAll, badRunsRP;
    goodRunList.clear();
    rpGoodRunList.clear();
    goodRunTrigger0.clear();
    goodRunTrigger1.clear();
    goodRunTrigger2.clear();

    setBranchAddresses();
    Long64_t lostPREvents = 0;

    if(DEBUGMODE) cout << "Starting with run loop in PlotGoodRun" << endl;
    for(int i = 0; i < tree->GetEntries(); i++){

        tree->GetEntry(i);

        if(mRunNumber == 0){
            continue;
        }

        if(nEventsJPsi < 0 || nEventsJPsi > 1e+6){
            if(DEBUGMODE) cout << "bad nEventsJPsi " << nEventsJPsi << " for run: " << mRunNumber << endl;
            badRunsAll.push_back(mRunNumber);
            continue;
        }

        hGoodRunListFlow->Fill(1);
        hNEventsJPsiFlow->SetBinContent(1, hNEventsJPsiFlow->GetBinContent(1) + nEventsJPsi);

        if(atLeast1JPsiTrigger == 0 && nEventsJPsi > 0){
            if(DEBUGMODE) cout << "Inconsistent JPsi trigger(" << atLeast1JPsiTrigger << ") and nEventsJPsi (" <<  nEventsJPsi <<") for run: " << mRunNumber << endl;
            badRunsAll.push_back(mRunNumber);
            continue;
        }

        hGoodRunListFlow->Fill(2);
        hNEventsJPsiFlow->SetBinContent(2, hNEventsJPsiFlow->GetBinContent(2) + nEventsJPsi);

        
        if(!inRangePhiBemc(bemcPhiAverage)){
            badRunsAll.push_back(mRunNumber);
            continue;
        }

        hGoodRunListFlow->Fill(3);
        hNEventsJPsiFlow->SetBinContent(3, hNEventsJPsiFlow->GetBinContent(3) + nEventsJPsi);

        if(!inRangeEtaBemc(bemcEtaAverage)){
            badRunsAll.push_back(mRunNumber);
            continue;
        }

        //check if it is in badeta phi runs
        if(find(badEtaPhi.begin(), badEtaPhi.end(), mRunNumber) != badEtaPhi.end()){
            badRunsAll.push_back(mRunNumber);
            continue;
        }

        hGoodRunListFlow->Fill(4);
        hNEventsJPsiFlow->SetBinContent(4, hNEventsJPsiFlow->GetBinContent(4) + nEventsJPsi);

        
        if(find(probRetainEventList.begin(), probRetainEventList.end(), mRunNumber) == probRetainEventList.end()){
            if(DEBUGMODE) cout << "Run " << mRunNumber << " not in p_retain event list. Number of JPsi events: " << nEventsJPsi << endl;
            lostPREvents +=nEventsJPsi;
            badRunsAll.push_back(mRunNumber);
            continue;
        }
        
        
        hGoodRunListFlow->Fill(5);
        hNEventsJPsiFlow->SetBinContent(5, hNEventsJPsiFlow->GetBinContent(5) + nEventsJPsi);
        
        goodRunList.push_back(mRunNumber);  
        
        if(nEventsJPsiTrigger0 > 0){
            goodRunTrigger0.push_back(mRunNumber);
        }
        if(nEventsJPsiTrigger1 > 0){
            goodRunTrigger1.push_back(mRunNumber);
        }
        if(nEventsJPsiTrigger2 > 0){
            goodRunTrigger2.push_back(mRunNumber);
        }

    
        if(RPsClose == 0 ){  
            badRunsRP.push_back(mRunNumber);
            continue;   
        } 
        rpGoodRunList.push_back(mRunNumber);

        hGoodRunListFlow->Fill(6);
        hNEventsJPsiFlow->SetBinContent(6, hNEventsJPsiFlow->GetBinContent(6) + nEventsJPsi);
    }
    if(DEBUGMODE) cout << "Finished run loop in PlotGoodRun" << endl;

    cout << "Good run list without RPs size: " << goodRunList.size() << endl;

    cout << "Good run list with RPs size: " << rpGoodRunList.size() << endl;

    if(DEBUGMODE) cout << "Lost JPsi events due to p_retain cut: " << lostPREvents << endl;

    if(DEBUGMODE) cout << "Creating good run list files..." << endl;
    createGoodRunList(goodRunList, "goodRunNoRP.list");

    createGoodRunList(rpGoodRunList, "goodRunWithRP.list");

    createGoodRunList(goodRunTrigger0, "goodRunTrigger570209.list");
    createGoodRunList(goodRunTrigger1, "goodRunTrigger570219.list");
    createGoodRunList(goodRunTrigger2, "goodRunTrigger570229.list");


    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Luminosity without the correction for analysis with RP: " << mProbRetainEvent->getLuminosity(rpGoodRunList, 1,0, 1,0) << " pb^-1" << endl;
    cout << "Corrected luminosity for analysis with RP: " << mProbRetainEvent->getLuminosity(rpGoodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB(), 1,0) << " pb^-1" << endl;
    cout << "Luminosity with the correction for analysis with RP corrected also to pile-up: " << mProbRetainEvent->getLuminosity(rpGoodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB(), mProbRetainEvent->getPileUpCorrectionA(), mProbRetainEvent->getPileUpCorrectionB()) << " pb^-1" << endl;
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Luminosity without the correction for analysis without RP: " << mProbRetainEvent->getLuminosity(goodRunList, 1,0, 1,0) << " pb^-1" << endl;
    cout << "Corrected luminosity for analysis without RP: " << mProbRetainEvent->getLuminosity(goodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB(), 1,0) << " pb^-1" << endl;
    cout << "Corrected luminosity for analysis without RP corrected also to pile-up: " << mProbRetainEvent->getLuminosity(goodRunList, mProbRetainEvent->getA(), mProbRetainEvent->getB(), mProbRetainEvent->getPileUpCorrectionA(), mProbRetainEvent->getPileUpCorrectionB()) << " pb^-1" << endl;
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;

    printBadRunsToFile(badRunsAll, "badRunsAll.list");

    printBadRunsToFile(badRunsRP, "badRunsRP.list");

    outFile->cd();

    handleHistograms(nameOfAnaGoodRunDir, "GoodRunData");

    //hGoodRunListFlow->Write("hGoodRunListFlow");
    TH1General("hGoodRunListFlow", hGoodRunListFlow);
    TH1General("hNEventsJPsiFlow", hNEventsJPsiFlow);
}

void PlotGoodRun::Finish(){

    if(outFile) outFile->Close();
    if(histFile) histFile->Close();

    if(DEBUGMODE) cout << "All histograms successfully saved to canvases..." << endl;
    if(DEBUGMODE)cout << "The output file is saved: " << outputPosition << endl;
}

void PlotGoodRun::Init(){
	//define the output file which will store all the canvases
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	if(!outFile || outFile->IsZombie() ){
		cerr << "ERROR in PlotGoodRun::Init(): Couldn't open output file. Leaving..." << endl;
        return;
	}

    outFile->cd();
    outFile->mkdir(nameOfAnaGoodRunDir);
    outFile->cd();
    
    if(DEBUGMODE) cout << "Opening histogram file in PlotGoodRun::Init()..." << endl;
    histFile = shared_ptr<TFile>( new TFile("histFile.root", "read") );
    
    if(!histFile || histFile->IsZombie() || !histFile->IsOpen()){
        cerr << "ERROR in PlotGoodRun::Init(): Could not get file with histograms. Leaving..." << endl;
        return;
    }
    
    if(DEBUGMODE) cout << "Connecting input tree for good run..." << endl;
    ConnectInputTree(inputPosition, nameOfAnaGoodRunTree, tree, bcgTree);
    if(!tree){
    	cerr << "ERROR in PlotGoodRun::Init(): Couldn't open tree with data. Returning." << endl;
    	return;
    }
    if(DEBUGMODE) cout << "Successfully connected input tree for good run." << endl;


    hGoodRunListFlow = new TH1D("hGoodRunListFlow", "Good run list flow", 6, 1, 7);
    hGoodRunListFlow->GetXaxis()->SetBinLabel(1, "All");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(2, "J/#psi trigger");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(3, "<#phi_{BEMC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(4, "<#eta_{BEMC}>");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(5, "p_{retain event}");
    hGoodRunListFlow->GetXaxis()->SetBinLabel(6, "RP in range");

    hNEventsJPsiFlow = new TH1D("hNEventsJPsiFlow", "NEvents JPsi flow", 6, 1, 7);
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(1, "All");
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(2, "J/#psi trigger");
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(3, "<#phi_{BEMC}>");
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(4, "<#eta_{BEMC}>");
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(5, "p_{retain event}");
    hNEventsJPsiFlow->GetXaxis()->SetBinLabel(6, "RP in range");

    if(DEBUGMODE) cout << "Finished PlotGoodRun::Init()..." << endl;
}
void PlotGoodRun::setBranchAddresses(){

    if(!tree){
        cerr << "ERROR in PlotGoodRun::setBranchAddresses(): Tree not loaded" << endl;
        return;
    }

    tree->SetBranchAddress("runNumber", &mRunNumber);
    tree->SetBranchAddress("AtLeast1JPsiTrigger", &atLeast1JPsiTrigger);
    tree->SetBranchAddress("RPsClose", &RPsClose);
    tree->SetBranchAddress("nEventsLumiFile", &nEventsLumiFile);
    tree->SetBranchAddress("nEventsJPsi", &nEventsJPsi);
    tree->SetBranchAddress("nEventsJPsiTrigger0", &nEventsJPsiTrigger0);
    tree->SetBranchAddress("nEventsJPsiTrigger1", &nEventsJPsiTrigger1);
    tree->SetBranchAddress("nEventsJPsiTrigger2", &nEventsJPsiTrigger2);
    tree->SetBranchAddress("nEventsZBVetoAll", &nEventsZBVetoAll);
    tree->SetBranchAddress("nEventsZBVetoPassed", &nEventsZBVetoPassed);
    tree->SetBranchAddress("instLumi", &instLumi);


    tree->SetBranchAddress("tpcEtaAverage", &tpcEtaAverage);
    tree->SetBranchAddress("bemcEtaAverage", &bemcEtaAverage);
    tree->SetBranchAddress("tpcPhiAverage", &tpcPhiAverage);
    tree->SetBranchAddress("bemcPhiAverage", &bemcPhiAverage);
}

vector<int> PlotGoodRun::getPREList(){ // this function has to exist because i need class ProbRetainEvent closed, so getting entries from tree works as i want it to
    
    // run over all runs and check if the probability of retaining event is in range
    mProbRetainEvent = new ProbRetainEvent(outFile, tree);

    if(DEBUGMODE) cout << "Running ProbRetainEvent" << endl;
    mProbRetainEvent->Make("linear");
    vector<int> probRetainEventList = mProbRetainEvent->getGoodRunList();

    if(probRetainEventList.empty()){
        cerr << "ERROR in PlotGoodRun::getPREList(): No good runs found. Exiting..." << endl;
        return {};
    }
    if(DEBUGMODE) cout << "Finished running ProbRetainEvent" << endl;

    // create also plots for tof mult veto and bbc e/w veto
    //mProbRetainEvent->plotTofMultVeto("linear");
    mProbRetainEvent->plotBbcEVeto("linear");
    mProbRetainEvent->plotBbcWVeto("linear");

    return probRetainEventList;
}


vector<int> PlotGoodRun::loadBadEtaPhiRuns(TString nameOfFile){

    ifstream file(nameOfFile);
    if(!file.is_open()){
        cerr << "ERROR in PlotGoodRun::loadBadEtaPhiRuns(): No file with name " << nameOfFile << endl;
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


void PlotGoodRun::printBadRunsToFile(vector<int> badRuns, TString outName){

    ofstream outFile(outName);
    if(!outFile.is_open()){
        cerr << "ERROR in PlotGoodRun::printBadRunsToFile(): Could not open output file to print bad runs. Leaving..." << endl;
        return;
    }

    for(auto run : badRuns){
        outFile << run << endl;
    }

    outFile.close();
}



void PlotGoodRun::createGoodRunList(vector<int> goodRunList, TString outName){

    // get all run numbers from a separate file
    TString allRuns = "/star/u/mvranovsk/star-upcDst/work/lists/fullRun17.list";
    ifstream runNumbersFile(allRuns);
    if(!runNumbersFile.is_open()){
        cerr << "ERROR in PlotGoodRun::createGoodRunList(): Could not open file with run numbers. Leaving..." << endl;
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
        cerr << "ERROR in PlotGoodRun::createAverageTracksPlots(): Empty hEtaBemc histogram" << endl;
    }

    tpcFitPars = make_pair(0,0);
    if(hEtaTpc->GetEntries() > 0){
        tpcFitPars = FitEtaDistributions(hEtaTpc, 1);
    }else{
        cerr << "ERROR in PlotGoodRun::createAverageTracksPlots(): Empty hEtaTpc histogram" << endl;
    }

    bemcPhiFitPars = make_pair(0,0);
    if(hPhiBemc->GetEntries() >0){
        bemcPhiFitPars = FitEtaDistributions(hPhiBemc, 3);
    }else{
        cerr << "ERROR in PlotGoodRun::createAverageTracksPlots(): Empty hPhiBemc histogram" << endl;
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

// functions below are mainly for running from CrossSectionMaker

void PlotGoodRun::loadGoodRunList(TString fileName, TString option){
    
    ifstream file(fileName);
    if(!file.is_open()){
        
        cerr << "ERROR in PlotGoodRun::loadGoodRunList(): Could not open good run list file. Leaving..." << endl;
        return;
    }

    string line;
    goodRunList.clear();
    rpGoodRunList.clear();
    while(getline(file, line)){
        if(line.empty())  continue;

        int runNumber = stoi(line.substr(line.find_last_of("/")+1, line.find_last_of(".") - line.find_last_of("/") - 1));
        goodRunList.push_back(runNumber);
        rpGoodRunList.push_back(runNumber);
    }

    if(DEBUGMODE) cout << "Initializing ProbRetainEvent with custom good run list..." << endl;
    mProbRetainEvent = new ProbRetainEvent(outFile, tree);
    if(DEBUGMODE) cout << "Running ProbRetainEvent with custom good run list..." << endl;
    mProbRetainEvent->runCustomList(goodRunList, option);
    
}