#include "CrossSectionMaker.h"

CrossSectionMaker::CrossSectionMaker(const TString AnaDir, const TString EmbeddingDir, const TString GoodRunDir, const TString ZeroBiasDir, TString outputFile) {

    mAnaDir = AnaDir;
    mEmbeddingDir = EmbeddingDir;
    mGoodRunDir = GoodRunDir;
    mZeroBiasDir = ZeroBiasDir;
    mOutFileName = outputFile;

    // initialize utility
    mUtil = new Util();

    // connect output file
    mOutFile = shared_ptr<TFile>(new TFile( mOutFileName.Data(), "RECREATE") );
    if(!mOutFile || mOutFile->IsZombie()){
        cerr << "ERROR: Could not open output file" << endl;
    }

    // initialize analysis class for data
    if(DEBUGMODE) cout << "Initiating analysis class..." << endl;
    mAna = new PlotAnaJPsi(AnaDir.Data(), mOutFile);
    mAna->Init();

    isRomanPots = mAna->isRomanPots();
    if(DEBUGMODE) cout << "This analysis is with RPs: " << isRomanPots << endl;


    // initialize embedding class
    if(DEBUGMODE) cout << "Initializing embedding class..." << endl;
    mEmbedding = new PlotEmbeddingJPsi(EmbeddingDir.Data(), mOutFile);
    mEmbedding->Init();

    // initialize goodrun class
    if(DEBUGMODE) cout << "Initializing good run class..." << endl;
    mGoodRun = new PlotGoodRun(GoodRunDir.Data(), mOutFile);
    mGoodRun->Init();

    if(DEBUGMODE) cout << "Initializing zero bias class..." << endl;
    mZeroBias = new PlotZeroBias(ZeroBiasDir.Data(), mOutFile);
    mZeroBias->Init();

    if(DEBUGMODE) cout << "Initializing ProbRetainEvent class..." << endl;
    TTree *mGoodRunTreeBcg = nullptr;
    mGoodRun->ConnectInputTree(GoodRunDir.Data(), nameOfAnaGoodRunTree.Data(), mGoodRunTree, mGoodRunTreeBcg);

    mRpSysError.val = 0.889;
    mRpSysError.errLow = 0.017;
    mRpSysError.errTop = 0.017;

}

CrossSectionMaker::~CrossSectionMaker(){
    mAna->Finish();
    mEmbedding->Finish();
    mGoodRun->Finish();
    mZeroBias->Finish();
    embMin->Finish();
    embMax->Finish();

    if(mOutFile) mOutFile->Close();

    if(mEmbedding) delete mEmbedding;
    if(mAna) delete mAna;
    if(mGoodRun) delete mGoodRun;
}

void CrossSectionMaker::Make(){


    if(DEBUGMODE) cout << "Starting CrossSectionMaker::Make()..." << endl;

    if(DEBUGMODE) cout << "Fitting J/Psi peak..." << endl;
    mAna->fitPeakJPsi("", "JPsiPeakFit", 0);

    if(!isRomanPots){

        if(DEBUGMODE) cout << "Running data systematic study..." << endl;
        //mAna->runSysStudy();

        if(DEBUGMODE) cout << "Running vertex Z study..." << endl;
        //mAna->vertexZStudy();
        
        if(DEBUGMODE) cout << "Running peak fitting study..." << endl;
        //mAna->peakFittingStudy();
    }

    if(DEBUGMODE) cout << "Making embedding plots..." << endl;
    mEmbedding->fitPeakJPsi();

    if(!isRomanPots){
        if(DEBUGMODE) cout << "Running embedding systematic study..." << endl;
        //mEmbedding->runSysStudy();
    }

    if(DEBUGMODE) cout << "Loading good run list..." << endl;
    if(isRomanPots){
        mGoodRun->loadGoodRunList(goodRunListWithRp, "linear");
    }else{
        mGoodRun->loadGoodRunList(goodRunListNoRp,"linear");
    }

    if(DEBUGMODE) cout << "Calculating trigger topology efficiency..." << endl;

    mZeroBias->calculateTriggerEfficiency("mee");

    if(DEBUGMODE) cout << "-------------------------------------------------------------------------------------------------" << endl;
    if(DEBUGMODE) cout << "-------------------------------------------------------------------------------------------------" << endl;
    if(DEBUGMODE) cout << "-------------------------------------------------------------------------------------------------" << endl;
    if(DEBUGMODE) cout << "Running analysis with Roman Pots: " << isRomanPots << endl;
    
    Value corrYield = getYieldCorrection(30, 0, 1.5, mAna->getCondition(), mEmbedding);
    mCorrectedYield = corrYield.val;
    if(DEBUGMODE) cout << "Corrected yield: " << mCorrectedYield << endl;

    if(corrYield.val <= 0){
        cerr << "ERROR in CrossSectionMaker::Make(): Could not get yield correction. Leaving..." << endl;
        return;
    }

    if(DEBUGMODE) cout << "Yield correction: " << corrYield.val << " + " << corrYield.errTop << " - " << corrYield.errLow << endl;
    
    saveYieldCorrection("YieldCorrectionMain");
    setYieldCorrFinal( corrYield.val);
    setYieldCorrErrTopFinal( corrYield.errTop );
    setYieldCorrErrLowFinal( corrYield.errLow );
    
    if(DEBUGMODE) cout << "Calculating cross section..." << endl;
    // calculate cross section
    double cross_section = CalculateCrossSection();

    if(DEBUGMODE) cout << "Calculating systematic errors..." << endl;
    
    pair<double,double> systematics = SysError();

    if(systematics.first == -1 || systematics.second == -1){
        cerr << "ERROR in CrossSectionMaker::Make(): One of the systematic errors is -1. Leaving..." << endl;
        return;
    }
    
    cout << "Final result: " << cross_section << " pb" << endl;
    cout << "Error low: " << systematics.first*cross_section << " pb" << endl;
    cout << "Error high: " << systematics.second*cross_section << " pb" << endl;
        

    if(DEBUGMODE) cout << "Finished CrossSectionMaker::Make()" << endl;

}

double CrossSectionMaker::CalculateCrossSection(){
    
    // get yields from data
    double yield = getYieldCorrFinal();
    
    double yieldGlobal = mAna->getYieldFinal()/mEmbedding->getEfficiencyFinal();
    //cout << "Yield corrected by pT spectrum of embedding: " << yield << endl;
    //cout << "Yield corrected globally: " << yieldGlobal << endl;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    if(DEBUGMODE) cout << "Trigger eff: " << triggerEff << endl;
    efficiencies.push_back( triggerEff );

    double RPEff = getRPEfficiency();
    if(DEBUGMODE) cout << "RP eff: " << RPEff << endl;
    efficiencies.push_back( RPEff );

    double vertexZEff = mAna->getVzStudyAverageEff();  // this cut is needed because embedding does not include this cut
    if(DEBUGMODE) cout << "Vz efficiency: " << vertexZEff << endl;
    efficiencies.push_back( vertexZEff );

    double vertexRecoEff = getVertexRecoEfficiency();
    if(DEBUGMODE) cout << "Vertex reco efficiency: " << vertexRecoEff << endl;
    efficiencies.push_back( vertexRecoEff );

    double lumi = luminosity();
    if(DEBUGMODE) cout << "Lumi: " << lumi << " pb^-1" << endl;

    double lumiUncorrected = luminosity(false, false);
    if(DEBUGMODE) cout << "Uncorrected lumi: " << lumiUncorrected << " pb^-1" << endl;

    double eff = efficiency(efficiencies);
    double CS = crossSection(yield, eff, lumi);
    if(DEBUGMODE) cout << "Final cross section: " << CS << " pb" << endl;
    setNominalCrossSection(CS);

    if(DEBUGMODE) cout << "Final cross section (corrected globally): " << crossSection(yieldGlobal, eff, lumi) << " pb" << endl;
    
    return CS;
    
}

pair<double,double> CrossSectionMaker::SysError(){
    // function which sums up all systematic errors (in fraction of the nominal cross section) separately for low and top errors and returns them

    /*
    if(isRomanPots){
        if(DEBUGMODE) cout << "Not running sys study for analysis with Roman Pots, returning precalculated values" << endl;
        return make_pair(0.101, 0.118);
    }else{
        if(DEBUGMODE) cout << "Not running sys study for analysis with Roman Pots, returning precalculated values" << endl;
    return make_pair(0.099,0.116 );
    }
    */

    vector<TString> sysErrorNames;
    vector<pair<double,double>> errors;
    // systematic error of cuts

    pair<double,double> combinedTPCErr = SysErrorTPC();
    errors.push_back(combinedTPCErr);
    sysErrorNames.push_back("Combined TPC");
    if(DEBUGMODE) cout << "Combined TPC sys error: " << combinedTPCErr.first << ", " << combinedTPCErr.second << endl;

    // 2 trigger efficiency sys errors: topology and veto
    pair<double,double> TEErr = CalculateTriggerTopologySysError();
    errors.push_back( TEErr );
    sysErrorNames.push_back("#epsilon_{topology}");
    if(DEBUGMODE) cout << "Trigger topology sys error: " << TEErr.first << ", " << TEErr.second << endl;

    pair<double, double> VEErr = CalculateTriggerVetoSysError();
    errors.push_back( VEErr);
    sysErrorNames.push_back("#epsilon_{veto}");
    if(DEBUGMODE) cout << "Trigger veto sys error: " << VEErr.first << ", " << VEErr.second << endl;

    // RP efficiency sys error
    pair<double, double> RPErr = CalculateRPEffSysError();
    errors.push_back(RPErr);
    sysErrorNames.push_back("#epsilon_{RP}");
    if(DEBUGMODE) cout << "RP efficiency sys error: " << RPErr.first << ", " << RPErr.second << endl;

    // yield extraction sys error
    pair<double,double> PFErr = CalculateYieldExtractionSysError();
    errors.push_back(PFErr);
    sysErrorNames.push_back("Yield extraction");
    if(DEBUGMODE) cout << "Yield extraction sys error: " << PFErr.first << ", " << PFErr.second << endl;
    
    // lumi sys error (van de Meer scans)
    pair<double,double> LErr = CalculateLumiSysError();
    errors.push_back(LErr);
    sysErrorNames.push_back("Van de Meer scans");
    if(DEBUGMODE) cout << "Lumi sys error: " << LErr.first << ", " << LErr.second << endl;

    // BEMC efficiency sys error (has to be at the end)
    pair<double,double> bemcErr = CalculateBemcSysError();
    errors.push_back(bemcErr);
    sysErrorNames.push_back("#epsilon_{BEMC}");
    if(DEBUGMODE) cout << "BEMC efficiency sys error: " << bemcErr.first << ", " << bemcErr.second << endl;

    pair<double,double> pileUpErr = CalculatePileUpSysError();
    errors.push_back(pileUpErr);
    sysErrorNames.push_back("#epsilon_{pile-up}");
    if(DEBUGMODE) cout << "Pile-up sys error: " << pileUpErr.first << ", " << pileUpErr.second << endl;


    TH1* hLow = new TH1D("sysErrorLowHist", "Sys Error Low", errors.size(), 0, errors.size());
    TH1* hTop = new TH1D("sysErrorTopHist", "Sys Error Top", errors.size(), 0, errors.size());

    // sum in percentages separately low and top
    double low = 0;
    double top = 0;
    int i = 1;
    for(auto err : errors){
        if(err.first == -1 || err.second == -1){
            cerr << "ERROR in CrossSectionMaker::SysError(): One of the systematic errors ("<< sysErrorNames[i-1] <<") is -1. Returning." << endl;
            return make_pair(-1,-1);
        }

        if(DEBUGMODE) cout << "Error " << sysErrorNames[i-1] << ": low " << err.first << ", top " << err.second << endl;
        low += pow(err.first,2);
        top += pow(err.second,2);
        hLow->SetBinContent(i, 1 - abs(err.first));
        hTop->SetBinContent(i, 1 + abs(err.second));
        i++;
    }
    double lowNoRp = sqrt(low - pow(RPErr.first,2));
    double topNoRp = sqrt(top - pow(RPErr.second,2));

    low = sqrt(low);
    top = sqrt(top);

    cout << "Total systematic error without RP efficiency: -" << lowNoRp*100 << "%, +" << topNoRp*100 << "%" << endl;
    cout << "Total systematic error with RP efficiency: -" << low*100 << "%, +" << top*100 << "%" << endl;

    createSysErrorPlot(hLow, hTop, sysErrorNames, "SysErrorOverview");

    // return pair
    return make_pair(low,top);
}

pair<double,double> CrossSectionMaker::SysErrorTPC(){
    // function which sums sys errors related to TPC 

    vector<TString> sysErrorNames;
    vector<pair<double,double>> errors;
    // systematic error of cuts
    pair<double,double> nHitsFitErr = CalculateSysErrorCuts(NHITSFIT);
    errors.push_back(nHitsFitErr);
    sysErrorNames.push_back("N^{fit}_{hits}");
    if(DEBUGMODE) cout << "nHitsFitErr: " << +nHitsFitErr.first << ", " << nHitsFitErr.second << endl;

    pair<double,double> nHitsDedxErr = CalculateSysErrorCuts(NHITSDEDX);
    errors.push_back(nHitsDedxErr);
    sysErrorNames.push_back("N^{dE/dx}_{hits}");
    if(DEBUGMODE) cout << "nHitsDedxErr: " << +nHitsDedxErr.first << ", " << nHitsDedxErr.second << endl;

    pair<double,double> etaErr = CalculateSysErrorCuts(ETA);
    errors.push_back(etaErr);
    sysErrorNames.push_back("#eta");
    if(DEBUGMODE) cout << "etaErr: " << +etaErr.first << ", " << etaErr.second << endl;

    pair<double,double> pidErr = CalculateSysErrorCuts(PID);
    errors.push_back(pidErr);
    sysErrorNames.push_back("#chi^{2}_{ee}");
    if(DEBUGMODE) cout << "pidErr: " << +pidErr.first << ", " << pidErr.second << endl;

    pair<double,double> nskErr = CalculateSysErrorCuts(NSIGMAKAON);
    errors.push_back(nskErr);
    sysErrorNames.push_back("n#sigma_{K}");
    if(DEBUGMODE) cout << "nsigmakaonErr: " << nskErr.first << ", " << nskErr.second << endl;

    pair<double,double> nspErr = CalculateSysErrorCuts(NSIGMAPROTON);
    errors.push_back(nspErr);
    sysErrorNames.push_back("n#sigma_{p}");
    if(DEBUGMODE) cout << "nsigmaprotonErr: " << nspErr.first << ", " << nspErr.second << endl;

    pair<double,double> dcaZErr = CalculateSysErrorCuts(DCAZINCM);
    errors.push_back(dcaZErr);
    sysErrorNames.push_back("DCA_{Z}");
    if(DEBUGMODE) cout << "dcaZErr: " << dcaZErr.first << ", " << dcaZErr.second << endl;

    pair<double,double> dcaXYErr = CalculateSysErrorCuts(DCAXYINCM);
    errors.push_back(dcaXYErr);
    sysErrorNames.push_back("DCA_{XY}");
    if(DEBUGMODE) cout << "dcaXYErr: " << dcaXYErr.first << ", " << dcaXYErr.second << endl;

    pair<double,double> VzErr = CalculateSysErrorCuts(VERTEXZ);
    errors.push_back(VzErr);
    sysErrorNames.push_back("V_{Z}");
    if(DEBUGMODE) cout << "VzErr: " << VzErr.first << ", " << VzErr.second << endl;

    // vertex z cut sys error
    pair<double,double> VZErr = CalculateVzSysError();
    errors.push_back(VZErr);
    sysErrorNames.push_back("#epsilon_{V_{Z}}");
    if(DEBUGMODE) cout << "Vz cut sys error: " << VZErr.first << ", " << VZErr.second << endl;

    //vertex reconstruction efficiency sys error
    pair<double,double> VRecoErr = CalculateVertexRecoSysError();
    errors.push_back(VRecoErr);
    sysErrorNames.push_back("#epsilon_{vertex reco}");
    if(DEBUGMODE) cout << "Vertex reco sys error: " << VRecoErr.first << ", " << VRecoErr.second << endl;

    TH1* hLow = new TH1D("tpcSysErrorLowHist", "TPC Sys Error Low", errors.size(), 0, errors.size());
    TH1* hTop = new TH1D("tpcSysErrorTopHist", "TPC Sys Error Top", errors.size(), 0, errors.size());

    // sum in percentages separately low and top
    double low = 0;
    double top = 0;
    int i = 1;
    for(auto err : errors){
        if(err.first == -1 || err.second == -1){
            cerr << "ERROR in CrossSectionMaker::SysError(): One of the systematic errors ("<< sysErrorNames[i-1] <<") is -1. Returning." << endl;
            return make_pair(-1,-1);
        }
        cout << "Error " << sysErrorNames[i-1] << ": low " << err.first << ", top " << err.second << endl;
        low += pow(err.first,2);
        top += pow(err.second,2);
        hLow->SetBinContent(i, 1 - abs(err.first));
        hTop->SetBinContent(i, 1 + abs(err.second));
        i++;
    }

    low = sqrt(low);
    top = sqrt(top);

    createSysErrorPlot(hLow, hTop, sysErrorNames, "CombinedTPCSysError");

    return make_pair(low,top);
}


void CrossSectionMaker::createSysErrorPlot(TGraph *gLow, TGraph *gTop, vector<TString> sysErrorNames, TString outName){
    
    if(DEBUGMODE) cout << "Now creating sys error plot in CrossSectionMaker::createSysErrorPlot()..." << endl;
    TCanvas *canvas = nullptr;
    mAna->CreateCanvas(&canvas, outName, 1200, 800);
    
    
    canvas->Clear();
    canvas->SetName("SysErrorOverview");
    //canvas->cd();


    mAna->SetTGraphStyle(gLow, kRed, 20);
    mAna->SetTGraphStyle(gTop, kBlue, 21);
    gLow->GetXaxis()->SetTitle("");
    gLow->GetYaxis()->SetTitle("Fraction of nominal cross section");
    gLow->GetXaxis()->SetLimits(0, gLow->GetN()+1);
    gLow->GetYaxis()->SetRangeUser(0.9, 1.1);

    TGraph *gNominal = new TGraph();
    gNominal->SetName("Nominal");

    gNominal->SetPoint(gNominal->GetN(), 0, 1.0);
    gNominal->SetPoint(gNominal->GetN(), gLow->GetN()+1, 1.0);
    mAna->SetTGraphStyle(gNominal, kBlack, 1);
    
    // set x axis labels
    gLow->GetXaxis()->SetNdivisions(gLow->GetN(), false);
    for(unsigned int i = 0; i < sysErrorNames.size(); i++){
        gLow->GetXaxis()->SetBinLabel(i+1, sysErrorNames[i].Data());
    }



    gLow->Draw("AP");
    gTop->Draw("same P");
    gNominal->Draw("same L");

    //draw legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.85);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(gLow, "Systematic error low", "p");
    legend->AddEntry(gTop, "Systematic error high", "p");
    legend->AddEntry(gNominal, "Nominal value", "l");

    legend->Draw("same");

    mAna->DrawSTARInternal();

    mOutFile->cd();
    canvas->Write(outName);
    canvas->Close();

    return;
}

void CrossSectionMaker::createSysErrorPlot(TH1 *&gLow, TH1 *&gTop, vector<TString> sysErrorNames, TString outName){
    
    if(DEBUGMODE) cout << "Now creating sys error plot in CrossSectionMaker::createSysErrorPlot()..." << endl;

    TCanvas *canvas = nullptr;
    mAna->CreateCanvas(&canvas, outName, 1200, 800);
    
    if(gLow->GetNbinsX() != gTop->GetNbinsX()){
        cerr << "ERROR in CrossSectionMaker::createSysErrorPlot(): gLow and gTop have different number of bins. Returning." << endl;
        return;
    }
    
    canvas->Clear();
    mAna->SetGPad();

    mAna->SetHistStyle(gLow, kRed, 20);
    mAna->SetHistStyle(gTop, kBlue, 21);
    gLow->GetXaxis()->SetTitleSize(0.07);
    gLow->GetXaxis()->SetLabelSize(0.06);
    gLow->GetXaxis()->SetLabelOffset(0.3);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.11);



    gLow->GetXaxis()->SetTitle("");
    gLow->GetYaxis()->SetTitle("Systematic Uncertainty");
    gLow->GetYaxis()->SetRangeUser(0.9, 1.1);

    TGraph *gNominal = new TGraph();
    gNominal->SetName("Nominal");

    gNominal->SetPoint(gNominal->GetN(), 0, 1.0);
    gNominal->SetPoint(gNominal->GetN(), gLow->GetNbinsX()+1, 1.0);
    mAna->SetTGraphStyle(gNominal, kBlack, 1);
    
    // set x axis labels
    for(unsigned int i = 0; i < sysErrorNames.size(); i++){
        gLow->GetXaxis()->SetBinLabel(i+1, sysErrorNames[i].Data());
    }

    if(DEBUGMODE) cout << "Set bin labels from sysErrorNames." << endl;

    gLow->Draw("P");
    gTop->Draw("same P");
    gNominal->Draw("same L");

    //draw legend
    TLegend *legend = new TLegend(0.6, 0.75, 0.42, 0.9);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(gLow, "Systematic uncertainty (low)", "p");
    legend->AddEntry(gTop, "Systematic uncertainty (top)", "p");
    legend->AddEntry(gNominal, "Nominal value", "l");

    legend->Draw("same");

    mAna->DrawSTARInternal();

    mOutFile->cd();
    canvas->Write(outName);
    canvas->Close();

}
pair<double,double> CrossSectionMaker::CalculateSysErrorCuts(int VAR){
    
    if(DEBUGMODE) cout << "Running systematic analysis of cut " << mUtil->nameOfVariable(VAR) << endl;

    // get yields from data
    double yieldTight = mAna->getYield(VAR, TIGHT);
    double yieldLoose = mAna->getYield(VAR, LOOSE);
    double yieldNominal = mAna->getYield(VAR,NOMINAL);
    
    cout << "Got the yields from analysis, tight: " << yieldTight << ", loose: " << yieldLoose << endl;
    // get reconstruction efficiency from embedding

    double recoEffLoose, recoEffTight, recoEffNominal;
    if(VAR == DCAZINCM || VAR == DCAXYINCM){ // dca is independent of track reconstruction efficiency
        recoEffTight = mEmbedding->getEfficiency(VAR, NOMINAL);
        recoEffLoose = mEmbedding->getEfficiency(VAR, NOMINAL);
        recoEffNominal = mEmbedding->getEfficiency(VAR, NOMINAL);
    }else{
        recoEffTight = mEmbedding->getEfficiency(VAR, TIGHT);
        recoEffLoose = mEmbedding->getEfficiency(VAR, LOOSE);
        recoEffNominal = mEmbedding->getEfficiency(VAR, NOMINAL);
    }
    cout << "reconstruction effs, tight: " << recoEffTight << ", loose: " << recoEffLoose << endl;

    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    cout << "Trigger eff: " << triggerEff << endl;
    vector<double> efficiencies;
    efficiencies.push_back( triggerEff );
    
    double RPEff = getRPEfficiency();
    cout << "RP eff: " << RPEff << endl;
    efficiencies.push_back( RPEff );

    double vertexRecoEff = getVertexRecoEfficiency();
    cout << "Vertex reco efficiency: " << vertexRecoEff << endl;
    efficiencies.push_back( vertexRecoEff );
    
    double lumi = luminosity();
    cout << "Lumi: " << lumi << endl;
    
    double vzEffTight, vzEffLoose, vzEff;
    vzEff = mAna->getVzStudyAverageEff();
    if(VAR == VERTEXZ){
        vzEffTight = mAna->getVzStudyEffTight();
        vzEffLoose = mAna->getVzStudyEffLoose();
    }else{
        vzEffTight = mAna->getVzStudyAverageEff();
        vzEffLoose = mAna->getVzStudyAverageEff();
    }
    cout << "Vertex z efficiency, tight: " << vzEffTight << ", loose: " << vzEffLoose << endl;
    
    double CSNominal = crossSection(yieldNominal, efficiency(efficiencies)*vzEff*recoEffNominal,lumi);
    cout << "Nominal cross section: " << CSNominal << endl;
    if(CSNominal == -1){
        return make_pair(-1,-1);
    }
    double effTight = efficiency(efficiencies)*vzEffTight*recoEffTight;
    double CSTight = crossSection(yieldTight, effTight, lumi);
    mSysErrorTight[VAR] = abs(CSTight - CSNominal)/CSNominal;
    cout << "Tight cross section: " << CSTight << ". Systematic error: " << mSysErrorTight[VAR] << endl;
    if(CSTight == -1){
        return make_pair(-1,-1);
    }

    double effLoose = efficiency(efficiencies)*vzEffLoose;
    double CSLoose = crossSection(yieldLoose, effLoose*recoEffLoose, lumi);
    mSysErrorLoose[VAR] = abs(CSLoose - CSNominal)/CSNominal;
    cout << "Loose cross section: " << CSLoose << ". Systematic error: " << mSysErrorLoose[VAR] << endl;
    if(CSLoose == -1){
        return make_pair(-1,-1);
    }
    
    return make_pair(mSysErrorTight[VAR], mSysErrorLoose[VAR]);
}


Value CrossSectionMaker::getYieldCorrection(int spectrumBins, double spectrumLow, double spectrumTop, TString condition, PlotEmbeddingJPsi *&embedding){
    // yield correction according to pT spectrum

    if(DEBUGMODE) cout << "Now calculating yield correction in CrossSectionMaker::getYieldCorrection() " << endl;

    // fit the peak and correct with embedding
    int nBins = 100;
    double lowerLim = 0.0;
    double upperLim = 5.0;
    if(isRomanPots){
        lowerLim = 2.5;
        upperLim = 3.5;
        nBins = 40;
    }
    TH1D* invMass = mAna->loadInvMassHist(nBins, lowerLim, upperLim, condition, true);
    if(!invMass || invMass->GetEntries() == 0){
        cerr << "ERROR: no invariant mass loaded for yield fitting and correction. " << endl;
        return Value();
    }

    if(DEBUGMODE) cout << "Loaded invMass with " << invMass->GetEntries() << " entries." << endl;

    FitJPsi *fit = new FitJPsi(invMass, "cb continuum");
    fit->setAlphaLoose(false);
    fit->setNLoose(false);
    fit->setFitRangeLow(0.4);
    fit->setFitRangeHigh(5.0);

    fit->fitPeak();
    fit->writeFitResult();
    fit->writeContinuumResult();

    if(DEBUGMODE) cout << "Fitted the peak." << endl;

    // get the pT spectrum from data
    condition += Form(" && invMass > %f && invMass < %f", fit->getLowLimitFit(), fit->getHighLimitFit());
    mAna->tree->Draw(Form("pt>>spectrum(%d,%f,%f)", spectrumBins, spectrumLow, spectrumTop), condition);


    TH1D *spectrum = (TH1D*)gPad->GetPrimitive("spectrum");
    if(!spectrum || spectrum->GetEntries() <= 0 ){
        cerr << "ERROR in CrossSectionMaker::getYieldCorrection(): Could not load spectrum for yield correction. The entries: " << spectrum->GetEntries() << endl;
        return Value();
    }
    
    mAna->bcgTree->Draw(Form("pt>>spectrumBcg(%d,%f,%f)", spectrumBins, spectrumLow, spectrumTop), condition);
    TH1D *spectrumBcg = (TH1D*)gPad->GetPrimitive("spectrumBcg");
    if( !spectrumBcg || spectrumBcg->GetEntries() <= 0){
        cerr << "ERROR in CrossSectionMaker::getYieldCorrection(): Could not load spectrum for yield correction. The entries: " << spectrumBcg->GetEntries() << endl;
        return Value();
    }
    spectrum->Add(spectrumBcg, -1.0);  // subtract background
    if(DEBUGMODE) cout << "Loaded pT spectrum with " << spectrum->GetEntries() << " entries." << endl;

    spectrum->SetName("pTSpectrumForCorrection");

    embedding->changeBinning(4, spectrumBins, spectrumLow, spectrumTop);  // 4 stands for pT of JPsi
    
    TGraphAsymmErrors* graph = embedding->reconstructionEfficiency(4, "extrapolatedPtRecoEff", false,true);


    if(!graph || graph == nullptr || graph->GetN() == 0){
        cerr << "ERROR in CrossSectionMaker::getYieldCorrection(): Could not load recostruction efficiency from embedding" << endl;
        return Value();
    }
    if(DEBUGMODE) cout << "Loaded reconstruction efficiency from embedding with " << graph->GetN() << " points." << endl;
    
    graph->SetName("recoEff_pTJPsi");

    mRawYield = fit->getYield();
    

    // fraction of signal in yield in range under the peak
    double frac = fit->getSOverB() <= 0 ? 1 : fit->getYield()/(fit->getYield() + fit->getYieldBackground());


    double correctedYield = 0;
    double correctedYieldErrTop = 0;
    double correctedYieldErrLow = 0;
    if(DEBUGMODE) cout << "Now running loop over spectrum for correction" << endl;
    for(int i = 1; i <= spectrum->GetNbinsX() ; ++i){ // loop over bins in the provided pT spectrum
      double binCenter = spectrum->GetBinCenter(i);
      for(int j = 0; j < graph->GetN(); j++){  // loop over points in the efficiency graph
         if(binCenter > (graph->GetX()[j]- graph->GetEXlow()[j]) && binCenter < (graph->GetX()[j]+ graph->GetEXhigh()[j])){
            //if(DEBUGMODE) cout << "Bin " << i << ", spectrum bin center: " << binCenter << ", graph bin center: " << graph->GetX()[j] << ", spectrum content = " << spectrum->GetBinContent(i) << ", graph y = " << graph->GetY()[j] << ", frac = " << frac << endl;
            correctedYield += spectrum->GetBinContent(i)*frac/graph->GetY()[j];
            correctedYieldErrTop += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYhigh()[j]/graph->GetY()[j], 2);
            correctedYieldErrLow += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYlow()[j]/graph->GetY()[j], 2);
            break;
         }
      }
    }

    //if(DEBUGMODE) cout << "Overflow bin: spectrum content = " << spectrum->GetBinContent(spectrum->GetNbinsX() + 1) << ", reconstruction eff for the last existing bin = " << graph->GetY()[graph->GetN()-1] << ", frac = " << frac << endl;
    correctedYield += spectrum->GetBinContent(spectrum->GetNbinsX() + 1)*frac/graph->GetY()[graph->GetN()-1];
    correctedYieldErrTop += pow(spectrum->GetBinError(spectrum->GetNbinsX() + 1)/spectrum->GetBinContent(spectrum->GetNbinsX() + 1), 2) + pow(graph->GetEYhigh()[graph->GetN()-1]/graph->GetY()[graph->GetN()-1], 2);
    correctedYieldErrLow += pow(spectrum->GetBinError(spectrum->GetNbinsX() + 1)/spectrum->GetBinContent(spectrum->GetNbinsX() + 1), 2) + pow(graph->GetEYlow()[graph->GetN()-1]/graph->GetY()[graph->GetN()-1], 2);


    correctedYield = round(correctedYield);
    correctedYieldErrTop = sqrt(correctedYieldErrTop);
    correctedYieldErrLow = sqrt(correctedYieldErrLow);


    if(DEBUGMODE) cout << "The raw yield from fitting: " << fit->getYield() << endl;
    if(DEBUGMODE) cout << "The corrected yield has been computed. Result is " << correctedYield << endl;
    if(DEBUGMODE) cout << "The corrected yield globally: " << fit->getYield() << "/" << embedding->getEfficiencyFinal() << " = " << fit->getYield()/embedding->getEfficiencyFinal() << endl;


    TPaveText *text = new TPaveText(0.7,0.52,0.85,0.49, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
    text->SetTextSize(0.03);
    text->SetFillColor(0);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    text->AddText(Form("Corrected yield = %.0f^{+%.0f}_{-%.0f}", correctedYield, correctedYieldErrTop, correctedYieldErrLow));
    text->Draw("same");

    TCanvas *fitCanvas = fit->getCanvas();
    if(!fitCanvas){
        cerr << "ERROR in CrossSectionMaker::getYieldCorrection(): Could not get fit canvas. Returning." << endl;
        return Value();
    }

    fitCanvas->SetName("CorrectedJPsi_canvas");
    fitCanvas->SetTitle("");
    fitCanvas->cd();
    mAna->DrawSTARpp510JPsi();

    if(isRomanPots){
        mAna->CreateText("With Roman Pots", 62,0.7,0.75,0.85,0.77);
    }else{
        mAna->CreateText("No Roman Pots", 62,0.7,0.75,0.85,0.77);
    }

    // save yield correction objects for later
    setYieldCorrectionCanvas(fitCanvas);
    setYieldCorrectionSpectrum(spectrum);
    setYieldCorrectionGraph(graph);

    Value result;
    result.val = correctedYield;
    result.errTop = correctedYieldErrTop;
    result.errLow = correctedYieldErrLow;

    return result;
}

void CrossSectionMaker::saveYieldCorrection(TString dir){
    if(!yieldCorrectionCanvas){
        cerr << "ERROR in CrossSectionMaker::saveYieldCorrection(): No yield correction canvas to save. Returning." << endl;
        return;
    }
    if(!hYieldCorrection){
        cerr << "ERROR in CrossSectionMaker::saveYieldCorrection(): No yield correction spectrum to save. Returning." << endl;
        return;
    }
    if(!gYieldCorrection){
        cerr << "ERROR in CrossSectionMaker::saveYieldCorrection(): No yield correction graph to save. Returning." << endl;
        return;
    }

    mOutFile->mkdir(dir);
    mOutFile->cd(dir);
    yieldCorrectionCanvas->SetName("YieldCorrectionCanvas");
    hYieldCorrection->SetName("YieldCorrectionSpectrum");
    gYieldCorrection->SetName("YieldCorrectionGraph");
    yieldCorrectionCanvas->Write();
    hYieldCorrection->Write();
    gYieldCorrection->Write();
    gDirectory->GetMotherDir()->cd();

}



double CrossSectionMaker::crossSection(double Y,double eff, double lumi, double deltaRap){

    if(eff*lumi*deltaRap == 0 ){
        cerr << "ERROR in CrossSectionMaker::crossSection(): Could not compute cross section, eff or lumi is 0" << endl;
        return -1;
    }

    return( Y/(eff*lumi*deltaRap) );
}

double CrossSectionMaker::efficiency(vector<double> effs){

    double prod = 1.0;
    for(auto e : effs){
        prod *= e;
    }
    return prod;
}

double CrossSectionMaker::luminosity(bool correctVeto, bool correctPileUp){
    double vetoA = 1.0;
    double vetoB = 0.0;
    double pileUpA = 1.0;
    double pileUpB = 0.0;

    if(correctVeto){
        vetoA = mGoodRun->mProbRetainEvent->getA();
        vetoB = mGoodRun->mProbRetainEvent->getB();
    }

    if(correctPileUp){
        pileUpA = mGoodRun->mProbRetainEvent->getPileUpCorrectionA();
        pileUpB = mGoodRun->mProbRetainEvent->getPileUpCorrectionB();
    }
    double lumi;
    if(isRomanPots){
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), vetoA, vetoB, pileUpA, pileUpB ); // A = 1.0, B = 0.0
    }else{
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), vetoA, vetoB, pileUpA, pileUpB ); // A = 1.0, B = 0.0
    }

    return lumi;

}


pair<double,double> CrossSectionMaker::CalculateRPEffSysError(){

    if(DEBUGMODE) cout << "Calculating systematic error of RP efficiency..." << endl;

    // get yields from data
    double yield = mCorrectedYield;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);

    //double RPEff = getRPEfficiency();
    Value RPEffError = getRPEffSysError();

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );
    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    double effLow = eff* (RPEffError.val - RPEffError.errLow);
    double effTop = eff* (RPEffError.val + RPEffError.errTop);

    double CSLow = crossSection(yield, effLow, lumi);
    double CSTop = crossSection(yield, effTop, lumi);
    double CS = crossSection(yield, eff*RPEffError.val, lumi);
    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    
    return make_pair(low, top);

}

pair<double,double> CrossSectionMaker::CalculateVertexRecoSysError(){
    if(DEBUGMODE) cout << "Calculating systematic error of vertex reconstruction efficiency..." << endl;

    // get yields from data
    double yield = mCorrectedYield;
    if(DEBUGMODE)   cout << "Corrected yield: " << yield << endl;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    if(DEBUGMODE) cout << "Trigger efficiency: " << triggerEff << endl;

    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);
    if(DEBUGMODE) cout << "RP efficiency: " << RPEff << endl;

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back( vzEff );
    if(DEBUGMODE) cout << "Vertex reconstruction efficiency: " << vzEff << endl;

    double vertexRecoEff = getVertexRecoEfficiency();
    double vertexRecoEffLow = getVertexRecoSysErrorLow();
    double vertexRecoEffTop = getVertexRecoSysErrorTop();
    if(DEBUGMODE) cout << "Vertex reconstruction efficiency: " << vertexRecoEff << " +/- " << vertexRecoEffLow << " (low), " << vertexRecoEffTop << " (top)" << endl;

    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    double effTop = eff* (vertexRecoEff + vertexRecoEffTop);
    double effLow = eff* (vertexRecoEff - vertexRecoEffLow);

    double CSLow = crossSection(yield, effTop, lumi);
    double CSTop = crossSection(yield, effLow, lumi);
    double CS = crossSection(yield, eff*vertexRecoEff, lumi);

    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    
    return make_pair(low, top);

}


pair<double,double> CrossSectionMaker::CalculateVzSysError(){
    
    cout << "Calculating systematic error of vertex z efficiency method..." << endl;

    // get yields from data
    double yield = mCorrectedYield;

    vector<double> efficiencies;
    efficiencies.clear();
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    if(DEBUGMODE) cout << "Trigger efficiency: " << triggerEff << endl;

    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);
    if(DEBUGMODE) cout << "RP efficiency: " << RPEff << endl;

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );
    if(DEBUGMODE) cout << "Vertex reconstruction efficiency: " << vertexRecoEff << endl;

    double vzEff = mAna->getVzStudyAverageEff();
    double vzError = mAna->getVzStudySysError(); // divide by 2 to get one side error
    if(DEBUGMODE) cout << "Vertex z efficiency: " << vzEff << " +/- " << vzError << endl;

    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    double effLow = eff* (vzEff - vzError);
    double effTop = eff* (vzEff + vzError);
    if(DEBUGMODE) cout << "Combined efficiency: " << eff << ", effLow: " << effLow << ", effTop: " << effTop << endl;

    double CSLow = crossSection(yield, effLow, lumi);
    double CSTop = crossSection(yield, effTop, lumi);
    double CS = crossSection(yield, eff*vzEff, lumi);

    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    
    return make_pair(low, top);
}

pair<double,double> CrossSectionMaker::CalculateYieldExtractionSysError(){
    
    cout << "Calculating systematic error of yield extraction..." << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);

    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );

    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    
    double errorPF = mAna->getPeakFittingSigma();

    double CSLow = crossSection(mCorrectedYield -errorPF/2, eff, lumi);
    double CSTop = crossSection(mCorrectedYield + errorPF/2, eff, lumi);

    double CS = crossSection(mCorrectedYield, eff, lumi);
    //errors in percentages
    double low = (CSLow - CS)/CS;
    double top = (CSTop - CS)/CS;
    
    return make_pair(low, top);


}
pair<double,double> CrossSectionMaker::CalculateBemcSysError(){

    if(DEBUGMODE) cout << "Calculating BEMS systematic error in CrossSectionMaker::CalculateBemcSysError()" << endl;

    TString path = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/EmbeddingJPsi/extreme/";
   
    embMin = new PlotEmbeddingJPsi((path + TString("outMin.root")).Data(), mOutFile);
    embMax = new PlotEmbeddingJPsi((path + TString("outMax.root")).Data(), mOutFile);

    embMin->Init();
    embMax->Init();

    if(DEBUGMODE) cout << "Fitting peak for max efficiency" << endl;
    embMax->fitPeakJPsi();
    double efficiencyMax = embMax->getEfficiencyFinal();

    if(DEBUGMODE) cout << "Fitting peak for min efficiency" << endl;
    embMin->fitPeakJPsi();
    double efficiencyMin = embMin->getEfficiencyFinal();

    double efficiencyNominal = mEmbedding->getEfficiencyFinal();

    if(DEBUGMODE) cout << "Now calculating cross section" << endl;
    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back( triggerEff );

    double vzEff = getVertexRecoEfficiency();
    efficiencies.push_back( vzEff );

    double RPEff = getRPEfficiency();
    efficiencies.push_back( RPEff );

    double vertexRecoEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back( vertexRecoEff );
    
    double lumi = luminosity();

    double eff = efficiency(efficiencies);
    double CSNominal  = crossSection(mAna->getYieldFinal() , eff*efficiencyNominal, lumi);
    double CSMin = crossSection(mAna->getYieldFinal() , eff*efficiencyMax, lumi);
    double CSMax = crossSection(mAna->getYieldFinal() , eff*efficiencyMin, lumi);

    //errors in percentages
    double low = abs(CSNominal - CSMin)/CSNominal;
    double top = abs(CSMax - CSNominal)/CSNominal;

    if(DEBUGMODE) cout << "Finished calculating BEMC systematic error in CrossSectionMaker::CalculateBemcSysError()" << endl;
    return make_pair(low, top);
}

pair<double,double> CrossSectionMaker::CalculateLumiSysError(){
    // this will be from van de meer scans
    if(DEBUGMODE) cout << "Starting  CrossSectionMaker::CalculateLumiSysError()" << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );

    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    
    double lumiLow = lumi + 0.05*lumi;
    double lumiTop = lumi - 0.05*lumi;

    double CSLow = crossSection(mCorrectedYield, eff, lumiLow);
    double CSTop = crossSection(mCorrectedYield, eff, lumiTop);
    double CS = crossSection(mCorrectedYield, eff, lumi);

    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    

    if(DEBUGMODE) cout << "Finished with CrossSectionMaker::CalculateLumiSysError()" << endl;
    return make_pair(low, top);

    
}
pair<double,double> CrossSectionMaker::CalculateTriggerVetoSysError(){
    // this systematic error comes from varying fit parameters for probability of retaining an event in their errors
    
    if(DEBUGMODE) cout << "Calculating systematic error of veto conditions in trigger..." << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );

    double eff = efficiency(efficiencies);

    double lumi, lumiLow, lumiTop;
    if(isRomanPots){
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiLow = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() - mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiTop = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() + mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
    }else{
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiLow = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() - mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiTop = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() + mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
    }

    if(DEBUGMODE) cout << "Luminosity: " << lumi << ", lumiLow: " << lumiLow << ", lumiTop: " << lumiTop << endl;
    if(DEBUGMODE) cout << "B value for veto calculation: " << mGoodRun->mProbRetainEvent->getB() << " +/- " << mGoodRun->mProbRetainEvent->getBError() << endl;

    double CSLow = crossSection(mCorrectedYield, eff, lumiTop);
    double CSTop = crossSection(mCorrectedYield, eff, lumiLow);
    double CS = crossSection(mCorrectedYield, eff, lumi);

    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    
    return make_pair(low, top);

}

pair<double,double> CrossSectionMaker::CalculateTriggerTopologySysError(){
    // this comes from choosing different widths of windows of inv mass 

    if(DEBUGMODE) cout << "Calculating systematic error of trigger topology..." << endl;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    double triggerEffError = mZeroBias->getTriggerEfficiencySysError()/2;  // returns 1 sigma, we want half sigma
    if(DEBUGMODE) cout << "Trigger efficiency: " << triggerEff << " +/- " << triggerEffError << endl;

    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);
    if(DEBUGMODE) cout << "RP efficiency: " << RPEff << endl;

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);
    if(DEBUGMODE) cout << "Vertex z efficiency: " << vzEff << endl;

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );
    if(DEBUGMODE) cout << "Vertex reconstruction efficiency: " << vertexRecoEff << endl;

    double lumi = luminosity();
    if(DEBUGMODE) cout << "Luminosity: " << lumi << endl;
    
    double eff = efficiency(efficiencies);
    
    double effLow = eff*(triggerEff - triggerEffError);
    double effTop = eff*(triggerEff + triggerEffError);
    if(DEBUGMODE) cout << "Combined efficiency: " << eff << ", effLow: " << effLow << ", effTop: " << effTop << endl;

    double CSLow = crossSection(mCorrectedYield, effTop, lumi);
    double CSTop = crossSection(mCorrectedYield, effLow, lumi);
    double CS = crossSection(mCorrectedYield, eff*triggerEff, lumi);
    if(DEBUGMODE) cout << "Cross sections: CS = " << CS << ", CSLow = " << CSLow << ", CSTop = " << CSTop << endl;
    //errors in percentages
    double low = (CSLow - CS)/CS;
    double top = (CSTop - CS)/CS;
    if(DEBUGMODE) cout << "Systematic errors: low = " << low << ", top = " << top << endl;

    
    return make_pair(low, top);

}


pair<double,double> CrossSectionMaker::CalculatePileUpSysError(){

    if(DEBUGMODE) cout << "Calculating systematic error for pile-up correction..." << endl;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vzEff = mAna->getVzStudyAverageEff();
    efficiencies.push_back(vzEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back( vertexRecoEff );

    double eff = efficiency(efficiencies);

    double lumi, lumiLow, lumiTop;
    if(isRomanPots){
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiLow = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() - mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA() - mGoodRun->mProbRetainEvent->getPileUpCorrectionAError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() - mGoodRun->mProbRetainEvent->getPileUpCorrectionBError() ); // A = 1.0, B = 0.0
        lumiTop = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() + mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA() + mGoodRun->mProbRetainEvent->getPileUpCorrectionAError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() + mGoodRun->mProbRetainEvent->getPileUpCorrectionBError()); // A = 1.0, B = 0.0
    }else{
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() ); // A = 1.0, B = 0.0
        lumiLow = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() - mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA() - mGoodRun->mProbRetainEvent->getPileUpCorrectionAError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() - mGoodRun->mProbRetainEvent->getPileUpCorrectionBError() ); // A = 1.0, B = 0.0
        lumiTop = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), 1, mGoodRun->mProbRetainEvent->getB() + mGoodRun->mProbRetainEvent->getBError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionA() + mGoodRun->mProbRetainEvent->getPileUpCorrectionAError(), mGoodRun->mProbRetainEvent->getPileUpCorrectionB() + mGoodRun->mProbRetainEvent->getPileUpCorrectionBError()); // A = 1.0, B = 0.0
    }
    if(DEBUGMODE) cout << "luminosity values: " << lumi << ", lumiLow: " << lumiLow << ", lumiTop: " << lumiTop << endl;
    if(DEBUGMODE) cout << "pile up correction A: " << mGoodRun->mProbRetainEvent->getPileUpCorrectionA() << " +/- " << mGoodRun->mProbRetainEvent->getPileUpCorrectionAError() << endl;
    if(DEBUGMODE) cout << "pile up correction B: " << mGoodRun->mProbRetainEvent->getPileUpCorrectionB() << " +/- " << mGoodRun->mProbRetainEvent->getPileUpCorrectionBError() << endl;

    if(DEBUGMODE) cout << "Luminosity: " << lumi << ", lumiLow: " << lumiLow << ", lumiTop: " << lumiTop << endl;

    double CSLow = crossSection(mCorrectedYield, eff, lumiTop);
    double CSTop = crossSection(mCorrectedYield, eff, lumiLow);
    double CS = crossSection(mCorrectedYield, eff, lumi);

    //errors in percentages
    double low = abs(CSLow - CS)/CS;
    double top = abs(CSTop - CS)/CS;
    
    return make_pair(low, top);
}


void CrossSectionMaker::rapidityDependencePlot(int nBins, double low, double top,vector<double> plotEdges){
    
    mOutFile->mkdir("RapidityDependence");
    mOutFile->cd("RapidityDependence");


    TGraphAsymmErrors *gRap = new TGraphAsymmErrors();
    gRap->SetName("rapidityGraph");


    for(unsigned int i = 0; i < plotEdges.size()-1 ; i++){
        TString condition = mAna->getCondition() + Form(" && abs(pairRapidity) > %.3f && abs(pairRapidity) < %.3f", plotEdges[i], plotEdges[i+1]);
        

        Value res = getYieldCorrection(30, 0, 1.5, condition , mEmbedding);
        if(DEBUGMODE) cout << "Rapidity bin: " << i << ", Yield Correction: " << res.val << ", condition: " << condition << endl;
        if(DEBUGMODE) cout << "Yield err top: " << res.errTop << ", yield err low: " << res.errLow << endl;
        
        mAna->DrawSTARInternal();
        mAna->CreateText(Form("%.1f < y < %.1f", plotEdges[i], plotEdges[i+1]), 62,0.7,0.75,0.85,0.77);
        saveYieldCorrection(Form("YieldCorrectionRap(%.1f,%.1f)", plotEdges[i], plotEdges[i+1]));

        double cs = crossSection(res.val, getVertexRecoEfficiency()*getRPEfficiency()*mAna->getVzStudyAverageEff()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]);

        double errTop = crossSection(res.val + res.errTop, getVertexRecoEfficiency()*getRPEfficiency()*mAna->getVzStudyAverageEff()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]) - cs;
        double errLow = cs - crossSection(res.val - res.errLow, getVertexRecoEfficiency()*getRPEfficiency()*mAna->getVzStudyAverageEff()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]);
        gRap->SetPoint(i, 0.5*(plotEdges[i]+plotEdges[i+1]), cs);
        gRap->SetPointError(i, 0.5*(plotEdges[i+1]-plotEdges[i]), 0.5*(plotEdges[i+1]-plotEdges[i]), errLow, errTop);

    }
    
    TCanvas *canvas = nullptr;
    mAna->CreateCanvas(&canvas, "RapidityDependence", 800, 600);

    TGraph* gMc = getMcNormalizedPlot(50, -1, 1);
    canvas->Clear();
    
    mAna->SetTGraphStyle(gRap, kBlue, 20);
    gRap->SetMarkerStyle(20);
    gRap->SetMarkerSize(2);
    gRap->SetLineColor(kBlue);
    gRap->SetLineWidth(2);
    gRap->SetLineColor(kBlue);
    gRap->GetYaxis()->SetRangeUser(0, 800);
    gRap->GetXaxis()->SetTitle("y [-]");
    gRap->GetYaxis()->SetTitle("#frac{d#sigma}{dy} [pb]");
    gRap->Draw("AEP");

    gMc->Draw("L");

    mAna->DrawSTARpp510JPsi();

       
    TLegend *l;
    mAna->CreateLegend(&l, 0.2,0.8,0.4,0.9);
    l->AddEntry(gRap, "Data", "EP");
    l->AddEntry(gMc, "STARlight MC", "LP");
    l->Draw("same");
    
    mOutFile->cd("RapidityDependence");
    gRap->Write();
    canvas->Write();
    mOutFile->cd();

}


TGraph* CrossSectionMaker::getMcNormalizedPlot(int nBins, double low, double top){

    cout << "Getting normalized MC rapidity distribution..." << endl;

    mEmbedding->starlightTree->Draw(Form("rapVM>>histMc(%d,%f,%f)", nBins, low, top));
    TH1D* histMc = (TH1D*)gPad->FindObject("histMc");
    if(!histMc){
        cerr << "ERROR in CrossSectionMaker::getMcNormalizedPlot(): Could not find histogram histMc. Returning." << endl;
        return nullptr;
    }

    double crossSectionMC = 49.815; // pb
    int nEntries = mEmbedding->starlightTree->GetEntries();

    TGraph *gMc = new TGraph();
    gMc->SetName("mcRapidityDistribution");
    gMc->SetTitle(" ; y [-]; d#sigma/dy [pb]");


    for(int i = 1; i <= histMc->GetNbinsX(); i++){
        double binWidth = histMc->GetBinWidth(i);
        double content = histMc->GetBinContent(i);
        //double err = histMc->GetBinError(i);
        double binCenter = histMc->GetBinCenter(i);

        if(content == 0) continue;

        double newContent = content * crossSectionMC / (nEntries * binWidth);

        gMc->SetPoint(i-1, binCenter, newContent);
        if(DEBUGMODE) cout << "Bin " << i << ": center = " << binCenter << ", content = " << content << ", new content = " << newContent << endl;

    }

    mEmbedding->SetTGraphStyle(gMc, kRed, 21);
    gMc->GetXaxis()->SetTitle("y [-]");
    gMc->GetYaxis()->SetTitle("d#sigma/dy [pb]");

    return gMc; 

}



