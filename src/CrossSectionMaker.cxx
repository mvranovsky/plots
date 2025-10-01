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
        cout << "Could not open output file" << endl;
    }

    // initialize analysis class for data
    cout << "Initiating analysis class..." << endl;
    mAna = new PlotAnaJPsi(AnaDir.Data(), mOutFile);
    mAna->Init();

    isRomanPots = mAna->isRomanPots();
    cout << "This analysis is with RPs: " << isRomanPots << endl;


    // initialize embedding class
    cout << "Initializing embedding class..." << endl;
    mEmbedding = new PlotEmbeddingJPsi(EmbeddingDir.Data(), mOutFile);
    mEmbedding->Init();

    // initialize goodrun class
    cout << "Initializing good run class..." << endl;
    mGoodRun = new PlotGoodRun(GoodRunDir.Data(), mOutFile);
    mGoodRun->Init();

    cout << "Initializing zero bias class..." << endl;
    mZeroBias = new PlotZeroBias(ZeroBiasDir.Data(), mOutFile);
    mZeroBias->Init();

    // initialize ProbRetainEvent class
    cout << "Initializing ProbRetainEvent class..." << endl;
    TTree *mGoodRunTreeBcg = nullptr;
    mGoodRun->ConnectInputTree(GoodRunDir.Data(), nameOfAnaGoodRunTree.Data(), mGoodRunTree, mGoodRunTreeBcg);

    if(!mGoodRunTree || mGoodRunTree->GetEntries() == 0){
        cerr << "Could not open good run tree. Returning." << endl;
        return;
    }

    //mProbRetainEvent = new ProbRetainEvent(mOutFile, mGoodRunTree);

}

CrossSectionMaker::~CrossSectionMaker(){
    mAna->Finish();
    mEmbedding->Finish();
    mGoodRun->Finish();
    mZeroBias->Finish();

    //if(mOutFile) mOutFile->Close();
    //if(mUtil) delete mUtil;
    //if(mProbRetainEvent) delete mProbRetainEvent;
    if(mEmbedding) delete mEmbedding;
    if(mAna) delete mAna;
    if(mGoodRun) delete mGoodRun;
}

void CrossSectionMaker::Make(){

    mAna->Make();

    mEmbedding->Make();

    mGoodRun->Make();

    mZeroBias->Make();

    //mProbRetainEvent->Make();
    cout << "-------------------------------------------------------------------------------------------------" << endl;
    cout << "-------------------------------------------------------------------------------------------------" << endl;
    cout << "-------------------------------------------------------------------------------------------------" << endl;
    cout << "Running analysis with Roman Pots: " << isRomanPots << endl;

    cout << "Yield correction: " << getYieldCorrection(30, 0, 1.5, mAna->getCondition() ) << endl;
    saveYieldCorrection("YieldCorrectionMain");
    setYieldCorrFinal( mCorrectedYield );
    setYieldCorrErrTopFinal( mCorrYieldErrTop );
    setYieldCorrErrLowFinal( mCorrYieldErrLow );


    // calculate cross section
    CalculateCrossSection();

    /*
    pair<double,double> systematics = SysError();
    
    cout << "Final result: " << crossSection << " pb" << endl;
    cout << "Error low: " << systematics.first*crossSection << " pb" << endl;
    cout << "Error high: " << systematics.second*crossSection << " pb" << endl;
    */

    if(!isRomanPots){
        cout << "Running rapidity dependence plot..." << endl;
        rapidityDependencePlot(40, 2.5,3.5, {-1, -0.3, 0.3, 1.0} );
    }

    cout << "Finished CrossSectionMaker::Make()" << endl;

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
    cout << "Trigger eff: " << triggerEff << endl;
    efficiencies.push_back( triggerEff );

    
    double RPEff = getRPEfficiency();
    cout << "RP eff: " << RPEff << endl;
    efficiencies.push_back( RPEff );

    double vertexRecoEff = getVertexRecoEfficiency();
    cout << "Vertex reco efficiency: " << vertexRecoEff << endl;
    efficiencies.push_back( vertexRecoEff );
    
    double lumi = luminosity();
    cout << "Lumi: " << lumi << " pb^-1" << endl;

    double eff = efficiency(efficiencies);
    double CS = crossSection(yield, eff, lumi);
    cout << "Final cross section: " << CS << " pb" << endl;
    setNominalCrossSection(CS);

    cout << "Final cross section (corrected globally): " << crossSection(yieldGlobal, eff, lumi) << " pb" << endl;
    
    return CS;
    
}

pair<double,double> CrossSectionMaker::SysError(){
    // function which sums up all systematic errors (in percentages of the nominal yield) separately for low and top errors and returns them

    vector<pair<double,double>> errors;
    // systematic error of cuts
    pair<double,double> nHitsFitErr = CalculateSysErrorCuts(NHITSFIT);
    errors.push_back(nHitsFitErr);

    pair<double,double> nHitsDedxErr = CalculateSysErrorCuts(NHITSDEDX);
    errors.push_back(nHitsDedxErr);
    
    pair<double,double> etaErr = CalculateSysErrorCuts(ETA);
    errors.push_back(etaErr);
    
    pair<double,double> pidErr = CalculateSysErrorCuts(PID);
    errors.push_back(pidErr);


    // how to calculate sys errors of cuts dcaZ, dcaXY, vertexZ?

    // BEMC efficiency sys error


    // 2 trigger efficiency sys errors: topology and veto
    pair<double,double> TEErr = CalculateTriggerTopologySysError();
    errors.push_back( TEErr );

    pair<double, double> VEErr = CalculateTriggerVetoSysError();
    errors.push_back( VEErr);

    // RP efficiency sys error
    pair<double, double> RPErr = CalculateRPEffSysError();
    errors.push_back(RPErr);

    // yield extraction sys error
    pair<double,double> PFErr = CalculateYieldExtractionSysError();
    errors.push_back(PFErr);

    // vertex reconstruction sys error
    pair<double,double> VRErr = CalculateVertexRecoSysError();
    errors.push_back(VRErr);

    // lumi sys error (van de Meer scans)
    pair<double,double> LErr = CalculateLumiSysError();
    errors.push_back(LErr);


    // sum in percentages separately low and top
    double low = 0;
    double top = 0;
    for(auto err : errors){
        low += pow(err.first,2);
        top += pow(err.second,2);
    }
    low = sqrt(low);
    top = sqrt(top);

    // return pair
    return make_pair(low,top);
}




pair<double,double> CrossSectionMaker::CalculateSysErrorCuts(int VAR){
    
    cout << "Running systematic analysis of cut " << mUtil->nameOfVariable(VAR) << endl;

    // get yields from data
    double yieldTight = mAna->getYield(VAR, TIGHT);
    double yieldLoose = mAna->getYield(VAR, LOOSE);
    
    cout << "Got the yields from analysis" << endl;
    // get reconstruction efficiency from embedding
    double recoEffTight = mEmbedding->getEfficiency(VAR, TIGHT);
    double recoEffLoose = mEmbedding->getEfficiency(VAR, LOOSE);
    cout << "Got the yields for sys study " << mUtil->nameOfVariable(VAR) << endl;

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

    double CSNominal = getNominalCrossSection();
    cout << "Nominal cross section: " << CSNominal << endl;

    double effTight = efficiency(efficiencies);
    double CSTight = crossSection(yieldTight, effTight*recoEffTight, lumi);
    mSysErrorTight[VAR] = abs(CSTight - CSNominal)/CSNominal;
    cout << "Tight cross section: " << CSTight << ". Systematic error: " << mSysErrorTight[VAR] << endl;


    double effLoose = efficiency(efficiencies);
    double CSLoose = crossSection(yieldLoose, effLoose*recoEffLoose, lumi);
    mSysErrorLoose[VAR] = abs(CSLoose - CSNominal)/CSNominal;
    cout << "Loose cross section: " << CSLoose << ". Systematic error: " << mSysErrorLoose[VAR] << endl;

    return make_pair(mSysErrorTight[VAR], mSysErrorLoose[VAR]);
}




double CrossSectionMaker::getYieldCorrection(int spectrumBins, double spectrumLow, double spectrumTop, TString condition){
    // yield correction according to pT spectrum


    // fit the peak and correct with embedding
    int nBins = 40;
    double lowerLim = 2.5;
    double upperLim = 3.5;
    TH1D* invMass = mAna->loadInvMassHist(nBins, lowerLim, upperLim, condition);
    if(!invMass || invMass->GetEntries() == 0){
        cerr << "ERROR: no invariant mass loaded for yield fitting and correction. " << endl;
        return -1;
    }

    FitJPsi *fit = new FitJPsi(invMass, "Poly1");
    fit->fitPeak();


    // get the pT spectrum from data
    condition += Form(" && invMass > %f && invMass < %f", fit->getLowLimitFit(), fit->getHighLimitFit());
    mAna->tree->Draw(Form("pt>>spectrum(%d,%f,%f)", spectrumBins, spectrumLow, spectrumTop), condition);

    TH1D *spectrum = (TH1D*)gPad->GetPrimitive("spectrum");
    if(!spectrum || spectrum->GetEntries() <= 0){
        cerr << "Error: Could not load spectrum for yield correction. " << endl;
        return -1;
    }

    mAna->bcgTree->Draw(Form("pt>>spectrumBcg(%d,%f,%f)", spectrumBins, spectrumLow, spectrumTop),condition);
    TH1D *spectrumBcg = (TH1D*)gPad->GetPrimitive("spectrumBcg");
    if(!spectrumBcg || spectrumBcg->GetEntries() <= 0){
        cerr << "Error: Could not load spectrum background for yield correction. " << endl;
        return -1;
    }
    spectrum->SetName("pTSpectrumForCorrectionWithoutBcg");
    spectrumBcg->SetName("pTSpectrumBcg");
    mOutFile->cd();
    spectrum->Write();
    spectrumBcg->Write();

    spectrum->Add(spectrumBcg, -1);

    spectrum->SetName("pTSpectrumForCorrection");


    // get the efficiency w.r.t. pT from embedding
    mEmbedding->changeBinning(4, spectrumBins, spectrumLow, spectrumTop);  // 4 stands for pT of JPsi
    TGraphAsymmErrors* recoEff = mEmbedding->reconstructionEfficiency(4, "extrapolatedPtRecoEff", true);
    
    if(!recoEff || recoEff->GetN() == 0){
        cerr << "ERROR: Could not load recostruction efficiency from embedding" << endl;
        return -1;
    }

    recoEff->SetName("recoEff_pTJPsi");
    mRawYield = fit->getYield();
    
    mCorrectedYield = fit->getCorrectedYield(spectrum, recoEff, mEmbedding->getEfficiencyFinal());

    mCorrYieldErrTop = fit->getCorrectedYieldErrTop();
    mCorrYieldErrLow = fit->getCorrectedYieldErrLow();

    TCanvas *fitCanvas = fit->getCanvas();
    if(!fitCanvas){
        cout << "Could not get fit canvas. Returning." << endl;
        return -1;
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
    setYieldCorrectionGraph(recoEff);

    return mCorrectedYield;
}

void CrossSectionMaker::saveYieldCorrection(TString dir){
    if(!yieldCorrectionCanvas){
        cerr << "No yield correction canvas to save. Returning." << endl;
        return;
    }
    if(!hYieldCorrection){
        cerr << "No yield correction spectrum to save. Returning." << endl;
        return;
    }
    if(!gYieldCorrection){
        cerr << "No yield correction graph to save. Returning." << endl;
        return;
    }

    mOutFile->cd();
    mOutFile->mkdir(dir);
    mOutFile->cd(dir);
    yieldCorrectionCanvas->SetName("YieldCorrectionCanvas");
    hYieldCorrection->SetName("YieldCorrectionSpectrum");
    gYieldCorrection->SetName("YieldCorrectionGraph");
    yieldCorrectionCanvas->Write();
    hYieldCorrection->Write();
    gYieldCorrection->Write();
    mOutFile->cd();

}



double CrossSectionMaker::crossSection(double Y,double eff, double lumi, double deltaRap){

    return( Y/(eff*lumi*deltaRap) );
}

double CrossSectionMaker::efficiency(vector<double> effs){

    double prod = 1.0;
    for(auto e : effs){
        prod *= e;
    }
    return prod;
}

double CrossSectionMaker::luminosity(){

    double lumi;
    if(isRomanPots){
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), mGoodRun->mProbRetainEvent->getA(), mGoodRun->mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }else{
        lumi = mGoodRun->mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), mGoodRun->mProbRetainEvent->getA(), mGoodRun->mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }

    return lumi;

}

pair<double,double> CrossSectionMaker::CalculateRPEffSysError(){

    cout << "Calculating systematic error of RP efficiency..." << endl;

    // get yields from data
    double yield = mCorrectedYield;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);

    
    double RPEff = getRPEfficiency();
    double RPEffError = getRPEffSysError();

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back(vertexRecoEff);

    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    double effLow = eff* (RPEff + RPEffError);
    double effTop = eff* (RPEff - RPEffError);

    double CSLow = crossSection(yield, effLow, lumi);
    double CSTop = crossSection(yield, effTop, lumi);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
    return make_pair(low, top);

}

pair<double,double> CrossSectionMaker::CalculateVertexRecoSysError(){
    
    cout << "Calculating systematic error of vertex reconstruction efficiency..." << endl;

    // get yields from data
    double yield = mCorrectedYield;

    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);

    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    double vertexRecoEffError = getVertexRecoEffSysError();

    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    double effLow = eff* (vertexRecoEff + vertexRecoEffError);
    double effTop = eff* (vertexRecoEff - vertexRecoEffError);

    double CSLow = crossSection(yield, effLow, lumi);
    double CSTop = crossSection(yield, effTop, lumi);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
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

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back(vertexRecoEff);
    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    
    double yieldPF = mAna->getPeakFittingMean();
    double errorPF = mAna->getPeakFittingSigma();

    if(abs(mRawYield- yieldPF) > errorPF){
        cerr << "The raw yield from correction is quite different from the peak fitting study. Try different binning for the correction..." << endl;
        cerr << "Raw yield: " << mRawYield << endl;
        cerr << "Yield from fitting study: " << yieldPF << endl;
        return make_pair(-999,-999);
    }

    double CSLow = crossSection(yieldPF -errorPF, eff, lumi);
    double CSTop = crossSection(yieldPF + errorPF, eff, lumi);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
    return make_pair(low, top);


}
pair<double,double> CrossSectionMaker::CalculateBemcSysError(){
    
    return make_pair(0,0);
}
pair<double,double> CrossSectionMaker::CalculateLumiSysError(){
    // this will be from van de meer scans
    cout << "Calculating systematic error of veto conditions in trigger..." << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back(vertexRecoEff);
    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    
    double lumiLow = lumi + 0.05*lumi;
    double lumiTop = lumi - 0.05*lumi;

    double CSLow = crossSection(mCorrectedYield, eff, lumiLow);
    double CSTop = crossSection(mCorrectedYield, eff, lumiTop);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
    return make_pair(low, top);

    
}
pair<double,double> CrossSectionMaker::CalculateTriggerVetoSysError(){
    // this systematic error comes from varying fit parameters for probability of retaining an event in their errors
    
    cout << "Calculating systematic error of veto conditions in trigger..." << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    efficiencies.push_back(triggerEff);
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back(vertexRecoEff);
    
    double lumi = luminosity();
    double lumiErr = mGoodRun->mProbRetainEvent->calculateSysError(isRomanPots ? mGoodRun->getRpGoodRunList() : mGoodRun->getGoodRunList());
    
    double eff = efficiency(efficiencies);
    
    double lumiLow = lumi + lumiErr;
    double lumiTop = lumi - lumiErr;

    double CSLow = crossSection(mCorrectedYield, eff, lumiLow);
    double CSTop = crossSection(mCorrectedYield, eff, lumiTop);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
    return make_pair(low, top);

}

pair<double,double> CrossSectionMaker::CalculateTriggerTopologySysError(){
    // this comes from choosing different widths of windows of inv mass 

    cout << "Calculating systematic error of trigger topology..." << endl;


    vector<double> efficiencies;
    // get trigger efficiency from good run
    double triggerEff = mZeroBias->getTriggerEfficiency();
    double triggerEffError = mZeroBias->getTriggerEfficiencySysError();
    
    double RPEff = getRPEfficiency();
    efficiencies.push_back(RPEff);

    double vertexRecoEff = getVertexRecoEfficiency();
    efficiencies.push_back(vertexRecoEff);
    
    double lumi = luminosity();
    
    double eff = efficiency(efficiencies);
    
    double effLow = eff*(triggerEff + triggerEffError);
    double effTop = eff*(triggerEff - triggerEffError);


    double CSLow = crossSection(mCorrectedYield, effLow, lumi);
    double CSTop = crossSection(mCorrectedYield, effTop, lumi);

    //errors in percentages
    double low = (CSLow - getNominalCrossSection())/getNominalCrossSection();
    double top = (CSTop - getNominalCrossSection())/getNominalCrossSection();
    
    return make_pair(low, top);

}



void CrossSectionMaker::rapidityDependencePlot(int nBins, double low, double top,vector<double> plotEdges){
    
    mOutFile->mkdir("RapidityDependence");


    TGraphAsymmErrors *gRap = new TGraphAsymmErrors();
    gRap->SetName("rapidityGraph");


    for(unsigned int i = 0; i < plotEdges.size()-1 ; i++){
        TString condition = mAna->getCondition() + Form(" && pairRapidity > %.3f && pairRapidity < %.3f", plotEdges[i], plotEdges[i+1]);
        

        double correctedYield = getYieldCorrection(30, 0, 1.5, condition );
        cout << "Rapidity bin: " << i << ", Yield Correction: " << correctedYield << ", condition: " << condition << endl;
        cout << "Yield err top: " << mCorrYieldErrTop << ", yield err low: " << mCorrYieldErrLow << endl;
        
        mAna->DrawSTARInternal();
        mAna->CreateText(Form("%.1f < y < %.1f", plotEdges[i], plotEdges[i+1]), 62,0.7,0.75,0.85,0.77);
        saveYieldCorrection(Form("YieldCorrectionRap(%.1f,%.1f)", plotEdges[i], plotEdges[i+1]));

        double cs = crossSection(correctedYield, getRPEfficiency()*getVertexRecoEfficiency()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]);

        double errTop = crossSection(correctedYield+mCorrYieldErrTop, getRPEfficiency()*getVertexRecoEfficiency()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]) - cs;
        double errLow = cs - crossSection(correctedYield-mCorrYieldErrLow, getRPEfficiency()*getVertexRecoEfficiency()*mZeroBias->getTriggerEfficiency(), luminosity(), plotEdges[i+1]-plotEdges[i]);
        gRap->SetPoint(i, 0.5*(plotEdges[i]+plotEdges[i+1]), cs);
        gRap->SetPointError(i, 0.5*(plotEdges[i+1]-plotEdges[i]), 0.5*(plotEdges[i+1]-plotEdges[i]), errLow, errTop);



    }
    
    TCanvas *canvas = nullptr;
    mAna->CreateCanvas(&canvas, "RapidityDependence", 800, 600);

    TGraph* gMc = getMcNormalizedPlot(100, -1, 1);
    canvas->Clear();
    
    mAna->SetTGraphStyle(gRap, kBlue, 20);
    gRap->SetMarkerStyle(20);
    gRap->SetMarkerSize(2);
    gRap->SetLineColor(kBlue);
    gRap->SetLineWidth(2);
    gRap->SetLineColor(kBlue);
    gRap->GetYaxis()->SetRangeUser(0, 800);
    gRap->Draw("AEP");

    gMc->Draw("L");

    mAna->DrawSTARpp510JPsi(0.2, 0.82,0.4,0.9);

       
    TLegend *l;
    mAna->CreateLegend(&l, 0.6,0.6,0.89,0.89);
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
        cerr << "Could not find histogram histMc. Returning." << endl;
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
        cout << "Bin " << i << ": center = " << binCenter << ", content = " << content << ", new content = " << newContent << endl;

    }

    mEmbedding->SetTGraphStyle(gMc, kRed, 21);
    gMc->GetXaxis()->SetTitle("y [-]");
    gMc->GetYaxis()->SetTitle("d#sigma/dy [pb]");

    return gMc; 

}



