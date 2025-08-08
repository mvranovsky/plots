#include "CrossSectionMaker.h"

CrossSectionMaker::CrossSectionMaker(const TString AnaDir, const TString EmbeddingDir, const TString GoodRunDir, TString outputFile) {

    mAnaDir = AnaDir;
    mEmbeddingDir = EmbeddingDir;
    mGoodRunDir = GoodRunDir;
    mOutFileName = outputFile;

    // initialize utility
    mUtil = new Util();


    // connect output file
	mOutFile = unique_ptr<TFile>(new TFile( mOutFileName.Data(), "WRITE") );
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

    // initialize ProbRetainEvent class
    cout << "Initializing ProbRetainEvent class..." << endl;
    TTree *mGoodRunTreeBcg = nullptr;
    mGoodRun->ConnectInputTree(GoodRunDir.Data(), nameOfAnaGoodRunTree.Data(), mGoodRunTree, mGoodRunTreeBcg);

    if(!mGoodRunTree || mGoodRunTree->GetEntries() == 0){
        cerr << "Could not open good run tree. Returning." << endl;
        return;
    }
    
    mProbRetainEvent = new ProbRetainEvent(mOutFile, mGoodRunTree);

}


void CrossSectionMaker::Make(){

    mAna->Make();

    mEmbedding->Make();

    mGoodRun->Make();

    mProbRetainEvent->Make();

    // calculate cross section
    double crossSection = CalculateCrossSection();

    CalculateSysError(NHITSFIT);
    CalculateSysError(NHITSDEDX);
    CalculateSysError(ETA);
    CalculateSysError(PID);

    /*
    CalculateSysError(VERTEXZ);
    CalculateSysError(DCAZINCM);
    CalculateSysError(DCAXYINCM);
    */
}


double CrossSectionMaker::CalculateCrossSection(){

    // get yields from data
    double yield = mAna->getYieldFinal();
    cout << "Final yield: " << yield << endl;

    // get reconstruction efficiency from embedding
    double recoEff = mEmbedding->getEfficiencyFinal();
    cout << "Reco eff: " << recoEff << endl;

    // get trigger efficiency from good run
    double triggerEff = getTopologyEfficiency();
    cout << "Trigger eff: " << triggerEff << endl;

    double RPEff = getRPEfficiency();
    cout << "RP eff: " << RPEff << endl;

    double lumi;
    if(isRomanPots){
        lumi = mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), mProbRetainEvent->getA(), mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }else{
        lumi = mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), mProbRetainEvent->getA(), mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }
    cout << "Lumi: " << lumi << " pb^-1" << endl;

    double deltaRap = 2;

    double eff = recoEff*RPEff*triggerEff;
    double CS = yield/(deltaRap*lumi*eff);
    cout << "Final cross section: " << CS << " pb" << endl;

    return CS;

}

void CrossSectionMaker::CalculateSysError(int VAR){

    cout << "Running systematic analysis of cut " << mUtil->nameOfVariable(VAR) << endl;

    // get yields from data
    double yieldNominal = mAna->getYield(VAR, NOMINAL);
    double yieldTight = mAna->getYield(VAR, TIGHT);
    double yieldLoose = mAna->getYield(VAR, LOOSE);

    cout << "Got the yields from analysis" << endl;
    // get reconstruction efficiency from embedding
    double recoEffNominal = mEmbedding->getEfficiency(VAR, NOMINAL);
    double recoEffTight = mEmbedding->getEfficiency(VAR, TIGHT);
    double recoEffLoose = mEmbedding->getEfficiency(VAR, LOOSE);
    cout << "Got the yields for sys study " << mUtil->nameOfVariable(VAR) << endl;

    // get trigger efficiency from good run
    double triggerEff = getTopologyEfficiency();
    cout << "Topology efficiency: " << triggerEff << endl;
    
    double RPEff = getRPEfficiency();
    cout << "RP eff: " << RPEff << endl;
    
    double lumi;
    if(isRomanPots){
        lumi = mProbRetainEvent->getLuminosity( mGoodRun->getRpGoodRunList(), mProbRetainEvent->getA(), mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }else{
        lumi = mProbRetainEvent->getLuminosity( mGoodRun->getGoodRunList(), mProbRetainEvent->getA(), mProbRetainEvent->getB() ); // A = 1.0, B = 0.0
    }
    cout << "Lumi: " << lumi << endl;

    double deltaRap = 2;

    double effNominal = recoEffNominal*RPEff*triggerEff;
    double CSNominal = yieldNominal/(deltaRap*lumi*effNominal);
    cout << "Nominal cross section: " << CSNominal << endl;

    double effTight = recoEffTight*RPEff*triggerEff;
    double CSTight = yieldTight/(deltaRap*lumi*effTight);
    mSysErrorTight[VAR] = abs(CSTight - CSNominal)/CSNominal;
    cout << "Tight cross section: " << CSTight << ". Systematic error: " << mSysErrorTight[VAR] << endl;


    double effLoose = recoEffLoose*RPEff*triggerEff;
    double CSLoose = yieldLoose/(deltaRap*lumi*effLoose);
    mSysErrorLoose[VAR] = abs(CSLoose - CSNominal)/CSNominal;
    cout << "Loose cross section: " << CSLoose << ". Systematic error: " << mSysErrorLoose[VAR] << endl;
}











CrossSectionMaker::~CrossSectionMaker() {
    if(mUtil) delete mUtil;
    if(mProbRetainEvent) delete mProbRetainEvent;
    if(mEmbedding) delete mEmbedding;
    if(mAna) delete mAna;
    if(mGoodRun) delete mGoodRun;
}









/*
vector<double> CrossSectionMaker::correctedYield(){
    double low = cbmean.getVal() - 4.0 * cbsigma.getVal();
    double high = cbmean.getVal() + 3.0 * cbsigma.getVal();
    
    TString condition = Form("invMass > %f && invMass < %f", low, high);
    tree->Draw("pt>>hist(60,0,1.5)", condition);
    TH1* hPt = (TH1*)gPad->GetPrimitive("hist");
    if(!hPt || hPt->GetEntries() == 0){
        cout << "Could not get pt histogram. Returning." << endl;
        return {0};
    }

    TFile *recEffFile = TFile::Open("AnaJPsiRecoEff.root", "READ");
    if(!recEffFile || recEffFile->IsZombie() || !recEffFile->IsOpen()){
        cout << "Could not open reconstruction efficiency file. Returning." << endl;
        return {0}; 
    }

    TGraphAsymmErrors* recEff = dynamic_cast<TGraphAsymmErrors*>(recEffFile->Get("recoEff_pTJPsi"));
    if(!recEff){
        cout << "Could not get reconstruction efficiency graph. Returning." << endl;
        recEffFile->Close();
        return {0};
    }

    if(recEff->GetN() != hPt->GetNbinsX()){
        cout << "Reconstruction efficiency graph has different number of bins than pT histogram. Returning." << endl;
        recEffFile->Close();
        return {0};
    }

    double fractionSignal = signal / (signal + background);
    
    double correctedYield = 0;
    double correctedError = 0;
    vector<double> vals;
    for(int i = 1; i <= hPt->GetNbinsX(); i++){
        double binContent = hPt->GetBinContent(i);
        double binError = hPt->GetBinError(i);

        // Get efficiency for this bin
        double eff = recEff->GetY()[i-1];
        double effError = (recEff->GetEYlow()[i-1] + recEff->GetEYhigh()[i-1]) / 2;
        if(eff > 0){
            // Calculate corrected yield for this bin
            reconstructionEff += eff;
            reconstructionEfficiencyError += pow(effError,2);
            vals.push_back(eff);
            correctedYield += binContent / eff * fractionSignal;
            correctedError += pow(binError / eff * fractionSignal, 2);

        }else{
            cout << "Efficiency is zero for bin " << i << ". Major problem!" << endl;
            return {0};
        }
    }
    
    
    correctedError = sqrt(correctedError);
    reconstructionEff /= hPt->GetNbinsX();
    
    double standardDeviation = 0.0;
    for(unsigned int i = 0; i < vals.size(); ++i){
        standardDeviation += pow(reconstructionEff  - vals[i], 2);
    }
    standardDeviation = sqrt(standardDeviation / vals.size());
    
    double varianceFromSpread = pow(standardDeviation, 2) / vals.size();
    reconstructionEfficiencyError = sqrt(reconstructionEfficiencyError + varianceFromSpread);
    
    return {correctedYield, correctedError};
}
*/


