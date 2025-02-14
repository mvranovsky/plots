#include "../include/PlotAnaJPsi.h"


PlotAnaJPsi::PlotAnaJPsi(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotAnaJPsi::Make(){

	int nBins = 40;
    double lowerLim = 2;
    double upperLim = 4;
    roofit_invMassFit_Misa(nBins, lowerLim, upperLim);
    cout << "Finished Misa's plot" << endl;

    invMassJPsi(nBins, lowerLim, upperLim);
    cout << "Finished JPsi plot" << endl;


    hists1D = GetAllTH1D();
    if(hists1D.size() == 0){
        cerr << "Couldn't load 1D histograms from file. Leaving..." << endl;
        return;
    }
    for (int i = 0; i < hists1D.size(); ++i){
        if(!hists1D[i].first){
            cerr << "Couldn't load histogram " << hists1D[i].second << ". Leaving..." << endl;
            return;
        }
        TH1DGeneral(hists1D[i].second, hists1D[i].first);
        cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
    }

    // create canvases for all TH2F
    hists2D = GetAllTH2F();
    if(hists2D.size() == 0){
        cerr << "Couldn't load 2D histograms from file. Leaving..." << endl;
        return;
    }
    for (int i = 0; i < hists2D.size(); ++i){
        if(!hists2D[i].first){
            cerr << "Couldn't load histogram " << hists2D[i].second << ". Leaving..." << endl;
            return;
        }
        TH2FGeneral(hists2D[i].second, hists2D[i].first);
        cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
    }


    // get control plots and save them as histograms
    TH1D *hSignalControl = (TH1D*)inFile->Get("hInvMassJPsi");
    TH1D* hBcgControl = (TH1D*)inFile->Get("hInvMassJPsiBcg");

    if(hSignalControl && hBcgControl){
        outFile->cd();
        hSignalControl->Write("hInvMassJPsi_Control");
        hBcgControl->Write("hInvMassJPsiBcg_Control");
    }else{
        cout << "Could not get control invariant mass histograms." << endl;
    }


    // save pT plots
    TH1D *hPtMissing = (TH1D*)inFile->Get("hPtMissing");
    TH1D *hPtMissingBcg = (TH1D*)inFile->Get("hPtMissingBcg");
    TH1D *hPhotonMomX = (TH1D*)inFile->Get("hPhotonMomX");
    TH1D *hPhotonMomY = (TH1D*)inFile->Get("hPhotonMomY");
    TH1D *hPhotonMomXBcg = (TH1D*)inFile->Get("hPhotonMomXBcg");
    TH1D *hPhotonMomYBcg = (TH1D*)inFile->Get("hPhotonMomYBcg");

    outFile->cd();
    hPtMissing->Write(hPtMissing->GetName());
    hPtMissingBcg->Write(hPtMissingBcg->GetName());
    hPhotonMomX->Write(hPhotonMomX->GetName());
    hPhotonMomY->Write(hPhotonMomY->GetName());
    hPhotonMomXBcg->Write(hPhotonMomXBcg->GetName());
    hPhotonMomYBcg->Write(hPhotonMomYBcg->GetName());

    Clear();
    cout << "All histograms successfully saved to canvases..." << endl;
    cout << "The output file is saved: " << outputPosition << endl;
    

}

void PlotAnaJPsi::Init(){
	//define the output file which will store all the canvases
	outFile = unique_ptr<TFile>( TFile::Open(outputPosition, "recreate") );

	if(!outFile || outFile->IsZombie() ){
		cerr << "Couldn't open output file. Leaving..." << endl;
        return;
	}


	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaJPsiTree, true);

    if(!inFile || inFile->IsZombie() ){
        cerr << "Could not open input file. Leaving..." << endl;
        return;
    }

    if(!tree || !bcgTree){
    	cerr << "Couldn't open tree or background tree with data. Returning." << endl;
    	return;
    }
}


void PlotAnaJPsi::invMassJPsi(int numBins, Double_t minRange, Double_t maxRange) {
    
    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("invMassJPsiCanvas"), "Fit Result", 1200, 800); 
    SetGPad();
    c->cd();
    TH1D *signal = new TH1D("invMassJPsi", "invMassJPsi", numBins, minRange, maxRange);
    TH1D *bcg = new TH1D("invMassJPsiBcg", "invMassJPsiBcg", numBins, minRange, maxRange);
    signal->Sumw2();
    bcg->Sumw2();

    signal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    signal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    
    TString cmd, condition;
   	tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
   	signal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

   	bcgTree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
   	bcg->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

   	if(signal->GetEntries() == 0 || bcg->GetEntries() == 0 ){
   		cout << "Could not load signal or background to histograms" << endl;
   		return;
   	}
    signal->Sumw2();
    bcg->Sumw2();
    signal->Add(bcg, -1.0);

    RooRealVar x1("x1", "Observable1", minRange,maxRange);
    //RooRealVar x2("x2", "Observable2", minRange,maxRange);

    // Convert TH1D to RooDataHist
    RooDataHist data1("data1", "dataset with x1", x1, Import(*signal));
    //RooDataHist data2("data2", "dataset with x2", x2, Import(*bcg));

    // Define the Crystal ball PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", 3.09 , 2.5, 3.5); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", 0.05,0.0,10);
    RooRealVar alpha1("alpha1", "Tail parameter 1", 1.3, 1.,10);
    RooRealVar n1("n1", "Power parameter 1", 1, 1., 12);
    RooCBShape CBfunc1("cb1", "Crystal Ball PDF 1", x1, mean1, sigma1, alpha1, n1);

    // Define the polynomial PDF for the background
    RooRealVar coef1_0("coef1_0", "constant term", 10, -50, 50);
    RooRealVar coef1_1("coef1_1", "linear term", -6, -50, 50);
    RooRealVar coef1_2("coef1_2", "quadratic term", 1, -50, 50);
    RooPolynomial poly1("poly1", "Background Polynomial", x1, RooArgList(coef1_0, coef1_1, coef1_2)); //prehodit


    // Coefficients for PDFs
    RooRealVar nsig1("nsig1", "signal 1", 150, 0, 2000);
    RooRealVar nbkg1("nbkg1", "Background 1",110, 0, 2000);
    RooAddPdf model1("model1", "Crystal ball plus polynomial 1", RooArgList(CBfunc1, poly1), RooArgList(nsig1, nbkg1)); //,poly   gaussFrac2

    // Fit the model to the data
    RooFitResult *fitResult = model1.fitTo(data1, RooFit::Save());


    RooPlot* frame1 = x1.frame();
    //RooPlot* frame2 = x2.frame();
    frame1->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitleOffset(1.8);
    data1.plotOn(frame1);
    model1.plotOn(frame1);
    model1.plotOn(frame1, Components(poly1), LineStyle(kDashed), LineColor(kRed));
    model1.plotOn(frame1, Components(CBfunc1), LineStyle(kDashed), LineColor(kBlue));
    frame1->SetMinimum(-12);   // Set lower bound
    frame1->SetMaximum(72); // Set upper bound
    frame1->Draw();
    //frame2->Draw("same");



    // create a legend with results of fit
    TPaveText* text_1 = new TPaveText(0.6,0.8, 0.9, 0.6,"brNDC");
    text_1 -> SetTextSize(textSize);
    text_1 -> SetTextAlign(11);
    text_1 -> SetFillColor(0);
    text_1 -> SetTextFont(62);
    TString meanFinal1 = Form("#mu_{1} = %.3f +/- %.3f",mean1.getVal(), mean1.getError() );
    TString sigmaFinal1 = Form("#sigma_{1} = %.3f +/- %.3f",sigma1.getVal(), sigma1.getError() );
    TString fitFinal1 = Form("#chi^{2}/NDF = %.0f/%d #approx %.1f",frame1->chiSquare(), fitResult->floatParsFinal().getSize(), frame1->chiSquare()/fitResult->floatParsFinal().getSize() ); 
    TString yieldText1 = Form("Yield = %.0f +/- %.0f",nsig1.getVal(),nsig1.getError());
    /*
    TString pol1_0 = Form("Pol_{0} = " + to_string( makeInt(coef1_0.getVal()) );
    TString pol1_1 = "Pol_{1} = " + to_string( makeInt(coef1_1 .getVal()) ); 
    TString nText1 = "n_{1} = " + to_string(n1.getVal()) + "+/-" + to_string( n1.getError() );
    TString alphaText1 = "#alpha_{1} = " + to_string(alpha1.getVal() ) + "+/-" + to_string(alpha1.getError());*/
    //text_1 -> AddText( pol1_0 );
    //text_1 -> AddText( pol1_1 );
    text_1 -> AddText(meanFinal1);
    text_1 -> AddText(sigmaFinal1);
    //text_1 -> AddText( nText1 );
    //text_1 -> AddText( alphaText1);
    text_1 -> AddText(fitFinal1);
    text_1 -> AddText(yieldText1);
    text_1 -> Draw("same");
    
    DrawSTARpp510JPsi(0.5, 0.9, 0.9, 0.9);
    //CreateText("at least 1 ToF track",0.23, 0.87, 0.43, 0.87 );    
    TLegend *leg1;
    CreateLegend(&leg1, 0.25,0.6,0.5,0.4);
    leg1->SetTextAlign(11);
    leg1->SetTextFont(textFont);
    leg1->SetTextSize(0.04);
    int j = frame1->numItems();
    cout << "numItems: " << j << endl;
    TString names[j];
    TObject* obj[j];
    cout << "Frame objects:\n";
    j = 0;
    for (int i=0; i<frame1->numItems(); i++) {
        TString obj_name=frame1->nameOf(i); 
        if (obj_name=="") 
            continue;
        cout << Form("%d. '%s'\n",i,obj_name.Data());
        names[j] = Form("%d. '%s'\n",i,obj_name.Data());
        obj[j] = frame1->findObject(obj_name.Data());
        ++j;
    }
    leg1 -> AddEntry(obj[0], "Data", "lep");
    leg1 -> AddEntry(obj[1], "Fit function", "l");
    leg1 -> AddEntry(obj[2], "Crystal ball", "l");
    leg1 -> AddEntry(obj[3], "Poly2", "l");
    leg1 -> Draw("same");


    outFile->cd();
    c->Write("invMassJPsi");
    signal->Write("invMassJPsiSignal");
    bcg->Write("invMassJPsiBcg");

}

double PlotAnaJPsi::GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data){
    double chiSquare = frame->chiSquare();  // This calculates chi-square per degree of freedom

    int numBins = frame->GetNbinsX();  // Number of bins used in the histogram
    int nParams = model.getParameters(data)->getSize();  // Number of floating parameters in the model
    int ndf = numBins - nParams;

    return chiSquare/ndf;
}


TString PlotAnaJPsi::convertToString(double val) {

    ostringstream streamA;
    streamA << fixed << setprecision(1) << val;
    TString formattedA = streamA.str();

    return formattedA;
}

int PlotAnaJPsi::makeInt(double val) {
	int result = val;
	return result;
}



void PlotAnaJPsi::roofit_invMassFit_Misa(int numBins, Double_t minRange, Double_t maxRange){ 

    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("invMassJPsiCanvas_Misa"), "Fit Result", 1200, 800); 
    SetGPad();

    TH1D *hSignal = new TH1D("invMassJPsi_Misa", "invMassJPsi", numBins, minRange, maxRange);
    TH1D *hBcg = new TH1D("invMassJPsiBcg_Misa", "invMassJPsiBcg", numBins, minRange, maxRange);

    hSignal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
    
    TString cmd, condition;
    tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
    hSignal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

    bcgTree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "");
    hBcg->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

    if(hSignal->GetEntries() == 0 || hBcg->GetEntries() == 0 ){
        cout << "Could not load signal or background to histograms" << endl;
        return;
    }

    hSignal->Add(hBcg, -1);
    hSignal->SetTitle("");
    // Declare observable x
    RooRealVar x("x","m_{e^{+}e^{-}} [GeV/c^{2}]",2,4) ;

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*hSignal));
    RooPlot* frame = x.frame(Title("Uncorrected invariant mass of e^{+}e^{-} pairs")) ;
    dh.plotOn(frame,DataError(RooAbsData::SumW2),Name("data")); 

    //define function to fit background - polynomial 2 degree
    RooRealVar a0("a0","a0",1,-1,20); 
    RooRealVar a1("a1","a1",-1,-10,10);
    RooRealVar a2("a2","a2",0.5,-1,1);
    RooPolynomial bkg("bkg","bkg",x,RooArgList(a0,a1,a2));
  
    // define function to fit peak in data - crystal ball
    RooRealVar cbmean("mean","mean",3.1,2.8,3.2);
    RooRealVar cbsigma("sigma","sigma",5,0.01,10);
    //RooRealVar cbn("cb n","cb_n",6,0.1,10);
    RooRealVar cbn("cb n","cb_n",2,1,10);
    RooRealVar cbalpha("cb alpha","cb_alpha",2,1,10);
    RooCBShape cb("cb","cb",x,cbmean,cbsigma,cbalpha,cbn) ;

    //combine function into model
    RooRealVar bkgfrac("bkgfrac","bkgfrac",0.5,0.,1.);
    //RooAddPdf model("model","model",RooArgList(cb,bkg),bkgfrac);
    //RooAddPdf model("model","g+a",RooArgList(cb,bkg), bkgfrac) ;
    RooRealVar nsig("nsig","signal events",150,0,10000);
    RooRealVar nbkg("nbkg","signal background events",100,-10000,10000);
    

    RooAddPdf model("model","model",RooArgList(cb,bkg),RooArgList(nsig,nbkg));

    //fit data with model
    RooFitResult* fitResult = model.fitTo(dh, RooFit::Save()); 


    //plotting
    model.plotOn(frame, Name("model"));
    //bkg.plotOn(frame);
    model.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed), Name("background")); 
    frame->Draw("hist E");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleOffset(0.8);
    gStyle->SetOptTitle(0);

    cerr << "Before integration" << endl;

    //integration - raw yield
    x.setRange("signal",2.8682,3.2462);
    double i_s, i_b, y;
    RooAbsReal *intPeakX = model.createIntegral(x,NormSet(x),Range("signal"));   //under peak of pol+gaus
    RooAbsReal *intBKGx = bkg.createIntegral(x,NormSet(x), Range("signal")) ;     //(2.908, 3.261); 
    i_s = intPeakX->getVal()*(nsig.getVal()+nbkg.getVal()); 
    i_b = intBKGx->getVal()*(nbkg.getVal()); 
    cout << "i_s = " << i_s << "\n";
    cout << "i_b = " << i_b << "\n";
    y = (i_s - i_b);
    cout << "yield = " << y << "\n";



    //lines around mass +3-4 sigma
    auto l1 = new TLine(2.8682, -5, 2.8682, 17);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->Draw("same hist E");

    auto l2 = new TLine(3.2462, -5, 3.2462, 17);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->Draw("same hist E");


    TPaveText *text = new TPaveText(0.65,0.42,0.8,0.9, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
    text->SetTextSize(0.03);
    text->SetFillColor(0);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    text->AddText("p + p #rightarrow p + J/#psi + p");
    text->AddText("#sqrt{s} = 510 GeV, 2017");
    //text->AddText("p_{T} range (0.6-1.5)");
    //text->AddText("#Delta #Phi > 1.6 && JPSIpT < 1.5"); 
    text->AddText("After RP cuts");
    text->AddText("This work");
    text->AddText("");
    text->AddText(Form("#mu = %.4f #pm %.4f",cbmean.getVal(),cbmean.getError()));
    text->AddText(Form("#sigma = %.4f #pm %.4f",cbsigma.getVal(),cbsigma.getError()));
    text->AddText(Form("#chi^{2}/NDF = %.0f/%d #approx %.2f",frame->chiSquare(), fitResult->floatParsFinal().getSize(), frame->chiSquare()/fitResult->floatParsFinal().getSize() ) ); 
    text->AddText(Form("Raw yield = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
    text->Draw("same hist E");



    TLegend *leg1 = new TLegend(0.25,0.68,0.4,0.9);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->AddEntry("data","Data", "LEP");
    leg1->AddEntry("model","Crystal Ball + Poly2","LP");
    leg1->AddEntry("background","Poly2", "LP");
    leg1->Draw("same hist E");


    outFile->cd();
    c->Write("invMassJPsi_MisaPlotting");
}
