#include "FitJPsi.h"


FitJPsi::FitJPsi(TH1D *&h,TString b){
   hist = (TH1D*)h->Clone("hSignal");
   bcg = b;

   fitRangeLow = hist->GetXaxis()->GetXmin();
   fitRangeHigh = hist->GetXaxis()->GetXmax();

   RooMsgService::instance().setSilentMode(true);
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
}

void FitJPsi::fitPeak() {

   // create a canvas that will hold both fits
   c = new TCanvas("peakFit", "Fit Result", 1200, 800); 

   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.03);
   gPad->SetBottomMargin(0.11);
   gPad->SetTickx();
   gPad->SetTicky(); 



   gStyle->SetOptStat("");   

   if(!hist || hist->GetEntries() == 0){
      cout << "ERROR in FitJPsi::fitPeak - Empty or no hist when trying to fit. Leaving." << endl;
      return;
   }

   // Set axis title
   hist->GetXaxis()->SetTitle("m_{ee} [GeV]");
   hist->SetTitle("");


   // Declare observable x
   x = new RooRealVar("x", "m_{ee} [GeV/c^{2}]", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

   // Create a binned dataset from histogram
   dh = new RooDataHist("dh", "dh", *x, Import(*hist));

   // Create a frame for plotting
   frame = x->frame(Title(""));

   dh->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));

   // Define the components of the model based on the chosen background function
   fitSignal = false, fitSignal2S = false;
   model = buildModel(fitBcg, fitSignal, fitSignal2S);


   // Check if model was created successfully before fitting
   if (!model) {
      std::cerr << "ERROR: FitJPsi::buildModel - Model was not initialized. Check option string: " << bcg << std::endl;
      return;
   }
   if(DEBUGMODE) cout << "Defining fit range..." << endl;
   defineFitRange();

   // Fit model to data
   if(DEBUGMODE) cout << "Starting fit..." << endl;
   fitResult = model->fitTo(*dh, PrintLevel(1), RooFit::Range("fitRange"), RooFit::Minimizer("Minuit2"), Strategy(2)/*, Hesse(true)*/, Save() );
   //fitResult = model->fitTo(*dh, Save(), PrintLevel(-1),RooFit::Verbose(false), RooFit::Warnings(false), RooFit::Range("fitRange"));

   if (!fitResult) {
      std::cerr << "ERROR: FitJPsi::buildModel - Fit result is null. Check if the model was created successfully." << std::endl;
      return;
   }
   if(DEBUGMODE) cout << "Fit completed." << endl;
   if(DEBUGMODE) fitResult->Print("v");

   // Print fit parameters
   // Calculate degrees of freedom
   NDF = hist->GetNbinsX() - (fitResult->floatParsFinal().getSize());

   //plotting
   model->plotOn(frame, Name("model"));
   if(fitBcg && visualizeError) {
      model->plotOn(frame, Name("bkgError"), RooFit::Range("fitRange"),VisualizeError(*fitResult, 1, false), Components("bkg"), FillStyle(3001), FillColor(kRed), DrawOption("F"));
      model->plotOn(frame, Components("bkg"), RooFit::Range("fitRange"),LineStyle(kDashed), LineColor(kRed), Name("bkg"));
   }else if(fitBcg ) {
      model->plotOn(frame, Components("bkg"), RooFit::Range("fitRange"),LineStyle(kDashed), LineColor(kRed), Name("bkg"));
   }
\
   if(fitSignal && visualizeError){
      model->plotOn(frame, RooFit::Range("fitRange"),Components("signal1S"), DrawOption("F"), LineStyle(kDashed), LineColor(kGreen), Name("signal"));
      model->plotOn(frame, Name("signalError"), RooFit::Range("fitRange"),VisualizeError(*fitResult, 1, false), Components("signal1S"), FillStyle(3001), FillColor(kGreen), DrawOption("F"));
   }else if(fitSignal){
      model->plotOn(frame, RooFit::Range("fitRange"),Components("signal1S"), LineStyle(kDashed), LineColor(kGreen), Name("signal"));
   }

   if(fitSignal2S && visualizeError){
      model->plotOn(frame, RooFit::Range("fitRange"),Components("signal2S"), DrawOption("F"), LineStyle(kDashed), LineColor(kMagenta), Name("signal2"));
      model->plotOn(frame, Name("signal2Error"), RooFit::Range("fitRange"),VisualizeError(*fitResult, 1, false), Components("signal2S"), FillStyle(3001), FillColor(kMagenta), DrawOption("F"));
   }else if(fitSignal2S){
      model->plotOn(frame, RooFit::Range("fitRange"),Components("signal2S"), LineStyle(kDashed), LineColor(kMagenta), Name("signal2"));
   }

   frame->Draw("hist E");
   frame->GetXaxis()->SetTitleSize(0.05);
   frame->GetYaxis()->SetTitleSize(0.05);
   frame->GetXaxis()->SetTitleOffset(0.8);
   frame->GetYaxis()->SetTitleOffset(0.8);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   if(DEBUGMODE) cout << "After fit plotting." << endl;
   
   leg1 = new TLegend(legPosX1, legPosY1, legPosX2, legPosY2);
   leg1->SetTextSize(0.03);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(0);
   if(mIsEmbedding){
      leg1->AddEntry("data", mEmbeddingDescription, "LEP");
   }else{
      leg1->AddEntry("data","Data", "LEP");
   }
   if(!(!fitSignal && !fitSignal2S && fitBcg)) leg1->AddEntry("model","Fit Model", "L"); // only add if at least one signal is fitted
   if(fitSignal) leg1->AddEntry("signal", modelNames[0], "L");
   if(fitBcg) leg1->AddEntry("bkg", modelNames[1], "L");
   if(fitSignal2S) leg1->AddEntry("signal2", modelNames[2], "L");
   //leg1->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), NDF, frame->chiSquare()/NDF ), "");
   
   if(showAdditionalContinuum) loadRooFitResult("fitResultRoofit.root");
   if(mShowDataPeak) loadDataPeak(mDataMean, mDataSigma);
   frame->Draw("hist F E");
   leg1->Draw("same");
   
   if(fitSignal) integrate(bcg);
   
   if(fitSignal)  writeSignalResult(fitSignal2S,0.7,0.58,0.85,0.78);
   
   if(fitSignal) drawLines();
   

   c->Update();
}

FitJPsi::~FitJPsi(){
   if(c) {
       c->Close();
       delete c;
   }
   if(fitResult) delete fitResult;
   if(model) delete model;
   if(frame) delete frame;
   if(hSignal) delete hSignal;
}

void FitJPsi::defineFitRange(){

   x->setRange("fitRange",fitRangeLow, fitRangeHigh);
}



RooAddPdf* FitJPsi::buildModel(bool &fitBcg, bool &fitSignal, bool &fitSignal2S){

   if(DEBUGMODE) cout << "Building model with option: " << bcg << endl;

   models.clear();
   modelNames.clear();
   fitBcg = false;
   // ----------------------------------------------------------------------------------------------
   if(bcg.Contains("cb") || bcg.Contains("CB") || bcg.Contains("CrystalBall") || bcg.Contains("crystalball") ){
      // Define signal model (Crystal Ball)
      cbmean   = new RooRealVar("mean", "mean", 3.0908, 2.8, 3.2);
      cbsigma  = new RooRealVar("sigma", "sigma", 0.050, 0.0, 0.0595);

      if(mAlphaLoose){  // for fitting alpha
         cout << "Setting alpha loose for Crystal Ball fit." << endl;
         cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41, 0.5, 5.0); // loose
      }else{            // for choosing specific alpha
         cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41); // fixed
         cbalpha->setConstant(true);
      }
      if(mNLoose){      // for fitting n
         cout << "Setting n loose for Crystal Ball fit." << endl;
         cbn      = new RooRealVar("cb_n", "cb_n", 1.93, 0.5, 5.0);         // loose
      }else{         // for choosing specific n 
         cbn      = new RooRealVar("cb_n", "cb_n", 1.93);
         cbn->setConstant(true);
      }

      cb = new RooCBShape("signal1S", "cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
      models.push_back(cb);
      modelNames.push_back("Crystal Ball (J/#psi)");
      fitSignal = true;
   }else if(bcg.Contains("gauss") || bcg.Contains("Gauss") || bcg.Contains("Gaussian") || bcg.Contains("gaussian")){
      // Define signal model (Gaussian)
      cbmean   = new RooRealVar("mean", "mean", 3.0908, 2.8, 3.2);
      cbsigma  = new RooRealVar("sigma", "sigma", 0.0543, 0.0, 0.065);
      gauss = new RooGaussian("signal1S", "gauss", *x, *cbmean, *cbsigma);
      models.push_back(gauss);
      modelNames.push_back("Gauss (J/#psi)");
      fitSignal = true;
   }else{ modelNames.push_back(""); }
   // ----------------------------------------------------------------------------------------------
   if (bcg.Contains("poly1") || bcg.Contains("Poly1")) {
      a0 = new RooRealVar("a0", "a0", 5, -10, 200);
      // include observable x in the variable list so expression 'a0' resolves in RooGenericPdf
      bkg = new RooGenericPdf("bkg", "bkg", "a0", RooArgList(*a0));
      models.push_back(bkg);
      modelNames.push_back("Polynomial 1");
      fitBcg = true;
   }else if (bcg.Contains("poly2") || bcg.Contains("Poly2"))  {
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      a1 = new RooRealVar("a1", "a1", -1, -200, 200);
      // include observable x in variable list so 'x' is recognized
      bkg = new RooGenericPdf("bkg", "bkg", "a0 + a1*x", RooArgList(*x, *a0, *a1));
      models.push_back(bkg);
      modelNames.push_back("Polynomial 2");
      fitBcg = true;
   }else if (bcg.Contains("poly3") || bcg.Contains("Poly3")) {
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      a1 = new RooRealVar("a1", "a1", -1, -200, 200);
      a2 = new RooRealVar("a2", "a2", 0.5, -200, 200);
      // include observable x so polynomial terms using x are valid
      bkg = new RooGenericPdf("bkg", "bkg", "a0 + a1*x + a2*x^2", RooArgList(*x, *a0, *a1, *a2));
      models.push_back(bkg);
      modelNames.push_back("Polynomial 3");
      fitBcg = true;
   }else if(bcg.Contains("continuum") || bcg.Contains("Continuum") || bcg.Contains("cont") ){ //custom function to fit the gamma+gamma -> e+ e- continuum background
      c1 = new RooRealVar("c1", "c1 parameter", 0.4, 0.0, 0.49);  //0.419
      c2 = new RooRealVar("c2", "c2 parameter", -0.4, -3.0, 3.0); //-0.01
      c3 = new RooRealVar("c3", "c3 parameter",0.03, -3.0, 1.0); //-0.45
      c4 = new RooRealVar("c4", "c4 parameter",0.00001, 0.0, 0.0001); //-0.45
      c4->setConstant(true);

      // Define the custom PDF. Include observable 'x' in the RooArgList so expressions referencing x compile.
      bkg = new RooGenericPdf("bkg", "Custom function", "sqrt((x - c1)^2 + c4) * exp( c2 *(x - c1)^2 + c3*x^3)", RooArgList(*x, *c1, *c2, *c3, *c4));
      models.push_back(bkg);
      modelNames.push_back("Continuum");
      fitBcg = true;
   } else{ modelNames.push_back(""); }
   // ----------------------------------------------------------------------------------------------
   if(bcg.Contains("2SG") || bcg.Contains("2sg")) {
      // Define signal model (Double Gaussian for 2S state)
      mean2S   = new RooRealVar("mean2S", "mean2S", 3.686, 3.5, 3.9);
      sigma2S  = new RooRealVar("sigma2S", "sigma2S", 0.0543, 0.0, 0.065);
      gauss2S = new RooGaussian("signal2S", "gauss2S", *x, *mean2S, *sigma2S);
      models.push_back(gauss2S);
      modelNames.push_back("Gauss (#psi (2S))");
      fitSignal2S = true;

   }else if(bcg.Contains("2sc") || bcg.Contains("2SC")) {
      // Define signal model (Double Gaussian for 2S state)
      mean2S   = new RooRealVar("mean2S", "mean2S", 3.686, 3.5, 3.9);
      sigma2S  = new RooRealVar("sigma2S", "sigma2S", 0.0543, 0.0, 0.065);
           // for fitting n
      cbn2S      = new RooRealVar("cb_n", "cb_n", 1.93, 0.5, 5.0);         // loose
      cbalpha2S  = new RooRealVar("cb_alpha", "cb_alpha", 1.41, 0.5, 5.0); // loose

      cb = new RooCBShape("signal2S", "cb", *x, *mean2S, *sigma2S, *cbalpha2S, *cbn2S);
      models.push_back(cb);
      modelNames.push_back("Crystal Ball (#psi(2S))");
      fitSignal2S = true;
   } else{ modelNames.push_back(""); }

   if(models.size() == 0 || models.size() > 3){
      cout << "ERROR: FitJPsi::buildModel - Invalid number of models. One can select at most 2 models, one signal (cb, gauss) and background (poly, continuum)." << endl;
      return nullptr;
   }


   RooAddPdf* mod;
   // Combine signal and background into model
   if(fitSignal && fitBcg && fitSignal2S) {
      nsig = new RooRealVar("nsig", "signal events", 400, 0, 20000);
      nsig2S = new RooRealVar("nsig2S", "2S signal events", 10, 0, 10000);
      nbkg = new RooRealVar("nbkg", "background events", 700, 0, 1000000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0], *models[1], *models[2]), RooArgList(*nsig , *nbkg,*nsig2S));
   } else if(fitBcg && fitSignal2S) {
      nsig2S = new RooRealVar("nsig2S", "2S signal events", 10, 0, 10000);
      nbkg = new RooRealVar("nbkg", "background events", 700, 0, 1000000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0], *models[1]), RooArgList( *nbkg, *nsig2S));
   } else if (fitSignal && fitBcg) {
      nsig = new RooRealVar("nsig", "signal events", 400, 0, 20000);
      nbkg = new RooRealVar("nbkg", "background events", 70, 0, 1000000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0], *models[1]), RooArgList(*nsig, *nbkg));
   } else if(fitSignal){
      nsig = new RooRealVar("nsig", "signal events", 1000, 0, 15000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0]), RooArgList(*nsig));
   }else if(fitBcg) {   
      nbkg = new RooRealVar("nbkg", "background events", 10, 0, 15000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0]), RooArgList(*nbkg));
   }else {
      cerr << "ERROR: FitJPsi::buildModel - No valid model components selected. Please choose at least one signal or background model." << endl;
      return nullptr;
   }

   return mod;
}

void FitJPsi::printCovarianceMatrix(){
   fitResult->covarianceMatrix().Print();
}


void FitJPsi::integrate(TString bcg){
   // 1) Define asymmetric integration window
   double low = cbmean->getVal() - 4.0 * cbsigma->getVal();
   double high = cbmean->getVal() + 3.0 * cbsigma->getVal();


   x->setRange("signalRange", low, high);
   if(DEBUGMODE) cout << "Integration range: " << low << " to " << high << endl;

   // 2) Calculate fractions inside the range for signal and background PDFs
   RooAbsReal* intSignal = cb->createIntegral(*x, NormSet(*x), RooFit::Range("signalRange"));
   double fracBackground = 0.0;
   RooAbsReal* intBackground = nullptr;
   if(fitBcg){
      intBackground = bkg->createIntegral(*x, NormSet(*x), RooFit::Range("signalRange"));
      fracBackground = intBackground->getVal();
      RooFormulaVar bkgYieldFormula("bkgYieldFormula", "@0 * @1", RooArgList(*nbkg, *intBackground));
      errBackground = bkgYieldFormula.getPropagatedError(*fitResult);
   }

   double fracSignal = intSignal->getVal();
   if(DEBUGMODE) cout << "Signal fraction in range: " << fracSignal << endl;
   if(DEBUGMODE) cout << "nsig: " << nsig->getVal() << endl;
   
   // 3) Calculate yields in the range (mu -4sigma, mu + 3sigma) by multiplying by the normalization parameters
   yieldSignal = nsig->getVal() * fracSignal;
   yieldBackground = (fracBackground > 0) ? nbkg->getVal() * fracBackground : 0.0;
   
   // Construct formula variables representing yields in the signal range
   RooFormulaVar sigYieldFormula("sigYieldFormula", "@0 * @1", RooArgList(*nsig, *intSignal));
   errSignal = sigYieldFormula.getPropagatedError(*fitResult);
   
   // 5) Calculate net yield and its error
   netYield = yieldSignal;
   errNetYield = errSignal;
   
   // 6) Calculate signal-to-background ratio and its error
   sOverB = (yieldBackground > 0) ? yieldSignal / yieldBackground : 0;
   
   errSOverB = 0;
   if (yieldBackground > 0 && yieldSignal > 0) {
      errSOverB = sOverB * sqrt(
         pow(errSignal / yieldSignal, 2) +
         pow(errBackground / yieldBackground, 2)
      );
   }
   
}



void FitJPsi::saveCanvas(TFile *&file, TString name ,TString dir){
   
   if(!file || file->IsZombie()){
      cerr << "ERROR in FitJPsi::saveCanvas - File is not open or is a zombie. Cannot save canvas." << endl;
      return;
   }

   if(name !=  ""){
      c->SetName(name);
   }

   
   file->cd();
   if(dir == ""){
       c->Write();
   }else{
       file->cd(dir);
       c->Write();
   }

}

void FitJPsi::saveCanvas(shared_ptr<TFile> &file, TString name ,TString dir){
   if(!file || file->IsZombie()){
      cerr << "File is not open or is a zombie. Cannot save canvas." << endl;
      return;
   }
   
   if(name !=  ""){
      c->SetName(name);
   }
   

   file->cd();
   if(dir == ""){
      c->Write();
   }else{
       file->cd(dir);
       c->Write();
   }
   
}


void FitJPsi::drawLines(){
   
   // find the maximum bin value of the histogram
   double maxVal = hist->GetMaximum();
   
   TLine *line = new TLine(getLowLimitFit(), 0, getLowLimitFit(), maxVal/2);
   line->SetLineColor(kBlack);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw("same");
   
   TLine *line2 = new TLine(getHighLimitFit(), 0, getHighLimitFit(), maxVal/2);
   line2->SetLineColor(kBlack);
   line2->SetLineStyle(2);
   line2->SetLineWidth(2);
   line2->Draw("same");
   
}

void FitJPsi::writeContinuumResult(double left, double bottom, double right, double top){

   TPaveText *cText = new TPaveText(left, bottom, right, top, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) 
   cText->SetTextSize(0.03);
   cText->SetFillStyle(0);
   cText->SetBorderSize(0);
   cText->SetFillColorAlpha(0, 0);
   cText->SetTextFont(42);
   // Set text alignment: horizontal=1 (left), vertical=3 (top) => 10*1 + 3 = 13
   cText->SetTextAlign(11);
   cText->AddText("f_{bkg}(x) = (x - c_{1}) exp(c_{2}(x - c_{1})^{2} + c_{3} x^{3})");
   cText->AddText(Form("c_{1} = %.2f #pm %.2f", c1->getVal(), c1->getError()));
   cText->AddText(Form("c_{2} = %.2f #pm %.2f", c2->getVal(), c2->getError()));
   cText->AddText(Form("c_{3} = %.3f #pm %.3f", c3->getVal(), c3->getError()));
   //cText->AddText(Form("c_{4} = %.4f #pm %.4f", c4->getVal(), c4->getError()));
   cText->Draw("same");

}

void FitJPsi::writeFitResult(){
   leg1->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare()*NDF, NDF, frame->chiSquare() ), "");
}


void FitJPsi::writeSignalResult(bool &fitSignal2S,double left, double bottom, double right, double top){
   if(fitSignal2S) bottom -= 0.12;
   TPaveText *text = new TPaveText(left, bottom, right, top, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
   text->SetTextSize(0.03);
   text->SetFillStyle(0);
   text->SetBorderSize(0);
   text->SetFillColorAlpha(0, 0);
   text->SetTextFont(42);
   text->SetTextAlign(12);
   text->AddText("");
   if(fitSignal && !fitSignal2S){
      text->AddText(Form("#mu = %.4f #pm %.4f",cbmean->getVal(),cbmean->getError()));
      text->AddText(Form("#sigma = %.5f #pm %.5f",cbsigma->getVal(),cbsigma->getError()));
      text->AddText(Form("Raw yield (J/#psi)= %.0f #pm %.0f", round(netYield), round( (errNetYield)) ));
   }else if(fitSignal && fitSignal2S){
      text->AddText(Form("#mu_{J/#psi} = %.4f #pm %.4f",cbmean->getVal(),cbmean->getError()));
      text->AddText(Form("#sigma_{J/#psi} = %.5f #pm %.5f",cbsigma->getVal(),cbsigma->getError()));
      text->AddText(Form("Raw yield (J/#psi)= %.0f #pm %.0f", round(netYield), round( (errNetYield)) ));
      if(fitBcg)  text->AddText(Form("S/B (J/#psi) = %.0f / %.0f #approx %.1f", yieldSignal, yieldBackground, sOverB ) );
      text->AddText(Form("#mu_{#psi(2S)} = %.4f #pm %.4f",mean2S->getVal(),mean2S->getError()));
      text->AddText(Form("#sigma_{#psi(2S)} = %.5f #pm %.5f",sigma2S->getVal(),sigma2S->getError()));
      text->AddText(Form("Raw yield (#psi (2S))= %.0f #pm %.0f", round(nsig2S->getVal()), round(nsig2S->getError()) ));
   }

   text->Draw("same");
}

void FitJPsi::saveRooFitResult(TString filename){
   if(!fitResult){
      cerr << "Fit result is null. Cannot save." << endl;
      return;
   }
   TFile f(filename, "RECREATE");
   fitResult->Write("fitResult");
   f.Close();
   cout << "Fit result saved to " << filename << endl;
}

void FitJPsi::loadRooFitResult(TString filename) {
    TFile f(filename, "READ");
    if (f.IsZombie()) {
        cerr << "File " << filename << " is a zombie. Cannot load fit result." << endl;
        return;
    }

    RooFitResult* loadedRes = dynamic_cast<RooFitResult*>(f.Get("fitResult"));
    if (!loadedRes) {
        cerr << "No fit result found in " << filename << endl;
        return;
    }

    // --- Extract parameter values from RooFitResult ---
    double c1_val = loadedRes->floatParsFinal().getRealValue("c1");
    double c2_val = loadedRes->floatParsFinal().getRealValue("c2");
    double c3_val = loadedRes->floatParsFinal().getRealValue("c3");

    // Define parameters. CRITICAL: These must be RooRealVar, not constants, 
    // IF you want VisualizeError to work on them. They hold the final values.
    RooRealVar* c1 = new RooRealVar("c1", "c1 parameter", c1_val);
    RooRealVar* c2 = new RooRealVar("c2", "c2 parameter", c2_val);
    RooRealVar* c3 = new RooRealVar("c3", "c3 parameter", c3_val);
    // c4_const is fine as it was fixed in the fit.
    RooRealVar* c4_const = new RooRealVar("c4_const", "c4 parameter", 0.00001);
    c4_const->setConstant(true);

    // --- Build continuum PDF ---
    // Use the *c1, *c2, *c3 parameters which match the names in the fit result
    RooGenericPdf* continuumFunc = new RooGenericPdf(
        "continuumFunc", "Custom function",
        "sqrt((x - c1)^2 + c4_const) * exp(c2*(x - c1)^2 + c3*x^3)",
        RooArgList(*x, *c1, *c2, *c3, *c4_const) // Pass the parameters, not the constants
    );

    // IMPORTANT: nbkg must be a RooRealVar initialized with the final value
    // If nbkg is globally defined as a RooRealVar, ensure it is set to the final value here
    double nbkg_val = loadedRes->floatParsFinal().getRealValue("nbkg");
    nbkg->setVal(nbkg_val); 
    
    RooAddPdf* modelWithContinuum = new RooAddPdf(
        "modelWithContinuum", "model with continuum",
        RooArgList(*continuumFunc), RooArgList(*nbkg)
    );

    const char* contRangeName = "contPlotRange";
    x->setRange(contRangeName, c1->getVal(), fitRangeHigh); // Use c1->getVal() for the lower bound

    // Plot the error band
    modelWithContinuum->plotOn(frame,
        RooFit::VisualizeError(*loadedRes, 1),
        RooFit::Range(contRangeName),
        RooFit::NormRange(contRangeName),
        RooFit::FillColor(kMagenta),
        RooFit::FillStyle(3001),
        RooFit::DrawOption("F"),
        RooFit::Name("modelWithContError")
    );

    // Plot the central line on top (essential for color in legend)
    modelWithContinuum->plotOn(frame,
        RooFit::Range(contRangeName),
        RooFit::NormRange(contRangeName),
        RooFit::LineColor(kMagenta),
        RooFit::Name("modelWithCont_line")
    );
    leg1->AddEntry("modelWithCont_line", "Fit to data (scaled)", "L");

    f.Close();
    
}

void FitJPsi::loadDataPeak(double mean, double sigma){

   datamean   = new RooRealVar("mean", "mean", mean, 2.8, 3.2);
   datamean->setConstant(true);
   datasigma  = new RooRealVar("sigma", "sigma", sigma, 0.0, 0.065);
   datasigma->setConstant(true);
   dataalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41); // fixed
   dataalpha->setConstant(true);
   datan      = new RooRealVar("cb_n", "cb_n", 1.93);
   datan->setConstant(true);


   datacb = new RooCBShape("signalData", "cb", *x, *datamean, *datasigma, *dataalpha, *datan);
   datamodel = new RooAddPdf("adddatamodel", "adddatamodel", RooArgList(*datacb), RooArgList(*nsig));
   
   datamodel->plotOn(frame,RooFit::Range("fitRange"),Components("signalData"), Name("addedDataPeak"), LineColor(kRed));
   RooCurve* curve = (RooCurve*)frame->findObject("addedDataPeak");

   // Add to legend by pointer, not by string
   if (curve) leg1->AddEntry(curve, "Data (scaled)", "L");   
}