#include "Fit.h"


Fit::Fit(TH1D *&h,TString b){
   hist = (TH1D*)h->Clone("hSignal");
   bcg = b;

   // setup the drawing style
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   fitRangeLow = hist->GetXaxis()->GetXmin();
   fitRangeHigh = hist->GetXaxis()->GetXmax();

   if(DEBUGMODE){
      RooMsgService::instance().setSilentMode(false);
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   }else {
      RooMsgService::instance().setSilentMode(true);
      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   }
}

void Fit::fitPeak() {

   // create a canvas that will hold both fits
   c = new TCanvas("peakFit", "", 1200, 800); 

   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.03);
   gPad->SetBottomMargin(0.11);
   gPad->SetTickx();
   gPad->SetTicky(); 

   gStyle->SetOptStat("");   

   if(!hist || hist->GetEntries() == 0){
      cout << "ERROR in Fit::fitPeak - Empty or no hist when trying to fit. Leaving." << endl;
      return;
   }

   // Set axis title
   hist->GetXaxis()->SetTitle("m [GeV/c^{2}]");
   hist->SetTitle("");


   // Declare observable x
   x = new RooRealVar("x", "m [GeV/c^{2}]", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

   // Create a binned dataset from histogram
   dh = new RooDataHist("dh", "dh", *x, Import(*hist));

   // Create a frame for plotting
   frame = x->frame(Title(""));

   // Plot data points with error bars
   dh->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));

   // Define the components of the model based on the chosen background function
   model = buildModel();

   // Check if model was created successfully before fitting
   if (!model) {
      std::cerr << "ERROR: Fit::buildModel - Model was not initialized. Check option string: " << bcg << std::endl;
      return;
   }
   
   // Define fit range in which the model will be fitted to the data
   if(DEBUGMODE) cout << "Defining fit range..." << endl;
   defineFitRange();

   // Fit model to data
   if(DEBUGMODE) cout << "Starting fit..." << endl;
   fitResult = model->fitTo(*dh, PrintLevel(1), RooFit::Range("fitRange"), RooFit::Minimizer("Minuit2"), Strategy(2)/*, Hesse(true)*/, Save() );

   if (!fitResult) {
      std::cerr << "ERROR: Fit::buildModel - Fit result is null. Check if the model was created successfully." << std::endl;
      return;
   }
   if(DEBUGMODE) cout << "Fit completed." << endl;
   
   // Print fit parameters
   if(DEBUGMODE) fitResult->Print("v");

   // Calculate approximate degrees of freedom
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


   frame->Draw("hist E");
   frame->GetXaxis()->SetTitleSize(0.05);
   frame->GetYaxis()->SetTitleSize(0.05);
   frame->GetXaxis()->SetTitleOffset(0.8);
   frame->GetYaxis()->SetTitleOffset(0.8);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   
   leg1 = new TLegend(legPosX1, legPosY1, legPosX2, legPosY2);
   leg1->SetTextSize(0.03);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(0);
   leg1->AddEntry("data","Data", "LEP");
   
   if(fitSignal && fitBcg) leg1->AddEntry("model","Fit Model", "L"); 
   if(fitSignal) leg1->AddEntry("signal", modelNames[0], "L");
   if(fitBcg) leg1->AddEntry("bkg", modelNames[1], "L");
   
   frame->Draw("hist F E");
   leg1->Draw("same");
   
   if(fitSignal) integrate(bcg);
   
   if(fitSignal)  writeSignalResult(0.7,0.58,0.85,0.78);
   
   if(fitSignal) drawLines();
   

   c->Update();
}

Fit::~Fit(){
   if(c) {
       c->Close();
       delete c;
   }
   if(fitResult) delete fitResult;
   if(model) delete model;
   if(frame) delete frame;
   if(hSignal) delete hSignal;
}

void Fit::defineFitRange(){

   x->setRange("fitRange",fitRangeLow, fitRangeHigh);
}



RooAddPdf* Fit::buildModel(){

   if(DEBUGMODE) cout << "Building model with option: " << bcg << endl;

   models.clear();
   modelNames.clear();
   fitBcg = false;
   // ----------------------------------------------------------------------------------------------
   if(bcg.Contains("cb") || bcg.Contains("CB") || bcg.Contains("CrystalBall") || bcg.Contains("crystalball") ){
      // Define signal model (Crystal Ball)
      cbmean   = new RooRealVar("mean", "mean", mMeanPar, 1.005, 1.05);
      cbsigma  = new RooRealVar("sigma", "sigma", mSigmaPar, 0.00001, 0.5);

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
      modelNames.push_back("Crystal Ball");
      fitSignal = true;
   }else if(bcg.Contains("gauss") || bcg.Contains("Gauss") || bcg.Contains("Gaussian") || bcg.Contains("gaussian")){
      // Define signal model (Gaussian)
      cbmean   = new RooRealVar("mean", "mean", mMeanPar, 1.0, 1.05);
      cbsigma  = new RooRealVar("sigma", "sigma", mSigmaPar, 0.00001, 0.5);
      gauss = new RooGaussian("signal1S", "gauss", *x, *cbmean, *cbsigma);
      models.push_back(gauss);
      modelNames.push_back("Gauss");
      fitSignal = true;
   }else{ modelNames.push_back(""); }
   
   // ----------------------------------------------------------------------------------------------
   if (bcg.Contains("poly1") || bcg.Contains("Poly1"))  {
      a1 = new RooRealVar("a1", "a1", mPolynomialPars[1], -20, 1000);
      // include observable x in variable list so 'x' is recognized
      bkg = new RooPolynomial("bkg", "bkg", *x, RooArgList(*a1));

      models.push_back(bkg);
      modelNames.push_back("Polynomial 1");
      fitBcg = true;
   }else if (bcg.Contains("poly2") || bcg.Contains("Poly2"))  {
      a1 = new RooRealVar("a1", "a1", mPolynomialPars[1], -20, 1000);
      a2 = new RooRealVar("a2", "a2", mPolynomialPars[2], -20, 1000);
      // include observable x so polynomial terms using x are valid
      bkg = new RooPolynomial("bkg", "bkg", *x, RooArgList(*a1, *a2));

      models.push_back(bkg);
      modelNames.push_back("Polynomial 2");
      fitBcg = true;
   } else{ modelNames.push_back(""); }


   if(models.size() == 0 || models.size() > 2){
      cout << "ERROR: Fit::buildModel - Invalid number of models. One can select at most 2 models, one signal (cb, gauss) and background (poly, continuum)." << endl;
      return nullptr;
   }


   RooAddPdf* mod;
   // Combine signal and background into model
   if (fitSignal && fitBcg) {
      nsig = new RooRealVar("nsig", "signal events", 40000, 0, 2000000);
      nbkg = new RooRealVar("nbkg", "background events", 70000, 0, 10000000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0], *models[1]), RooArgList(*nsig, *nbkg));
   } else if(fitSignal){
      nsig = new RooRealVar("nsig", "signal events", 1000, 0, 1500000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0]), RooArgList(*nsig));
   }else if(fitBcg) {   
      nbkg = new RooRealVar("nbkg", "background events", 10, 0, 1500000);
      mod = new RooAddPdf("model", "model", RooArgList(*models[0]), RooArgList(*nbkg));
   }else {
      cerr << "ERROR: Fit::buildModel - No valid model components selected. Please choose at least one signal or background model." << endl;
      return nullptr;
   }

   return mod;
}


void Fit::integrate(TString bcg){
   // Define symmetric integration window
   double low = cbmean->getVal() - 3.0 * cbsigma->getVal();
   double high = cbmean->getVal() + 3.0 * cbsigma->getVal();


   x->setRange("signalRange", low, high);
   if(DEBUGMODE) cout << "Integration range: " << low << " to " << high << endl;

   // Calculate fractions inside the range for signal and background PDFs
   RooAbsReal* intSignal;
   if(cb) {
      intSignal = cb->createIntegral(*x, NormSet(*x), RooFit::Range("signalRange"));
   }else{
      intSignal = gauss->createIntegral(*x, NormSet(*x), RooFit::Range("signalRange"));
   }
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
   
   // Calculate yields in the range (mu -3sigma, mu + 3sigma) by multiplying by the normalization parameters
   yieldSignal = nsig->getVal() * fracSignal;
   yieldBackground = (fracBackground > 0) ? nbkg->getVal() * fracBackground : 0.0;
   
   // Construct formula variables representing yields in the signal range
   RooFormulaVar sigYieldFormula("sigYieldFormula", "@0 * @1", RooArgList(*nsig, *intSignal));
   errSignal = sigYieldFormula.getPropagatedError(*fitResult);
   
   // Calculate net yield and its error
   netYield = yieldSignal;
   errNetYield = errSignal;
   
   // Calculate signal-to-background ratio and its error
   sOverB = (yieldBackground > 0) ? yieldSignal / yieldBackground : 0;
   
   errSOverB = 0;
   if (yieldBackground > 0 && yieldSignal > 0) {
      errSOverB = sOverB * sqrt(
         pow(errSignal / yieldSignal, 2) +
         pow(errBackground / yieldBackground, 2)
      );
   }
   
}

void Fit::saveCanvas(TFile *&file, TString name ,TString dir){
   
   if(!file || file->IsZombie()){
      cerr << "ERROR in Fit::saveCanvas - File is not open or is a zombie. Cannot save canvas." << endl;
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

void Fit::saveCanvas(shared_ptr<TFile> &file, TString name ,TString dir){
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

void Fit::drawLines(){
   
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


void Fit::writeFitResult(){
   leg1->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare()*NDF, NDF, frame->chiSquare() ), "");
}


void Fit::writeSignalResult(double left, double bottom, double right, double top){

   TPaveText *text = new TPaveText(left, bottom, right, top, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
   text->SetTextSize(0.03);
   text->SetFillStyle(0);
   text->SetBorderSize(0);
   text->SetFillColorAlpha(0, 0);
   text->SetTextFont(42);
   text->SetTextAlign(12);
   text->AddText("");
   if(fitSignal){
      text->AddText(Form("#mu = %.4f #pm %.4f",cbmean->getVal(),cbmean->getError()));
      text->AddText(Form("#sigma = %.5f #pm %.5f",cbsigma->getVal(),cbsigma->getError()));
      text->AddText(Form("Raw yield = %.0f #pm %.0f", round(netYield), round( (errNetYield)) ));
   }

   text->Draw("same");
}

