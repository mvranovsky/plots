#include "FitJPsi.h"


FitJPsi::FitJPsi(TH1D *&h,TString b){
   hist = (TH1D*)h->Clone("hSignal");
   bcg = b;

   fitRangeLow = hist->GetXaxis()->GetXmin();
   fitRangeHigh = hist->GetXaxis()->GetXmax();

   //RooMsgService::instance().setSilentMode(true);
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
      cout << "Empty or no hist when trying to fit. Leaving." << endl;
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
   bool fitBcg = false, fitSignal = false, fitSignal2S = false;
   model = buildModel(fitBcg, fitSignal, fitSignal2S);


   // Check if model was created successfully before fitting
   if (!model) {
      std::cerr << "ERROR: FitJPsi::buildModel - Model was not initialized. Check option string: " << bcg << std::endl;
      return;
   }

   defineFitRange();

   // Fit model to data
   fitResult = model->fitTo(*dh, PrintLevel(1), RooFit::Range("fitRange"), RooFit::Minimizer("Minuit2"), Strategy(2)/*, Hesse(true)*/, Save() );
   //fitResult = model->fitTo(*dh, Save(), PrintLevel(-1),RooFit::Verbose(false), RooFit::Warnings(false), RooFit::Range("fitRange"));

   cout << "Fit completed." << endl;
   if (!fitResult) {
      std::cerr << "ERROR: FitJPsi::buildModel - Fit result is null. Check if the model was created successfully." << std::endl;
      return;
   }

   // Calculate degrees of freedom
   NDF = hist->GetNbinsX() - (fitResult->floatParsFinal().getSize());
   mNLoose ? NDF += 0 : NDF += 1;
   mAlphaLoose ? NDF += 0 : NDF += 1;

   //plotting
   model->plotOn(frame, Name("model"));
   if(fitBcg) model->plotOn(frame, Components("bkg"), LineStyle(kDashed), LineColor(kRed), Name("bkg"));
   if(fitSignal) model->plotOn(frame, Components("signal1S"), LineStyle(kDashed), LineColor(kGreen), Name("signal"));
   if(fitSignal2S) model->plotOn(frame, Components("signal2S"), LineStyle(kDashed), LineColor(kMagenta), Name("signal2"));


   frame->Draw("hist E");
   frame->GetXaxis()->SetTitleSize(0.05);
   frame->GetYaxis()->SetTitleSize(0.05);
   frame->GetXaxis()->SetTitleOffset(0.8);
   frame->GetYaxis()->SetTitleOffset(0.8);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   cout << "After fit plotting." << endl;
   
   if(fitSignal) integrate(bcg);
   
   if(fitSignal)  writeSignalResult(fitSignal2S,0.7,0.58,0.85,0.78);
   
   if(fitSignal) drawLines();
   
   leg1 = new TLegend(legPosX1, legPosY1, legPosX2, legPosY2);
   leg1->SetTextSize(0.03);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(0);
   leg1->AddEntry("data","Data", "LEP");
   leg1->AddEntry("model","Model Fit", "L");
   if(fitSignal) leg1->AddEntry("signal", modelNames[0], "L");
   if(fitBcg) leg1->AddEntry("bkg", modelNames[1], "L");
   if(fitSignal2S) leg1->AddEntry("signal2", modelNames[2], "L");
   //leg1->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), NDF, frame->chiSquare()/NDF ), "");
   leg1->Draw("same");
   

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

   cout << "Building model with option: " << bcg << endl;

   models.clear();
   modelNames.clear();
   fitBcg = false;
   // ----------------------------------------------------------------------------------------------
   if(bcg.Contains("cb") || bcg.Contains("CB") || bcg.Contains("CrystalBall") || bcg.Contains("crystalball") ){
      // Define signal model (Crystal Ball)
      cbmean   = new RooRealVar("mean", "mean", 3.0908, 2.8, 3.2);
      cbsigma  = new RooRealVar("sigma", "sigma", 0.0543, 0.0, 0.065);

      if(mAlphaLoose){  // for fitting alpha
         cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41, 0.5, 5.0); // loose
      }else{            // for choosing specific alpha
         cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41); // fixed
         cbalpha->setConstant(true);
      }
      if(mNLoose){      // for fitting n
         cbn      = new RooRealVar("cb_n", "cb_n", 1.93, 0.5, 5.0);         // loose
      }else{         // for choosing specific n 
         cbn      = new RooRealVar("cb_n", "cb_n", 1.93);
         cbn->setConstant(true);
      }

      cb = new RooCBShape("signal1S", "cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
      models.push_back(cb);
      modelNames.push_back("Crystal Ball ( J/#psi )");
      fitSignal = true;
   }else if(bcg.Contains("gauss") || bcg.Contains("Gauss") || bcg.Contains("Gaussian") || bcg.Contains("gaussian")){
      // Define signal model (Gaussian)
      cbmean   = new RooRealVar("mean", "mean", 3.0908, 2.8, 3.2);
      cbsigma  = new RooRealVar("sigma", "sigma", 0.0543, 0.0, 0.065);
      gauss = new RooGaussian("signal1S", "gauss", *x, *cbmean, *cbsigma);
      models.push_back(gauss);
      modelNames.push_back("Gauss ( J/#psi )");
      fitSignal = true;
   }else{ modelNames.push_back(""); }
   // ----------------------------------------------------------------------------------------------
   if (bcg.Contains("poly1") || bcg.Contains("Poly1")) {
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      // include observable x in the variable list so expression 'a0' resolves in RooGenericPdf
      bkg = new RooGenericPdf("bkg", "bkg", "a0", RooArgList(*x, *a0));
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
      c2 = new RooRealVar("c2", "c2 parameter", 0.7, -3.0, 3.0); //-0.01
      c3 = new RooRealVar("c3", "c3 parameter",-0.7, -3.0, 1.0); //-0.45
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
      modelNames.push_back("Crystal Ball (#psi (2S))");
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
      nbkg = new RooRealVar("nbkg", "background events", 700, 0, 1000000);
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


double FitJPsi::getCorrectedYield(TH1D *spectrum, TGraphAsymmErrors* graph, double averageCorrection) { //spectrum should be background subtracted, because so is inv mass


   if(netYield <= 0 || errNetYield <= 0) {
      cerr << "Error: The fit has not been run or failed. Before getting corrected yield, please run the fit to get the raw yield." << endl;
      return -1;
   }

   
   double frac = sOverB <= 0 ? 1 : yieldSignal/(yieldSignal + yieldBackground);

   
   correctedYield = 0;
   correctedYieldErrTop = 0;
   correctedYieldErrLow = 0;
   for(int i = 1; i <= spectrum->GetNbinsX() + 1 ; ++i){
      
      if(i == (spectrum->GetNbinsX() + 1)){  //embedding does not reach beyond 1.5 GeV and neither does spectrum, but there are overflow bins, just use the last value
         //cout << "Bin " << i << ": spectrum content = " << spectrum->GetBinContent(i) << ", graph y = " << graph->GetY()[graph->GetN()-1] << ", frac = " << frac << endl;
         correctedYield += spectrum->GetBinContent(i)*frac/graph->GetY()[graph->GetN()-1];
         correctedYieldErrTop += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYhigh()[graph->GetN()-1]/graph->GetY()[graph->GetN()-1], 2);
         correctedYieldErrLow += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYlow()[graph->GetN()-1]/graph->GetY()[graph->GetN()-1], 2);
         break;
      }

      double binCenter = spectrum->GetBinCenter(i);
      for(int j = 0; j < graph->GetN(); j++){

         
         if(binCenter > (graph->GetX()[j]- graph->GetEXlow()[j]) && binCenter < (graph->GetX()[j]+ graph->GetEXhigh()[j])){
            //cout << "Bin " << i << ": spectrum content = " << spectrum->GetBinContent(i) << ", graph y = " << graph->GetY()[j] << ", frac = " << frac << endl;
            correctedYield += spectrum->GetBinContent(i)*frac/graph->GetY()[j];
            correctedYieldErrTop += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYhigh()[j]/graph->GetY()[j], 2);
            correctedYieldErrLow += pow(spectrum->GetBinError(i)/spectrum->GetBinContent(i), 2) + pow(graph->GetEYlow()[j]/graph->GetY()[j], 2);
            break;
         }
      }
   }
   correctedYield = round(correctedYield);
   correctedYieldErrTop = sqrt(correctedYieldErrTop);
   correctedYieldErrLow = sqrt(correctedYieldErrLow);
   
   /*
   cout << "The raw yield from fitting: " << yieldSignal << endl;
   cout << "The corrected yield has been computed. Result is " << correctedYield << endl;
   cout << "The corrected yield globally: " << yieldSignal << "/" << averageCorrection << " = " << yieldSignal/averageCorrection << endl;
   */
   
   TPaveText *text = new TPaveText(0.7,0.52,0.85,0.49, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
   text->SetTextSize(0.03);
   text->SetFillColor(0);
   text->SetTextFont(42);
   text->SetTextAlign(12);
   text->AddText(Form("Corrected yield = %.0f^{+%.0f}_{-%.0f}", correctedYield, correctedYieldErrTop, correctedYieldErrLow));
   text->Draw("same");
   
   
   return correctedYield;
}


void FitJPsi::integrate(TString bcg){
   // 1) Define asymmetric integration window
   double low = cbmean->getVal() - 4.0 * cbsigma->getVal();
   double high = cbmean->getVal() + 3.0 * cbsigma->getVal();


   x->setRange("signalRange", low, high);
   cout << "Integration range: " << low << " to " << high << endl;

   // 2) Calculate fractions inside the range for signal and background PDFs
   RooAbsReal* intSignal = cb->createIntegral(*x, NormSet(*x), Range("signalRange"));
   double fracBackground = 0.0;
   RooAbsReal* intBackground = nullptr;
   if(bcg == "poly1" || bcg == "poly2" || bcg == "poly3" || bcg == "Poly1" || bcg == "Poly2" || bcg == "Poly3"){
      intBackground = bkg->createIntegral(*x, NormSet(*x), Range("signalRange"));
      fracBackground = intBackground->getVal();
      RooFormulaVar bkgYieldFormula("bkgYieldFormula", "@0 * @1", RooArgList(*nbkg, *intBackground));
      errBackground = bkgYieldFormula.getPropagatedError(*fitResult);
   }

   double fracSignal = intSignal->getVal();
   cout << "Signal fraction in range: " << fracSignal << endl;
   cout << "nsig: " << nsig->getVal() << endl;
   
   // 3) Calculate yields in the range by multiplying by the normalization parameters
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
   cText->AddText(Form("c_{1} = %.4f #pm %.4f", c1->getVal(), c1->getError()));
   cText->AddText(Form("c_{2} = %.4f #pm %.4f", c2->getVal(), c2->getError()));
   cText->AddText(Form("c_{3} = %.4f #pm %.4f", c3->getVal(), c3->getError()));
   //cText->AddText(Form("c_{4} = %.4f #pm %.4f", c4->getVal(), c4->getError()));
   cText->Draw("same");

}

void FitJPsi::writeFitResult(){
   leg1->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), NDF, frame->chiSquare()/NDF ), "");
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
   //text->AddText(Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), NDF, frame->chiSquare()/NDF ) );
   text->AddText(Form("#mu_{J/#psi} = %.4f #pm %.4f",cbmean->getVal(),cbmean->getError()));
   text->AddText(Form("#sigma_{J/#psi} = %.4f #pm %.4f",cbsigma->getVal(),cbsigma->getError()));
   text->AddText(Form("Raw yield (J/#psi )= %.0f #pm %.0f", round(netYield), round( (errNetYield)) ));
   if(fitBcg)  text->AddText(Form("S/B (J/#psi ) = %.0f / %.0f #approx %.1f", yieldSignal, yieldBackground, sOverB ) );
   if(fitSignal2S){
      text->AddText(Form("#mu_{#psi(2S)} = %.4f #pm %.4f",mean2S->getVal(),mean2S->getError()));
      text->AddText(Form("#sigma_{#psi(2S)} = %.4f #pm %.4f",sigma2S->getVal(),sigma2S->getError()));
      text->AddText(Form("Raw yield (#psi (2S))= %.0f #pm %.0f", round(nsig2S->getVal()), round(nsig2S->getError()) ));
   }
   text->Draw("same");
}


void FitJPsi::addContinuumFunction(double c_1, double c_2, double c_3){
   // Define the custom PDF. Include observable 'x' in the RooArgList so expressions referencing x compile.
   RooRealVar *c1_const = new RooRealVar("c1_const", "c1 parameter", c_1);
   RooRealVar *c2_const = new RooRealVar("c2_const", "c2 parameter", c_2);
   RooRealVar *c3_const = new RooRealVar("c3_const", "c3 parameter", c_3);
   RooRealVar *c4_const = new RooRealVar("c4_const", "c4 parameter",0.00001, 0.0, 0.0001); //-0.45

      // Define the custom PDF. Include observable 'x' in the RooArgList so expressions referencing x compile.
   RooGenericPdf *continuumFunc = new RooGenericPdf("continuumFunc", "Custom function", "sqrt((x - c1_const)^2 + c4_const) * exp( c2_const *(x - c1_const)^2 + c3_const*x^3)", RooArgList(*x, *c1_const, *c2_const, *c3_const, *c4_const));
   RooAddPdf *modelWithContinuum = new RooAddPdf("modelWithContinuum", "model with continuum", RooArgList(*continuumFunc), RooArgList(*nbkg));
   
   x->setRange("fitRange",c_1, fitRangeHigh);
   modelWithContinuum->plotOn(frame, RooFit::LineColor(kMagenta), Name("modelWithCont"), RooFit::Range("fitRange"));
   leg1->AddEntry("modelWithCont", "Fit from data (scaled)", "L");

   frame->Draw("hist E");
   leg1->Draw("same");

}