#include "FitJPsi.h"


FitJPsi::FitJPsi(TH1D *&h,TString b){
   hist = (TH1D*)h->Clone("hSignal");
   bcg = b;

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

   // Define background model
   bool fitBcg = false;

   if (bcg == "poly1" || bcg == "Poly1") {
      fitBcg = true;
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      bkg = new RooPolynomial("bkg", "bkg", *x, RooArgList(*a0));
   } 
   else if (bcg == "poly2" || bcg == "Poly2") {
      fitBcg = true;
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      a1 = new RooRealVar("a1", "a1", -1, -200, 200);
      bkg = new RooPolynomial("bkg", "bkg", *x, RooArgList(*a0, *a1));
   } 
   else if (bcg == "poly3" || bcg == "Poly3") {
      fitBcg = true;
      a0 = new RooRealVar("a0", "a0", 10, 0, 200);
      a1 = new RooRealVar("a1", "a1", -1, -200, 200);
      a2 = new RooRealVar("a2", "a2", 0.5, -200, 200);
      bkg = new RooPolynomial("bkg", "bkg", *x, RooArgList(*a0, *a1, *a2));
   }

   // Define signal model (Crystal Ball)
   cbmean   = new RooRealVar("mean", "mean", 3.0908, 2.8, 3.2);
   cbsigma  = new RooRealVar("sigma", "sigma", 0.0543, 0.0, 0.065);

   if(mAlphaLoose){
      cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41, 0.5, 5.0); // loose
   }else{
      cbalpha  = new RooRealVar("cb_alpha", "cb_alpha", 1.41, 1.410, 1.411); // fixed
   }

   if(mNLoose){
      cbn      = new RooRealVar("cb_n", "cb_n", 1.93, 0.5, 5.0);         // loose
   }else{
      cbn      = new RooRealVar("cb_n", "cb_n", 1.93, 1.930, 1.931);         // fixed
   }

   cb = new RooCBShape("cb", "cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);

   // Combine signal and background into model
   nsig = new RooRealVar("nsig", "signal events", 1500, 0, 15000);
   nbkg = new RooRealVar("nbkg", "background events", 10, 0, 15000);


   if (fitBcg) {
      model = new RooAddPdf("model", "model", RooArgList(*cb, *bkg), RooArgList(*nsig, *nbkg));
   } else {
      model = new RooAddPdf("model", "model", RooArgList(*cb), RooArgList(*nsig));
   }


   // Check if model was created successfully before fitting
   if (!model) {
      std::cerr << "ERROR: Model was not initialized. Check background string: " << bcg << std::endl;
      return;
   }

   // Fit model to data
   fitResult = model->fitTo(*dh, RooFit::Save(), RooFit::PrintLevel(-1));

   if (!fitResult) {
      std::cerr << "ERROR: Fit result is null. Check if the model was created successfully." << std::endl;
      return;
   }

   // Calculate degrees of freedom
   NDF = hist->GetNbinsX() - (fitResult->floatParsFinal().getSize());
   mNLoose ? NDF += 0 : NDF += 1;
   mAlphaLoose ? NDF += 0 : NDF += 1;

   //plotting
   model->plotOn(frame, Name("model"));
   if(fitBcg){
      model->plotOn(frame, Components("bkg"), LineStyle(kDashed), LineColor(kRed), Name("background"));
   }
   //model->plotOn(frame, Components("cb"), LineStyle(kDashed), LineColor(kGreen), Name("signal"));


   frame->Draw("hist E");
   frame->GetXaxis()->SetTitleSize(0.05);
   frame->GetYaxis()->SetTitleSize(0.05);
   frame->GetXaxis()->SetTitleOffset(0.8);
   frame->GetYaxis()->SetTitleOffset(0.8);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   integrate(bcg);

   TPaveText *text = new TPaveText(0.7,0.52,0.85,0.78, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
   text->SetTextSize(0.03);
   text->SetFillStyle(0);
   text->SetBorderSize(0);
   text->SetFillColorAlpha(0, 0);
   text->SetTextFont(42);
   text->SetTextAlign(12);
   text->AddText("");
   text->AddText(Form("#mu = %.4f #pm %.4f",cbmean->getVal(),cbmean->getError()));
   text->AddText(Form("#sigma = %.4f #pm %.4f",cbsigma->getVal(),cbsigma->getError()));
   text->AddText(Form("#chi^{2}/NDF = %.2f/%d #approx %.1f",frame->chiSquare(), NDF, frame->chiSquare()/NDF ) ); 
   text->AddText(Form("Raw yield = %.0f #pm %.0f", round(netYield), round( (errNetYield)) ));
   if(fitBcg)  
      text->AddText(Form("S/B = %.0f / %.0f #approx %.1f", yieldSignal, yieldBackground, sOverB ) );
   text->Draw("same hist E");

   drawLines();

   TLegend *leg1 = new TLegend(0.25,0.58,0.4,0.8);
   leg1->SetTextSize(0.03);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(0);
   leg1->AddEntry("data","Data", "LEP");
   if(!fitBcg){
      leg1->AddEntry("model", "Crystal Ball", "LP");
   }else{
      leg1->AddEntry("model",TString("Crystal Ball + ")+ bcg,"LP");
      leg1->AddEntry("background",bcg, "LP");
   }
   leg1->Draw("same hist E");
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

void FitJPsi::fitContinuum(){

   
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   // create a canvas that will hold both fits
   c = new TCanvas("peakFit", "Fit Result", 1200, 800); 

   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.03);
   gPad->SetTopMargin(0.03);
   gPad->SetBottomMargin(0.11);
   gPad->SetTickx();
   gPad->SetTicky(); 

   gStyle->SetOptStat(0);   

   if(!hist || hist->GetEntries() == 0){
      cout << "Empty or no hist when trying to fit. Leaving." << endl;
      return;
   }

   // Set axis title
   hist->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");
   hist->SetTitle("");
   hist->Draw("EP");

   TF1 *f = new TF1("f", "[0] * exp( [1]*x )", minimumContinuum, maximumContinuum);
   f->SetParameter(0, 300);
   f->SetParameter(1, -0.3);
   f->SetLineColor(kRed);
   f->SetLineWidth(2);
   f->SetLineStyle(1);
   hist->Fit(f, "R");
   

   TPaveText *text = new TPaveText(0.72,0.52,0.92,0.78, "NDC"); //in plot text (x_beggining, y_beggining, x_end, y_end) .
   text->SetTextSize(0.03);
   text->SetFillStyle(0);
   text->SetBorderSize(0);
   text->SetFillColorAlpha(0, 0);
   text->SetTextFont(42);
   text->SetTextAlign(12);
   text->AddText("f(x) = A #times exp( #lambda x )");
   text->AddText(Form("    A = %.0f #pm %.0f",f->GetParameter(0),f->GetParError(0)));
   text->AddText(Form("    #lambda = %.2f #pm %.2f",f->GetParameter(1),f->GetParError(1)));
   text->AddText(Form("chi^{2}/NDF = %.2f/%d #approx %.1f", f->GetChisquare(), f->GetNDF() ,  f->GetChisquare()/f->GetNDF() ) ); 
   text->Draw("same");

   TLegend *leg1 = new TLegend(0.35,0.75,0.5,0.85);
   leg1->SetTextSize(0.03);
   leg1->SetFillStyle(0);
   leg1->SetBorderSize(0);
   leg1->AddEntry(hist,"Data", "LEP");
   leg1->AddEntry(f, "Exponential Fit", "LP");
   leg1->Draw("same");
   c->Update();


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