#include "../include/Plot.h"

Plot::Plot(TFile *mOutFile, const string mInputList, const char* filePath){
   inFile = mOutFile;
   outputPosition = filePath;
   inputPosition = mInputList;
   //nameOfTree = treeName;
   //plots = plotsFromManager;
   //mUtil = new Util();
}

Plot::~Plot(){
   //if(mUtil) delete mUtil;
}


void Plot::CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth, int canvasHeight)
{
   *canvas = new TCanvas(canvasName,canvasName,canvasWidth,canvasHeight);
   gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   //gPad->SetLogy(0);

   // setup the drawing style
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetNumberContours(99);
   gStyle->SetPalette(55);
}


void Plot::SetHistStyle(TH1* hist, Int_t color, Int_t markStyle)
{
   hist->GetYaxis()->SetTitle(yAxisTitle);
   hist->SetStats(false);
   hist->GetXaxis()->SetTitleFont(textFont);
   hist->GetXaxis()->SetTitleFont(textFont);
   hist->GetXaxis()->SetLabelFont(textFont);
   hist->GetYaxis()->SetLabelFont(textFont);
   hist->GetXaxis()->SetLabelSize(labelSize);
   hist->GetYaxis()->SetLabelSize(labelSize);
   hist->GetXaxis()->SetTitleSize(labelSize);
   hist->GetYaxis()->SetTitleSize(labelSize);
   hist->GetXaxis()->SetTitleOffset(1.);
   hist->GetYaxis()->SetTitleOffset(1.4);
   //hist->GetYaxis()->SetRangeUser(0, 0.4); 
   hist->SetLineColor(color);
   hist->SetLineStyle(1);
   hist->SetLineWidth(1);  
   hist->SetMarkerSize(2);
   hist->SetMarkerColor(color);
   hist->SetMarkerStyle(markerStyle);
}//SetHistStyle

void Plot::SetTH2Style(TH2* hist)
{
   hist->SetStats(false);
   hist->GetXaxis()->SetTitleFont(textFont);
   hist->GetYaxis()->SetTitleFont(textFont);
   hist->GetZaxis()->SetTitleFont(textFont);
   hist->GetXaxis()->SetLabelFont(textFont);
   hist->GetYaxis()->SetLabelFont(textFont);
   hist->GetZaxis()->SetLabelFont(textFont);
   hist->GetXaxis()->SetLabelSize(labelSize);
   hist->GetYaxis()->SetLabelSize(labelSize);
   hist->GetZaxis()->SetLabelSize(labelSize);
   hist->GetXaxis()->SetTitleSize(labelSize);
   hist->GetYaxis()->SetTitleSize(labelSize);
   hist->GetZaxis()->SetTitleSize(labelSize);
   hist->GetXaxis()->SetTitleOffset(0.8);
   hist->GetYaxis()->SetTitleOffset(0.9);
   hist->GetZaxis()->SetTitleOffset(0.9);
}//SetTH2Style

void Plot::CreateLegend(TLegend **legend, double xl, double yl, double xr, double yr)
{
   *legend = new TLegend(xl, yl,xr, yr);
   (*legend)->SetFillStyle(0);
   (*legend)->SetBorderSize(0);
   (*legend)->SetTextAlign(11);
   (*legend)->SetTextSize(textSize);
   (*legend)->SetTextFont(42);
   (*legend)->SetMargin(0.1);   
}

void Plot::DrawSTARInternal(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize);
   textSTAR -> SetFillColor(0);
   textSTAR -> SetTextFont(62);
   textSTAR->AddText("STAR Internal");
   textSTAR -> Draw("same");
}

void Plot::DrawSTARpp510(double xl, double yl, double xr, double yr, double textSizeRel){
   
   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl, xr, yr,"brNDC");
   textpp510 -> SetTextSize(textSize + textSizeRel);
   textpp510 -> SetTextAlign(11);
   textpp510 -> SetFillColor(0);
   textpp510 -> SetTextFont(62);
   textpp510->AddText(ppSTAR);
   textpp510 -> Draw("same");
}


void Plot::CreateText(TString writtenText, double xl, double yl, double xr, double yr)
{
   text = new TPaveText(xl, yl, xr, yr,"brNDC");
   text -> SetTextSize(textSize);
   text -> SetFillColor(0);
   text -> SetTextFont(textFont);
   text -> SetTextAlign(11);   
   text -> AddText(writtenText);
   text -> Draw("same");
}

 
void Plot::SetGPad(bool isLogY, double left, double right, double bottom, double top)
{
   //gPad->SetMargin(xl,yl,xr,yr); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetLeftMargin(left);
   gPad->SetRightMargin(right);
   gPad->SetTopMargin(top);
   gPad->SetBottomMargin(bottom);
   gPad->SetTickx();
   gPad->SetTicky(); 
   if (isLogY)
      gPad->SetLogy();
   gStyle->SetOptStat("");   
}

void Plot::DrawFiducial()
{
   const Int_t n = 100;
   Double_t x[n], y[n];
   Double_t tmp;
   for(int i = 0; i < n; ++i){
      x[i] = 0.175 + (0.27*i)/n;
      //tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
      tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
      y[i] = -sqrt(abs(tmp));
   }
   TGraph* gr = new TGraph(n,x,y);
   gr->SetLineWidth(4);
   gr->Draw("same");


   TLine *left02 = new TLine(-0.27,-0.4,0.445,-0.4);
   SetLineStyle(left02);
   left02->Draw("same");


   TLine *left01 = new TLine(-0.27,-0.8,-0.27,-0.4);
   SetLineStyle(left01);
   left01->Draw("same");

   left01 = new TLine(-0.27,-0.8,0.185,-0.8);
   SetLineStyle(left01);
   left01->Draw("same");          
   // UP
   left02 = new TLine(-0.27,0.4,0.445,0.4);
   SetLineStyle(left02);
   left02->Draw("same");


   left01 = new TLine(-0.27,0.4,-0.27,0.8);
   SetLineStyle(left01);
   left01->Draw("same");

   left01 = new TLine(-0.27,0.8,0.185,0.8);
   SetLineStyle(left01);
   left01->Draw("same");

   for(int i = 0; i < n; ++i){
      x[i] = 0.175 + (0.27*i)/n;
      //tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
      tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
      y[i] = sqrt(abs(tmp));
   }
   gr = new TGraph(n,x,y);
   gr->SetLineWidth(4);
   gr->Draw("same"); 
}

void Plot::SetLineStyle(TLine* line)
{
   line->SetLineStyle(1);
   line->SetLineColor(1);
   line->SetLineWidth(4);
}//SetLineStyle

int Plot::fitGaussPol2(TH1D **histToFit, Int_t binWidth, Double_t minRange, Double_t maxRange, Double_t pol0, Double_t pol1, Double_t pol2, Double_t amp,Double_t mean, Double_t sigma) {
   //func is function of signal(gauss) and background(pol2)
   func = new TF1("func","pol2(0) + gaus(3)",minRange,maxRange);
   func->SetParNames("p_{0}", "p_{1}", "p_{2}", "A","#mu","#sigma");
   func->SetParameter(0, pol0);
   func->SetParameter(1, pol1);
   func->SetParameter(2, pol2);
   func->SetParameter(3, amp); //amplitude //Lambda
   func->SetParameter(4, mean); //mean
   func->SetParameter(5, sigma); //sigma
   func->SetLineColor(kBlue);
   func->SetLineWidth(2);
   (*histToFit)->Fit("func", "", "", minRange, maxRange);   //Lambda   
   func->Draw("same");

   Double_t upperLim, lowerLim;
   pol0 = func->GetParameter(0);
   pol1 = func->GetParameter(1);
   pol2 = func->GetParameter(2);
   amp = func->GetParameter(3);
   mean = func->GetParameter(4);
   sigma = func->GetParameter(5);
   upperLim = mean + 3*sigma;
   lowerLim = mean - 3*sigma;

   bcg = new TF1("bcg","pol2", minRange, maxRange); //Lambda
   bcg->SetParameter(0,pol0);
   bcg->SetParameter(1,pol1);
   bcg->SetParameter(2,pol2);
   bcg->SetLineColor(kRed);
   bcg->SetLineWidth(2);
   bcg->Draw("same");


   Double_t i_s, i_b;
   int y;
   i_s = func->Integral(lowerLim, upperLim);  //under peak of pol+gaus
   i_b = bcg->Integral(lowerLim, upperLim);  //under peak of pol from pol+gaus
   cout << "i_s = " << i_s << "\n";    
   cout << "i_b = " << i_b << "\n";    
   y = (i_s - i_b) / binWidth;             //yield calculation must be renormalised to bin size in GeV (here bin is 50 MeV)

   return y;
}




bool Plot::defineAnalysis(){
   //define which analysis ran and what histograms are needed to be plotted

   if (runAnaV0){
      nameOfTree =nameOfAnaV0Tree;
      plots = plotsV0;
      return true;
   } else if (runAnaBP){
      nameOfTree = nameOfAnaBPTree; 
      plots = plotsBP;
      return true;
   } else if (runAnaV0Control){
      nameOfTree = nameOfAnaV0ControlTree;
      plots = plotsV0Control;
      return true;  
   } else if(runAnaV0SingleState){
      nameOfTree = nameOfAnaV0SingleStateTree;
      plots = plotsV0;
      return true;  
   } else {
      cout << "All analyses set to false. Returning..." << endl;
      return false;
   }
}//defineAnalysis




bool Plot::ConnectInputTree(const string& input) {
   
   int nInputFiles;
   //cout << "Input from: " << input.c_str() << endl;

   TFile *inputFile;
   chain = new TChain(nameOfTree);

   //cout << "Input file: " << input.c_str() << endl;
   if(input.find(".root") != string::npos){
         cout << "Input from root file: "<< input << endl;
      inputFile = TFile::Open(input.c_str(), "read");
      if(!inputFile){
         cout<< "Couldn't open input root file..."<<endl;
         return false;
      } 
      chain->AddFile(input.c_str());
      nInputFiles = 1;
   } 
   else if(input.find(".list") != string::npos ){
      cout << "Input from root list: " << input.c_str() << endl;
      ifstream instr(input.c_str());
      if (!instr.is_open()){
         cout<< "Couldn't open: "<<input.c_str()<<endl;
         return false;
      }
      int lineId = 0;
      string line;
      TTree *currentTree;
       while(getline(instr, line)) {
         if(line.empty())
            continue;
         inputFile = TFile::Open(line.c_str(), "read");
         if(!inputFile){
            cout << "Couldn't open: " << line.c_str() << endl;
            return false;
         } 
         currentTree = dynamic_cast<TTree*>( inputFile->Get(nameOfTree) );
         if(currentTree){
            chain->AddFile(line.c_str());
         } else {
            cout << "Name of tree: " << nameOfTree << endl;
            cout << "Couldn't open .root file with name: " << line.c_str() << endl;
         }  
         lineId++;
         nInputFiles++;

      }//while

      instr.close();
      tree = dynamic_cast<TTree*>( chain );
   }//else if

   inputFile->Close();
   //might take a lot of memory
   
   cout << "Input from " << nInputFiles << " files..." << endl;
   return true;
}




vector<pair<TH2F*, TString>> Plot::GetAllTH2F() {
    // Open the ROOT file
   if (!inFile || inFile->IsZombie()) {
      cerr << "Error opening file1 to get TH2F." << endl;
      return {};
   }

   // Vector to store pointers to TH1D histograms
   vector<pair<TH2F*,TString>> histograms;

   // Iterate over all keys in the file
   TKey *key;
   TIter next(inFile->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      if (TH2F *h1 = dynamic_cast<TH2F*>(obj)) {
         // If the object is a TH2F histogram, add it to the vector
         histograms.push_back(make_pair(h1, h1->GetName()));
      }
   }

   return histograms;
}


vector<pair<TH1D*, TString>> Plot::GetAllTH1D() {
    // Open the ROOT file
   if (!inFile || inFile->IsZombie()) {
      cerr << "Error opening inFile to get TH1D." << endl;
      return {};
   }

   // Vector to store pointers to TH1D histograms
   vector<pair<TH1D*,TString>> histograms;

   // Iterate over all keys in the file
   TKey *key;
   TIter next(inFile->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      if (TH1D *h1 = dynamic_cast<TH1D*>(obj)) {
         // If the object is a TH2F histogram, add it to the vector
         histograms.push_back(make_pair(h1, h1->GetName()));
      }
   }

   return histograms;
}



void Plot::TH1DGeneral(TString nameOfHist,TH1D* hist) {

   if (!hist){
      cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }
   
   if (nameOfHist = "") {
      CreateCanvas(&canvas, hist->GetName(), widthTypical, heightTypical );      
   }else{
      CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
   }

   canvas->Clear();
   SetHistStyle(hist, kBlack, markerStyleTypical);
   hist->Draw("same E hist");
   //canvas->SetLogy();


   if(strcmp(hist->GetName(), "hAnalysisFlow") == 0){
      canvas->SetLogy();
      SetGPad(false, 0.14, 0.05,0.11,0.06);
      DrawSTARpp510(0.6,0.85,0.9,0.85, 0.02);

     hist->GetXaxis()->SetBinLabel(3, TString("pre-selection"));
     hist->GetXaxis()->SetBinLabel(4, TString("Good quality tracks"));
     hist->GetXaxis()->SetBinLabel(7, TString("|V_{Z}| < 60 cm"));
     hist->GetXaxis()->SetBinLabel(8, TString("Q_{tot} = 0"));


   }else{
      SetGPad(false, 0.14, 0.05,0.11,0.06);
      DrawSTARpp510(0.6,0.85,0.9,0.85, 0.02);
   }

   canvas->Update();

   outFile->cd();
   canvas->Write();
   //hist->Write();
}


void Plot::TH2FGeneral(TString nameOfHist , TH2F* hist){
   CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }
   canvas->Clear();
   SetTH2Style(hist);
   hist->Draw("COLZ");
   DrawSTARpp510();

   outFile->cd();
   canvas->Write();
   //hist->Write();
}

void Plot::Clear(){
   if(inFile){
      inFile->Close();
   }

   if(outFile){
      outFile->Close();
   }
}