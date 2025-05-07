#include "../include/Plot.h"

Plot::Plot(const string mInputList, const char* filePath){
   outputPosition = filePath;
   inputPosition = mInputList;
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

void Plot::DrawSTARpp510JPsi(double xl, double yl, double xr, double yr, double textSizeRel){
   
   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl, xr, yr,"brNDC");
   textpp510 -> SetTextSize(textSize + textSizeRel);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(4000);
   textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
   textpp510 -> SetTextFont(62);
   textpp510->AddText(ppSTARJPsi[0]);
   textpp510->AddText(ppSTARJPsi[1]);
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


bool Plot::ConnectInputTree(const string& input, TString nameOfTree, bool alsoBcgTree) {
   
   int nInputFiles;
   //cout << "Input from: " << input.c_str() << endl;

   chain = new TChain(nameOfTree);
   if(alsoBcgTree){
      bcgChain = new TChain(nameOfTree + TString("_Bcg"));
   }

   //cout << "Input file: " << input.c_str() << endl;
   if(input.find(".root") != string::npos){
         cout << "Input from root file: "<< input << endl;
      unique_ptr<TFile> inputFile(new TFile(input.c_str(), "read") );
      if(!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()){
         cout<< "Couldn't open input root file..."<<endl;
         return false;
      } 
      chain->AddFile(input.c_str());
      nInputFiles = 1;
      if(alsoBcgTree)
         bcgChain->AddFile(input.c_str());
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
         unique_ptr<TFile> inputFile(new TFile(line.c_str(), "read") ); 
         if(!inputFile){
            cout << "Couldn't open: " << line.c_str() << endl;
            return false;
         } 
         currentTree = dynamic_cast<TTree*>( inputFile->Get(nameOfTree) );
         if(currentTree){
            chain->AddFile(line.c_str());
            if(alsoBcgTree)
               bcgChain->AddFile(line.c_str());
         } else {
            cout << "Name of tree: " << nameOfTree << endl;
            cout << "Couldn't open .root file with name: " << line.c_str() << endl;
         }  
         lineId++;
         nInputFiles++;
         //all histograms should be added together from the input files, saved to .root file, and then the .root file should be opened with pointer inFile. So far i am pointing to the last inputFile, because I use only one. Fuck off i will finish this one day ;)
      }//while

      instr.close();
      tree = dynamic_cast<TTree*>( chain );
      bcgTree = dynamic_cast<TTree*>( bcgChain );

   }//else if

   cout << "Input from " << nInputFiles << " files..." << endl;
   return true;
}


bool Plot::handleHistograms(){

   vector<pair<TH1*,TString>> hists1D;
   vector<pair<TH2*,TString>> hists2D;

   hists1D = GetAllTH1();
   if(hists1D.size() == 0){
      cerr << "Couldn't load 1D histograms from file. Leaving..." << endl;
      return false;
   }
   for (int i = 0; i < hists1D.size(); ++i){
      if(!hists1D[i].first){
         cerr << "Couldn't load histogram " << hists1D[i].second << ". Leaving..." << endl;
         return false;
      }
      TH1General(hists1D[i].second, hists1D[i].first);
      cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
   }

   // create canvases for all TH2F
   hists2D = GetAllTH2();
   if(hists2D.size() == 0){
      cerr << "Couldn't load 2D histograms from file. Leaving..." << endl;
      return false;
   }
   for (int i = 0; i < hists2D.size(); ++i){
      if(!hists2D[i].first){
         cerr << "Couldn't load histogram " << hists2D[i].second << ". Leaving..." << endl;
         return false;
      }
      TH2General(hists2D[i].second, hists2D[i].first);
      cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
   }

   return true;
}



vector<pair<TH2*, TString>> Plot::GetAllTH2() {

   // Vector to store pointers to TH1D histograms
   vector<pair<TH2*,TString>> histograms;

   // Iterate over all keys in the file
   TKey *key;
   TIter next(histFile->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      if (obj->InheritsFrom(TH2D::Class()) || obj->InheritsFrom(TH2F::Class()) || obj->InheritsFrom(TH2I::Class())) {
         // If the object is a TH2F histogram, add it to the vector
         TH2 *h2 = dynamic_cast<TH2*>(obj);
         histograms.push_back(make_pair(h2, h2->GetName()));
         cout << "Just added 2D histogram: " << h2->GetName() << endl;

      }
   }

   return histograms;
}


vector<pair<TH1*, TString>> Plot::GetAllTH1() {

   // Vector to store pointers to TH1D histograms
   vector<pair<TH1*,TString>> histograms;
   // Iterate over all keys in the file
   TKey *key;
   TIter next(histFile->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      //cout << "Name of object: " << obj->GetName() << endl;
      if (obj->InheritsFrom(TH1D::Class()) || obj->InheritsFrom(TH1F::Class()) || obj->InheritsFrom(TH1I::Class())) {
         // If the object is a TH1 histogram, add it to the vector
         TH1 *h1 = dynamic_cast<TH1*>(obj);
         histograms.push_back(make_pair(h1, h1->GetName()));
         cout << "Just added 1D histogram: " << h1->GetName() << endl;
      }
   }
   return histograms;
}



void Plot::TH1General(TString nameOfHist,TH1*& hist) {

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
      hist->GetXaxis()->SetLabelOffset(0.01);
   }else if(strcmp(hist->GetName(), "hNTpcTracks") == 0 || strcmp(hist->GetName(), "hNTofTracks") == 0 || strcmp(hist->GetName(), "hPointingAngle") == 0 ){
      canvas->SetLogy();
   }else if(strcmp(hist->GetName(), "hTrackQualityFlow") == 0){
      canvas->SetLogy();
      SetGPad(false, 0.14,0.1, 0.11,0.06);
   }else{
      SetGPad(false, 0.14, 0.05,0.11,0.06);
   }
   DrawSTARpp510JPsi(0.6,0.85,0.93,0.93, 0.01);

   if(strcmp(hist->GetName(),"hNTracksRP") == 0){
      canvas->SetLogy();
   }

   canvas->Update();

   outFile->cd();
   canvas->Write();
   //hist->Write();
}


void Plot::TH2General(TString nameOfHist , TH2*& hist){
   CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if(nameOfHist.Contains("hNSigma")){
      hist->SetMinimum(0.98);
   }


   canvas->Clear();
   SetTH2Style(hist);

   hist->Draw("COLZ");
   DrawSTARpp510JPsi(0.53,0.85,0.82,0.93, 0.01);
   if(nameOfHist.Contains("hRPcorr") ){
      DrawFiducial();
   }

   outFile->cd();
   canvas->Write();
   //hist->Write();
}

