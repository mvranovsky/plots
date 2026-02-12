#include "Plot.h"

Plot::Plot(const std::string& mInputList, const char* filePath)
    : outputPosition(filePath), inputPosition(mInputList), outFile(std::make_shared<TFile>(outputPosition.c_str(), "RECREATE"))
{
}

Plot::Plot(const std::string& mInputList, const std::shared_ptr<TFile>& file)
    : outputPosition(file ? file->GetName() : ""), inputPosition(mInputList),outFile(file)
{
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

/*
void Plot::SetTH1Style(TH1*& hist, Int_t color, Int_t markStyle)
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
}//SetTH1Style


void Plot::SetTH1Style(TH1D*& hist, Int_t color, Int_t markStyle)
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
}//SetTH1Style

void Plot::SetTH1Style(TH1F*& hist, Int_t color, Int_t markStyle)
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
}//SetTH1Style

void Plot::SetTH1Style(TH1I*& hist, Int_t color, Int_t markStyle)
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
}//SetTH1Style
*/
void Plot::SetTH2Style(TH2*& hist)
{

// Define color gradient (blue → green → yellow)
   const Int_t NRGBs = 3;
   const Int_t NCont = 255;
   Double_t stops[NRGBs]  = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]    = { 0.00, 0.00, 1.00 };
   Double_t green[NRGBs]  = { 0.00, 1.00, 1.00 };
   Double_t blue[NRGBs]   = { 1.00, 1.00, 0.00 };

   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);

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


void Plot::SetTH2Style(TH2D*& hist)
{

// Define color gradient (blue → green → yellow)
   const Int_t NRGBs = 3;
   const Int_t NCont = 255;
   Double_t stops[NRGBs]  = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]    = { 0.00, 0.00, 1.00 };
   Double_t green[NRGBs]  = { 0.00, 1.00, 1.00 };
   Double_t blue[NRGBs]   = { 1.00, 1.00, 0.00 };

   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);

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



void Plot::SetTH2Style(TH2F*& hist)
{

// Define color gradient (blue → green → yellow)
   const Int_t NRGBs = 3;
   const Int_t NCont = 255;
   Double_t stops[NRGBs]  = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]    = { 0.00, 0.00, 1.00 };
   Double_t green[NRGBs]  = { 0.00, 1.00, 1.00 };
   Double_t blue[NRGBs]   = { 1.00, 1.00, 0.00 };

   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);

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

void Plot::SetTH2Style(TH2I*& hist)
{

// Define color gradient (blue → green → yellow)
   const Int_t NRGBs = 3;
   const Int_t NCont = 255;
   Double_t stops[NRGBs]  = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]    = { 0.00, 0.00, 1.00 };
   Double_t green[NRGBs]  = { 0.00, 1.00, 1.00 };
   Double_t blue[NRGBs]   = { 1.00, 1.00, 0.00 };

   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);

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

void Plot::SetTGraphStyle(TGraph*& graph, Int_t color, Int_t markStyle)
{
   graph->SetLineColor(color);
   graph->SetLineWidth(2);
   graph->SetMarkerColor(color);
   graph->SetMarkerStyle(markStyle);
   graph->SetMarkerSize(1.5);
   graph->GetXaxis()->SetTitleFont(textFont);
   graph->GetYaxis()->SetTitleFont(textFont);
   graph->GetXaxis()->SetLabelFont(textFont);
   graph->GetYaxis()->SetLabelFont(textFont);
   graph->GetXaxis()->SetLabelSize(labelSize);
   graph->GetYaxis()->SetLabelSize(labelSize);
   graph->GetXaxis()->SetTitleSize(labelSize);
   graph->GetYaxis()->SetTitleSize(labelSize);
}
void Plot::SetTGraphStyle(TGraphErrors*& graph, Int_t color, Int_t markStyle)
{
   graph->SetLineColor(color);
   graph->SetLineWidth(2);
   graph->SetMarkerColor(color);
   graph->SetMarkerStyle(markStyle);
   graph->SetMarkerSize(1.5);
   graph->GetXaxis()->SetTitleFont(textFont);
   graph->GetYaxis()->SetTitleFont(textFont);
   graph->GetXaxis()->SetLabelFont(textFont);
   graph->GetYaxis()->SetLabelFont(textFont);
   graph->GetXaxis()->SetLabelSize(labelSize);
   graph->GetYaxis()->SetLabelSize(labelSize);
   graph->GetXaxis()->SetTitleSize(labelSize);
   graph->GetYaxis()->SetTitleSize(labelSize);
}

void Plot::SetTGraphStyle(TGraphAsymmErrors*& graph, Int_t color, Int_t markStyle)
{
   graph->SetLineColor(color);
   graph->SetLineWidth(2);
   graph->SetMarkerColor(color);
   graph->SetMarkerStyle(markStyle);
   graph->SetMarkerSize(1.5);
   graph->GetXaxis()->SetTitleFont(textFont);
   graph->GetYaxis()->SetTitleFont(textFont);
   graph->GetXaxis()->SetLabelFont(textFont);
   graph->GetYaxis()->SetLabelFont(textFont);
   graph->GetXaxis()->SetLabelSize(labelSize);
   graph->GetYaxis()->SetLabelSize(labelSize);
   graph->GetXaxis()->SetTitleSize(labelSize);
   graph->GetYaxis()->SetTitleSize(labelSize);
}

void Plot:: CreateLegend(TLegend **legend, double xl, double yl, double xr, double yr)
{
   *legend = new TLegend(xl, yl,xr, yr);
   (*legend)->SetFillStyle(0);
   (*legend)->SetBorderSize(0);
   (*legend)->SetTextAlign(11);
   (*legend)->SetTextSize(textSize);
   (*legend)->SetTextFont(42);
   (*legend)->SetMargin(0.1);   
}


void Plot::DrawSTARTag(double xl, double yl, double xr, double yr, double textSizeRel){

   TPaveText *text1;
   text1 = new TPaveText(xl, yl, xr, yr,"brNDC");
   text1 -> SetTextSize(textSize+0.01);
   text1 -> SetBorderSize(0);
   text1 -> SetFillStyle(0);
   text1 -> SetFillColorAlpha(kWhite, 0.0);
   text1 -> SetTextFont(72);
   text1 -> SetTextAlign(33);
   text1->AddText(plotDescription);
   text1 -> Draw("same");
 
   TPaveText *text2;
   text2 = new TPaveText(xl, yl-0.05, xr, yr-0.05,"brNDC");
   text2 -> SetTextSize(textSize );
   text2 -> SetTextAlign(33);
   text2 -> SetFillStyle(0);
   text2 -> SetFillColorAlpha(kWhite, 0.0);  // white background with 0% opacity
   text2 -> SetBorderSize(0);
   text2 -> SetTextFont(62);
   text2 -> AddText(process);
   text2 -> AddText(energyOfCollision);
   text2 -> Draw("same");
}



void Plot::CreateText(TString writtenText, int textF,double xl, double yl, double xr, double yr)
{
   text = new TPaveText(xl, yl, xr, yr,"brNDC");
   text -> SetTextSize(textSize);
   text -> SetFillColorAlpha(0,0);
   text -> SetTextFont(textF);
   text -> SetTextAlign(11);   
   text -> AddText(writtenText);
   text -> Draw("same");
}

void Plot::SetGPad(bool isLogY, double left, double right, double bottom, double top) // set standard margins
{
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


void Plot::SetLineStyle(TLine*& line)
{
   line->SetLineStyle(1);
   line->SetLineColor(1);
   line->SetLineWidth(3);
}//SetLineStyle

void Plot::CreateLine(double xl, double yl, double xr, double yr)
{
   line = new TLine(xl, yl, xr, yr);
   line->SetLineStyle( lineStyle );
   line->SetLineColor( lineColor );
   line->SetLineWidth( lineWidth );
}//CreateLine

void Plot::CreateDashedLine(double xl, double yl, double xr, double yr)
{
   l = new TLine(xl, yl, xr, yr);
   l->SetLineStyle(2);
   l->SetLineColor( kBlack );
   l->SetLineWidth( 3 );
   l->Draw("same");
}//CreateLine


bool Plot::ConnectInputTree(const string& input, TString nameOfTree, TTree *&TREE) {

   chain = new TChain(nameOfTree);   

   vector<TFile*> inputFiles;
   //cout << "Input file: " << input.c_str() << endl;
   if(input.find(".root") != string::npos){   //!!!! something is wrong
         if(DEBUGMODE) cout << "Input from root file: "<< input << endl;
      TFile *inputFile = new TFile(input.c_str(), "read");
      if(!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()){
         cerr << "ERROR in Plot::ConnectInputTree: Couldn't open input root file..."<<endl;
         return false;
      } 

      inputFiles.push_back(inputFile);  // Store the file so it doesn't get deleted

      TTree* currentTree = dynamic_cast<TTree*>(inputFile->Get(nameOfTree));
      if (currentTree) {
         chain->AddFile(input.c_str());
      }
   } 
   else if(input.find(".list") != string::npos ){
      if(DEBUGMODE) cout << "Input from root list: " << input.c_str() << endl;
      ifstream instr(input.c_str());
      if (!instr.is_open()){
         cerr << "ERROR in Plot::ConnectInputTree: Couldn't open input root list file..."<<endl;
         return false;
      }
      int lineId = 0;
      string line;
      TTree *currentTree;
       while(getline(instr, line)) {
         if(line.empty())
            continue;
         if(line.find(".root") == string::npos)
            continue;
         
         TFile *inputFile = new TFile(line.c_str(), "read");
         if(!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()){
            cerr << "ERROR in Plot::ConnectInputTree: Couldn't open file: " << line.c_str() << endl;
            continue;
         }
         inputFiles.push_back(inputFile);  // Store the file so it doesn't get deleted

         currentTree = dynamic_cast<TTree*>( inputFile->Get(nameOfTree) );
         if(currentTree){
            chain->AddFile(line.c_str());
         } else {
            if(DEBUGMODE) cout << "Name of tree: " << nameOfTree << endl;
            if(DEBUGMODE) cout << "Couldn't open .root file with name: " << line.c_str() << endl;
         }  
         lineId++;
         //all histograms should be added together from the input files, saved to .root file, and then the .root file should be opened with pointer inFile. So far i am pointing to the last inputFile, because I use only one. Fuck off i will finish this one day ;)
      }//while

      instr.close();
      
   }//else if
   
   TREE = dynamic_cast<TTree*>( chain );

   if(DEBUGMODE) cout << "Input from " << inputFiles.size() << " files..." << endl;
   return true;
}


bool Plot::saveHistograms(TString dir){

   vector<pair<TH1*,TString>> hists1D;
   vector<pair<TH2*,TString>> hists2D;


   hists1D = GetAllTH1();
   if(hists1D.size() == 0){
      cerr << "ERROR in Plot::saveHistograms: Couldn't load 1D histograms from file. Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists1D.size(); ++i){
      if(!hists1D[i].first){
         cerr << "ERROR in Plot::saveHistograms: Couldn't load histogram " << hists1D[i].second << " in dir "<< dir << ". Leaving..." << endl;
         return false;
      }
      TH1General( hists1D[i].first,hists1D[i].second, dir);
      if(DEBUGMODE) cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
   }

   
   // create canvases for all TH2F
   hists2D = GetAllTH2();
   if(hists2D.size() == 0){
      cerr << "ERROR in Plot::saveHistograms: Couldn't load 2D histograms from file in dir "<< dir << ". Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists2D.size(); ++i){
      if(!hists2D[i].first){
         cerr << "ERROR in Plot::saveHistograms: Couldn't load histogram " << hists2D[i].second << " in dir "<< dir << ". Leaving..." << endl;
         return false;
      }
      TH2General( hists2D[i].first,hists2D[i].second, dir);
      if(DEBUGMODE) cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
   }
   
   
   return true;
}




vector<pair<TH2*, TString>> Plot::GetAllTH2() {
   vector<pair<TH2*, TString>> histograms;
   /*
   // Navigate to the directory in the ROOT file
   TDirectory* targetDir = histFile.get(); // root file is a TDirectory
   if (!dir.IsNull()) {
      targetDir = dynamic_cast<TDirectory*>(histFile->Get(dir));
      if (!targetDir) {
         cerr << "ERROR in Plot::GetAllTH2: Directory not found: " << dir << endl;
         return histograms;
      }
   }
   */

    map<TString, int> saved;
    // Iterate over all keys in the directory
    TKey *key;
    TIter next(histFile->GetListOfKeys());

    while ((key = static_cast<TKey*>(next()))) {
        TObject *obj = key->ReadObj();
        if (!obj) continue;

        if (obj->InheritsFrom(TH2::Class())) {
            TH2 *h2 = dynamic_cast<TH2*>(obj);
            if (h2 && h2->GetEntries() > 0) {
               if(saved.find(h2->GetName()) == saved.end()){
                  saved[h2->GetName()] = 1;
                  //h2->Write();
               }else{
                  if(DEBUGMODE) cout << "Histogram " << h2->GetName() << " already saved. Skipping." << endl;
               }
               histograms.push_back(make_pair(h2, h2->GetName()));
               if(DEBUGMODE) cout << "Added 2D histogram: " << h2->GetName() << endl;
            }
        }

        // Clean up if ReadObj creates new objects
        // Do not delete h2 here; ownership is typically with the directory unless cloned
    }

    return histograms;
}


vector<pair<TH1*, TString>> Plot::GetAllTH1() {
    vector<pair<TH1*, TString>> histograms;
   /*
   // Navigate to the directory in the ROOT file
   TDirectory* targetDir = histFile.get(); // root file is a TDirectory
   if (!dir.IsNull()) {
      targetDir = dynamic_cast<TDirectory*>(histFile->Get(dir));
      if (!targetDir) {
         cerr << "ERROR in Plot::GetAllTH1(): Directory not found: " << dir << endl;
         return histograms;
      }
   }else{
      cerr << "ERROR in Plot::GetAllTH1(): Directory name is empty. Please provide a valid directory name. Leaving..." << endl;
      return histograms;
   }
   */

    // Iterate over all keys in the directory
    TKey *key;
    TIter next(histFile->GetListOfKeys());

    while ((key = static_cast<TKey*>(next()))) {
        TObject *obj = key->ReadObj();
        if (!obj) continue;

        if (obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom(TH2::Class())) { // Ensure it's a 1D histogram
            TH1 *h1 = dynamic_cast<TH1*>(obj);
            if (h1 && h1->GetEntries() > 0) {
               histograms.push_back(make_pair(h1, h1->GetName()));
               if(DEBUGMODE) cout << "Added 1D histogram: " << h1->GetName() << endl;
            }
        }

        // Clean up if ReadObj creates new objects
        // Do not delete h2 here; ownership is typically with the directory unless cloned
    }

    return histograms;
}


void Plot::TH1General(TH1*& hist,TString nameOfHist,TString dir ) {
   
   if (!hist){
      cerr << "ERROR in Plot::TH1General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName(), widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
   }
   canvas->Clear();
   
   SetGPad(false, 0.14, 0.05,0.11,0.06);
   
   SetTH1Style(hist, kBlue, markerStyle);
   
   DrawSTARTag();

   canvas->Update();

   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
}

void Plot::TH1General(TH1D*& hist,TString nameOfHist,TString dir ) {
   
   if (!hist){
      cerr << "ERROR in Plot::TH1General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName(), widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
   }
   canvas->Clear();
   
   SetGPad(false, 0.14, 0.05,0.11,0.06);
   
   SetTH1Style(hist, kBlue, markerStyle);
   
   DrawSTARTag();

   canvas->Update();

   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
}

void Plot::TH1General(TH1F*& hist,TString nameOfHist,TString dir ) {
   
   if (!hist){
      cerr << "ERROR in Plot::TH1General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName() , widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist , widthTypical, heightTypical );
   }
   canvas->Clear();
   
   SetGPad(false, 0.14, 0.05,0.11,0.06);
   
   SetTH1Style(hist, kBlue, markerStyle);
   
   DrawSTARTag();

   canvas->Update();

   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
}

void Plot::TH1General(TH1I*& hist,TString nameOfHist,TString dir ) {
   
   if (!hist){
      cerr << "ERROR in Plot::TH1General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName() , widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist , widthTypical, heightTypical );
   }
   canvas->Clear();
   
   SetGPad(false, 0.14, 0.05,0.11,0.06);
   
   SetTH1Style(hist, kBlue, markerStyle);
   
   DrawSTARTag();

   canvas->Update();

   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
}


void Plot::TH2General(TH2*& hist, TString nameOfHist , TString dir){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, hist->GetName(), widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist, widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "ERROR in Plot::TH2General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }


   SetTH2Style(hist);

   DrawSTARTag();

   c->Update();
   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   c->Write();
   //hist->Write();
   c->Close();
   c->Clear();
   delete c;
}


void Plot::TH2General(TH2I*& hist, TString nameOfHist , TString dir){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, hist->GetName() , widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist , widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "ERROR in Plot::TH2General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }


   SetTH2Style(hist);

   DrawSTARTag();

   c->Update();
   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   c->Write();
   //hist->Write();
   c->Close();
   c->Clear();
   delete c;
}

void Plot::TH2General(TH2D*& hist, TString nameOfHist , TString dir){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, hist->GetName() , widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist , widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "ERROR in Plot::TH2General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }


   SetTH2Style(hist);

   DrawSTARTag();

   c->Update();
   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   c->Write();
   //hist->Write();
   c->Close();
   c->Clear();
   delete c;
}

void Plot::TH2General(TH2F*& hist, TString nameOfHist , TString dir){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, hist->GetName() , widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist, widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "ERROR in Plot::TH2General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }


   SetTH2Style(hist);

   DrawSTARTag();

   c->Update();
   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   c->Write();
   //hist->Write();
   c->Close();
   c->Clear();
   delete c;
}


void Plot::loadMCTree( TString filename, TString treename ){


   if(DEBUGMODE) cout << "Connecting MC file..." << endl;


   MCFile = shared_ptr<TFile>(new TFile(filename, "READ"));
   if(!MCFile || MCFile->IsZombie() || !MCFile->IsOpen()){
      cerr << "ERROR in Plot::loadMCTree(): Could not get file with MC data. Returning." << endl;
      return;
   }
   
   if(DEBUGMODE) cout << "Checking MC tree..." << endl;

   MCTree = (TTree*)MCFile->Get(treename);
   if(!MCTree || MCTree->IsZombie() || MCTree->GetEntries() <= 0){
      cerr << "ERROR in Plot::loadMCTree(): Could not get MC tree. Returning." << endl;
      return;
   }

   if(DEBUGMODE) cout << "Successfully loaded MC tree" << endl;

}
