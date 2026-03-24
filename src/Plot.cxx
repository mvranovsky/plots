#include "Plot.h"

Plot::Plot(const string mInputList, const char* filePath){
   outputPosition = filePath;
   inputPosition = mInputList;
   outFile = shared_ptr<TFile>(new TFile(outputPosition, "RECREATE"));
}

Plot::Plot(const string mInputList, shared_ptr<TFile> file) : outFile(file), inputPosition(mInputList) {
   outputPosition = file->GetName();
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


void Plot::SetHistStyle(TH1*& hist, Int_t color, Int_t markStyle)
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

void Plot::DrawSTARInternal(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize+0.01);
   textSTAR -> SetFillStyle(0);
   textSTAR -> SetFillColorAlpha(kWhite,0);
   textSTAR -> SetBorderSize(0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl - 0.05, xr, yr - 0.05,"brNDC");
   textpp510 -> SetTextSize(textSize);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite,0);
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(62);
   textpp510 -> SetTextAlign(33);
   textpp510->AddText("p+p #sqrt{s} = 510 GeV");
   textpp510 -> Draw("same");
}

void Plot::DrawSTARInternalZB(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize+0.01);
   textSTAR -> SetFillStyle(0);
   textSTAR -> SetFillColorAlpha(kWhite,0);
   textSTAR -> SetBorderSize(0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl - 0.07, xr, yr - 0.05,"brNDC");
   textpp510 -> SetTextSize(textSize);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite,0);
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(62);
   textpp510 -> SetTextAlign(33);
   textpp510->AddText("RP_zerobias data");
   textpp510->AddText("p+p #sqrt{s} = 510 GeV");
   textpp510 -> Draw("same");
}


void Plot::DrawSTARpp510(double xl, double yl, double xr, double yr, double textSizeRel){
   
   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl, xr, yr,"brNDC");
   textpp510 -> SetTextSize(textSize + textSizeRel);
   textpp510 -> SetTextAlign(11);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite,0);
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(62);
   textpp510->AddText(ppSTAR);
   textpp510 -> Draw("same");
}

void Plot::DrawSTARpp510JPsi(double xl, double yl, double xr, double yr, double textSizeRel){

   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize+0.01);
   textSTAR -> SetBorderSize(0);
   textSTAR -> SetFillStyle(0);
   textSTAR -> SetFillColorAlpha(kWhite, 0.0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");
 
   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl-0.05, xr, yr-0.05,"brNDC");
   textpp510 -> SetTextSize(textSize );
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // white background with 0% opacity
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(62);
   if(noRomanPots){
      textpp510->AddText(ppSTARJPsi[2]);
      textpp510->AddText(ppSTARJPsi[3]);
   }else{
      textpp510->AddText(ppSTARJPsi[0]);
      textpp510->AddText(ppSTARJPsi[1]);
   }
   textpp510 -> Draw("same");
}

void Plot::DrawEmbeddingpp510JPsi(double xl, double yl, double xr, double yr, double textSizeRel){
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize+0.01);
   textSTAR -> SetBorderSize(0);
   textSTAR -> SetFillStyle(0);
   textSTAR -> SetFillColorAlpha(kWhite, 0.0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl-0.05, xr, yr-0.05,"brNDC");
   textpp510 -> SetTextSize(textSize + textSizeRel);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite, 0.0);
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(42);
   textpp510->AddText("Zerobias Embedding");
   textpp510->AddText("MC: p + p #rightarrow J/#psi");
   textpp510 -> Draw("same");
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


   left02 = new TLine(-0.27,-0.4,0.445,-0.4);
   SetLineStyle(left02);
   left02->Draw("same");


   left01 = new TLine(-0.27,-0.8,-0.27,-0.4);
   SetLineStyle(left01);
   left01->Draw("same");

   left03 = new TLine(-0.27,-0.8,0.185,-0.8);
   SetLineStyle(left03);
   left03->Draw("same");          
   // UP
   left04 = new TLine(-0.27,0.4,0.445,0.4);
   SetLineStyle(left04);
   left04->Draw("same");


   left05 = new TLine(-0.27,0.4,-0.27,0.8);
   SetLineStyle(left05);
   left05->Draw("same");

   left06 = new TLine(-0.27,0.8,0.185,0.8);
   SetLineStyle(left06);
   left06->Draw("same");

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

void Plot::DrawElipse(double centerX, double centerY, double radius1, double radius2)
{
   circle = new TEllipse(centerX, centerY, radius1, radius2);
   // 3. Set visual style (optional)
   circle->SetFillStyle(0);   // Important: 0 means no fill (transparent)
   circle->SetLineColor(2);   // Set the line color (e.g., 2 is red)
   circle->SetLineWidth(3);   // Set the line thickness 
   circle->Draw("SAME"); // Use "SAME" to draw it on the same canvas as the histogram
}

void Plot::DrawEtaVtxZFiducial(){

   l1 = new TLine(-0.9, 100, 0.5, 100);
   SetLineStyle(l1);
   l1->Draw("same");

   l2 = new TLine(0.5, 100, 0.9, 0);
   SetLineStyle(l2);
   l2->Draw("same");

   l3 = new TLine(0.9, 0, 0.9, -100);
   SetLineStyle(l3);
   l3->Draw("same");

   l4 = new TLine(0.9, -100, -0.5, -100);
   SetLineStyle(l4);
   l4->Draw("same");

   l5 = new TLine(-0.5, -100, -0.9, 0);
   SetLineStyle(l5);
   l5->Draw("same");

   l6 = new TLine(-0.9, 0, -0.9, 100);
   SetLineStyle(l6);
   l6->Draw("same");

}

bool Plot::ConnectInputTree(const string& input, TString nameOfTree, TTree *&TREE, TTree *&BCGTREE) {

   chain = new TChain(nameOfTree);

   bcgChain = new TChain(nameOfTree + TString("_Bcg"));
   

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
         bcgChain->AddFile(input.c_str());
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
            bcgChain->AddFile(line.c_str());
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
   BCGTREE = dynamic_cast<TTree*>( bcgChain );

   if(!BCGTREE || BCGTREE->IsZombie() || BCGTREE->GetEntries() == 0) {
      if(DEBUGMODE) cerr << "WARNING in Plot::ConnectInputTree: Did not load background tree. Continuing without it..." << endl;
   }

   if(DEBUGMODE) cout << "Input from " << inputFiles.size() << " files..." << endl;
   return true;
}


bool Plot::handleHistograms(TString dir, TString addOn){

   vector<pair<TH1*,TString>> hists1D;
   vector<pair<TH2*,TString>> hists2D;


   hists1D = GetAllTH1(dir);
   if(hists1D.size() == 0){
      cerr << "ERROR in Plot::handleHistograms: Couldn't load 1D histograms from file in dir "<< dir << ". Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists1D.size(); ++i){
      if(!hists1D[i].first){
         cerr << "ERROR in Plot::handleHistograms: Couldn't load histogram " << hists1D[i].second << " in dir "<< dir << ". Leaving..." << endl;
         return false;
      }
      TH1General(hists1D[i].second, hists1D[i].first, dir, addOn);
      if(DEBUGMODE) cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
   }

   
   // create canvases for all TH2F
   hists2D = GetAllTH2(dir);
   if(hists2D.size() == 0){
      cerr << "ERROR in Plot::handleHistograms: Couldn't load 2D histograms from file in dir "<< dir << ". Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists2D.size(); ++i){
      if(!hists2D[i].first){
         cerr << "ERROR in Plot::handleHistograms: Couldn't load histogram " << hists2D[i].second << " in dir "<< dir << ". Leaving..." << endl;
         return false;
      }
      TH2General(hists2D[i].second, hists2D[i].first, dir, addOn);
      if(DEBUGMODE) cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
   }
   
   
   return true;
}




vector<pair<TH2*, TString>> Plot::GetAllTH2(TString dir) {
    vector<pair<TH2*, TString>> histograms;

    // Navigate to the directory in the ROOT file
    TDirectory* targetDir = histFile.get(); // root file is a TDirectory
    if (!dir.IsNull()) {
        targetDir = dynamic_cast<TDirectory*>(histFile->Get(dir));
        if (!targetDir) {
            cerr << "ERROR in Plot::GetAllTH2: Directory not found: " << dir << endl;
            return histograms;
        }
    }

    map<TString, int> saved;
    // Iterate over all keys in the directory
    TKey *key;
    TIter next(targetDir->GetListOfKeys());

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


vector<pair<TH1*, TString>> Plot::GetAllTH1(TString dir) {
    vector<pair<TH1*, TString>> histograms;

    // Navigate to the directory in the ROOT file
    TDirectory* targetDir = histFile.get(); // root file is a TDirectory
    if (!dir.IsNull()) {
        targetDir = dynamic_cast<TDirectory*>(histFile->Get(dir));
        if (!targetDir) {
            cerr << "ERROR in Plot::GetAllTH1(): Directory not found: " << dir << endl;
            return histograms;
        }
    }

    // Iterate over all keys in the directory
    TKey *key;
    TIter next(targetDir->GetListOfKeys());

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


void Plot::TH1General(TString nameOfHist,TH1*& hist,TString dir, TString addOn ) {
   
   if (!hist){
      cerr << "ERROR in Plot::TH1General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName() + addOn, widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist  + addOn, widthTypical, heightTypical );
   }
   canvas->Clear();
   
   SetGPad(false, 0.14, 0.05,0.11,0.06);
   
   SetHistStyle(hist, kBlue, markerStyleTypical);
   TH1Specifics(nameOfHist, hist, addOn );
   
   DrawSTARInternal();

   canvas->Update();

   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
}

void Plot::TH1Specifics(TString nameOfHist, TH1*& hist, TString addOn){ // function called for every TH1, set specifics for each plot here


   if(strcmp(hist->GetName(), "hAnalysisFlow") == 0 ) {
      canvas->SetLogy(1);
      SetGPad(false, 0.14, 0.05,0.11,0.06);
      hist->GetXaxis()->SetLabelOffset(0.01);
      hist->SetMarkerSize(1.2);      
      hist->SetMarkerColor(kBlack);
      hist->SetLineColor(kBlack);
      gStyle->SetPaintTextFormat("4.2f"); // Controls text format (e.g., 4 wide, 2 decimals)
      if(addOn.Contains("NoRP")){
         hist->GetYaxis()->SetRangeUser( 1e+3, 3e+10); 
      }else if(addOn.Contains("WithRP")){
         hist->GetYaxis()->SetRangeUser( 300, 3e+10);
      }else if(addOn.Contains("ZBData")){
         hist->GetYaxis()->SetRangeUser( 1e+2, 3e+10);
         hist->GetXaxis()->SetBinLabel(2, "ZB trigger");
         hist->GetXaxis()->SetBinLabel(3, "BEMC tracks");
         hist->GetXaxis()->SetBinLabel(7, "J/#psi trigger");
      }else if(addOn.Contains("BEData")){
         gPad->SetBottomMargin(0.14);
         hist->GetYaxis()->SetRangeUser( 1e+4, 3e+10);
      }
      hist->Draw("hist TEXT");
      return;
   }

   if(strcmp(hist->GetName(),"hGoodRunListFlow") == 0){
      SetGPad(false, 0.14, 0.05,0.11,0.06);
      hist->GetXaxis()->SetLabelOffset(0.01);
      hist->SetMarkerSize(1.2);      
      hist->SetMarkerColor(kBlack);
      hist->SetLineColor(kBlack);
      gStyle->SetPaintTextFormat("4.2f"); // Controls text format (e.g., 4 wide, 2 decimals)
      hist->GetYaxis()->SetRangeUser( 1400, 3000);
      hist->Draw("hist TEXT");
      return;
   }

   if(strcmp(hist->GetName(), "hEClusters") == 0){
      hist->Rebin(2);
   }
   
   hist->Draw("EP");
   TH1D* hClone = (TH1D*)hist->Clone("hClone");
   hClone->Draw("same hist"); 

   
   if(strcmp(hist->GetName(), "hNTpcTracks") == 0 || strcmp(hist->GetName(), "hNTofTracks") == 0 || strcmp(hist->GetName(), "hPointingAngle") == 0 ){
      canvas->SetLogy(1);
   }
   /*
   if(strcmp(hist->GetName(), "hTrackQualityFlow") == 0){
      canvas->SetLogy(1);
      SetGPad(false, 0.14,0.08, 0.13,0.06);
      hist->GetXaxis()->SetLabelOffset(1.1);
      DrawSTARInternal( 0.75, 0.89, 0.89,0.93);
      
   }
   */
   
   if(strcmp(hist->GetName(), "hDeltaDipAngle") == 0 ){
      hist->GetXaxis()->SetNdivisions(505);
      hist->GetYaxis()->SetRangeUser(hist->GetMinimum()*0.9, hist->GetMaximum()*1.2);
      hist->GetXaxis()->SetBinLabel(6, "V_{Z} - #eta Cut");
   }

   
   if(strcmp(hist->GetName(),"hNTracksRP") == 0){
      canvas->SetLogy(1);
   }

   if(nameOfHist == "hPosZ"){  // z position of vertex
      int maxVal  = hist->GetMaximum();
      CreateDashedLine(-vertexZInCm[0], maxVal*0.5, -vertexZInCm[0], 0);
      CreateDashedLine(vertexZInCm[0], maxVal*0.5, vertexZInCm[0], 0);
   }
   
   if(strcmp(hist->GetName(),"hDcaZ") == 0){ // dca z
      canvas->SetLogy(1);
      int maxVal  = hist->GetMaximum();
      CreateDashedLine(-dcaZInCm[0], maxVal*0.7, -dcaZInCm[0], 0);
      CreateDashedLine(dcaZInCm[0], maxVal*0.7, dcaZInCm[0], 0);
   }


   if(strcmp(hist->GetName(),"hDcaXY") == 0){ // dca xy
      //canvas->SetLogy(1);
      int maxVal  = hist->GetMaximum();
      CreateDashedLine(dcaXYInCm[0], maxVal*0.7, dcaXYInCm[0], 0);
   }

      
   if(strcmp(hist->GetName(),"hNfitHits") == 0){ // dca z
      int maxVal  = hist->GetMaximum();
      CreateDashedLine(nHitsFit[0], maxVal*0.7, nHitsFit[0], 0);
   }
      
   if(strcmp(hist->GetName(),"hNhitsDEdx") == 0){ // dca z
      int maxVal  = hist->GetMaximum();
      CreateDashedLine(nHitsDedx[0], maxVal*0.7, nHitsDedx[0], 0);
   }

   if(nameOfHist.Contains("hDelta")){
      int maxVal  = hist->GetMaximum();
      hist->GetYaxis()->SetRangeUser(0, maxVal*1.3);
   }


   if(strcmp(hist->GetName(),"hDeltaPhiAbs") == 0){
      int maxVal  = hist->GetMaximum();

      CreateDashedLine(TMath::Pi() - TMath::Pi()/6, maxVal*0.7, TMath::Pi() - TMath::Pi()/6 , 0);
      CreateDashedLine(TMath::Pi() + TMath::Pi()/6, maxVal*0.7, TMath::Pi() + TMath::Pi()/6 , 0);

   }   

   if(strcmp(hist->GetName(),"hPIDChiee") == 0){
      int maxVal  = hist->GetMaximum();
      hist->GetYaxis()->SetRangeUser(0, maxVal*1.3);
      CreateDashedLine( chiSquareE[0], maxVal*0.7, chiSquareE[0], 0);
   }

   
   if(strcmp(hist->GetName(),"hBranchRP") == 0){
      hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
   }


}

void Plot::TH2Specifics(TString nameOfHist, TH2*& hist, TString addOn){ // function called for every TH2, set specifics for each plot here
   
   if(nameOfHist.Contains("hEClusterTrack")){
      gPad->SetLogz(1);
   }


   hist->Draw("COLZ");

   if(nameOfHist.Contains("hRPcorr") ){
      DrawFiducial();
   }
   
   if(nameOfHist.Contains("hNSigma")){
      gPad->SetLogz(1);
      hist->SetMinimum(0.98);
   }

   // circle in chi square
   if(nameOfHist.Contains("hNSigmaEE")){
      DrawElipse(0,0, 3, 3);
      gPad->SetLogz(1);
      hist->SetMinimum(0.98);
   }


   //eta vtx z
   if(nameOfHist.Contains("hEtaVtxZ")){
      //gPad->SetLogz(1);
      hist->SetMinimum(0.98);
      DrawEtaVtxZFiducial();
   }

   /*
   if(strcmp(hist->GetName(),"") == 0){
      
   }

   */
   
}

void Plot::TH2General(TString nameOfHist , TH2*& hist, TString dir, TString addOn){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, hist->GetName() + addOn, widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist + addOn, widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "ERROR in Plot::TH2General: Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }


   SetTH2Style(hist);

   TH2Specifics(nameOfHist, hist, addOn );

   //DrawSTARpp510JPsi(0.53,0.85,0.82,0.93, 0.01);
   DrawSTARInternal(0.7, 0.89, 0.83, 0.93);

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


TString Plot::convertToString(double val) {

   ostringstream streamA;
   streamA << fixed << setprecision(1) << val;
   TString formattedA = streamA.str();

   return formattedA;
}

int Plot::makeInt(double val) {
	int result = val;
	return result;
}


void Plot::nSigmaCorrPlot(int particles, TString condition) {
   // 1 = electron, 2 = pion, 3 = kaon, 4 = proton


   if(particles == 1){
      tree->Draw("nSigmaTPCelectronMinus:nSigmaTPCelectronPlus>>histnsig(40,-4,4,40,-4,4)", condition);
      tree->Draw("chiSquareelectron>>chiPlot(50,0,50)", condition );
   }else if(particles == 2){
      tree->Draw("nSigmaTPCpionPlus:nSigmaTPCpionMinus>>histnsig(40,-4,4,40,-4,4)", condition );
      tree->Draw("chiSquarekion>>chiPlot(50,0,50)", condition );
   }else if(particles == 3){
      tree->Draw("nSigmaTPCkaonPlus:nSigmaTPCkaonMinus>>histnsig(40,-4,4,40,-4,4)", condition );
      tree->Draw("chiSquarekaon>>chiPlot(50,0,50)", condition );
   }else if(particles == 4){
      tree->Draw("nSigmaTPCprotonPlus:nSigmaTPCprotonMinus>>histnsig(40,-4,4,40,-4,4)", condition );
      tree->Draw("chiSquareproton>>chiPlot(50,0,50)", condition );
   }else{
      cout << "Unknown particle type. Returning." << endl;
      return;
   }
   TH2* nSigma = dynamic_cast<TH2*>(gDirectory->FindObject("histnsig"));
   if(nSigma && nSigma->GetEntries() > 0){
      nSigma->GetXaxis()->SetTitle("n#sigma_{e} [-]");
      nSigma->GetYaxis()->SetTitle("n#sigma_{e} [-]");
      gStyle->SetOptStat("euo"); // Enable underflow and overflow bins
      TH2General(TString("nSigmaJPsi"), nSigma);
      cout << "Obtained and saved nSigma histogram." << endl;
   }else{
      cout << "Could not get nSigma histogram." << endl;
   }
   
   TH1* chiSquare = dynamic_cast<TH1*>(gDirectory->FindObject("chiPlot"));
   if(chiSquare && chiSquare->GetEntries() > 0){
      chiSquare->GetXaxis()->SetTitle("#chi^{2}_{ee} [-]");
      chiSquare->GetYaxis()->SetTitle("counts");
      TH1General(TString("chiSquareJPsi"), chiSquare);
      cout << "Obtained and saved chiSquare histogram." << endl;
   }else{
      cout << "Could not get chiSquare histogram." << endl;
   }
}



TH1D* Plot::loadInvMassHist(int numBins, Double_t minRange, Double_t maxRange, TString c, bool removeBcg){ 

   TRandom3 randGen(0); // seed with 0 for different results each time
   int seed = randGen.Integer(100);
   int seed_c = randGen.Integer(100);

   // create a canvas that will hold both fits
   CreateCanvas(&canvas, "invMassJPsi" + TString::Format("_%d", seed_c), widthTypical, heightTypical);
   SetGPad();

   TString nameSignal = "invMassJPsi" + TString::Format("_%d", seed);;
   TString nameBcg = "invMassJPsiBcg" + TString::Format("_%d", seed);;
   TH1D *hSignal = new TH1D(nameSignal, "invMassJPsi", numBins, minRange, maxRange);
   TH1D *hBcg = new TH1D(nameBcg, "invMassJPsiBcg", numBins, minRange, maxRange);

   hSignal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");

   TString cmd;
   tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), c);
   hSignal->Add((TH1D*)(gPad->GetPrimitive("hist")));
   if(removeBcg){
      bcgTree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), c);
      hBcg->Add((TH1D*)(gPad->GetPrimitive("hist")));
   }

   if(!hSignal || hSignal->GetEntries() == 0 ){
      cout << "ERROR in Plot::loadInvMassHist - Could not load signal to histograms" << endl;
      return nullptr;
   }


   if(hBcg->GetEntries() == 0 ){
      return hSignal;
   }else{
      hSignal->Add(hBcg, -1); // subtract background from signal
      return hSignal;
   }
    
}


TString Plot::getCondition(TString var, int j ){

   if(var.Contains("tight") || var.Contains("TIGHT") || var.Contains("Tight") ){
      j = 1;
   }else if(var.Contains("loose") || var.Contains("LOOSE") || var.Contains("Loose")){
      j = 2;
   }

   TString c = "";
   if(var.Contains("nHitsFit") || var.Contains("NHITSFIT") || var.Contains("nhitsfit")){
      c += Form("nHitsFit0 > %d && nHitsFit1 > %d && ", nHitsFit[j], nHitsFit[j]);
   }else{
      c += Form("nHitsFit0 > %d && nHitsFit1 > %d && ", nHitsFit[0], nHitsFit[0]);
   }

   if(var.Contains("nHitsDedx") || var.Contains("NHITSDEDX") || var.Contains("nhitsdedx")){
      c += Form("nHitsDEdx0 > %d && nHitsDEdx1 > %d && ", nHitsDedx[j], nHitsDedx[j]);
   }else{
      c += Form("nHitsDEdx0 > %d && nHitsDEdx1 > %d && ", nHitsDedx[0], nHitsDedx[0]);
   }

   if(var.Contains("etaHadron") || var.Contains("etahadron") || var.Contains("ETAHADRON") || var.Contains("eta")){
      c += Form("abs(etaHadron0) < %.1f && abs(etaHadron1) < %.1f && ", etaHadron[j],etaHadron[j]);
   }else{
      c += Form("abs(etaHadron0) < %.1f && abs(etaHadron1) < %.1f && ", etaHadron[0], etaHadron[0]);
   }

   if(var.Contains("chiSquare") || var.Contains("PID") || var.Contains("chiSquareEE") || var.Contains("chiSquareE")){
      c += Form("chiSquareelectron < %.d && ", chiSquareE[j]);
   }else{
      c += Form("chiSquareelectron < %.d && ", chiSquareE[0]);
   }
   
   if(!var.Contains("nonsigma")){

      if(var.Contains("nSigmaProton") || var.Contains("NSIGMAPROTON") || var.Contains("nsigmaproton")){
         c+= Form("abs(nSigmaTPCprotonPlus) > %.1f && abs(nSigmaTPCprotonMinus) > %.1f && ", nSigmaP[j], nSigmaP[j]);
      }else{
         c+= Form("abs(nSigmaTPCprotonPlus) > %.1f && abs(nSigmaTPCprotonMinus) > %.1f && ", nSigmaP[0], nSigmaP[0]);
         
      }
      
      if(var.Contains("nSigmaKaon") || var.Contains("NSIGMAKAON") || var.Contains("nsigmakaon")){
         c+= Form("abs(nSigmaTPCkaonPlus) > %.1f && abs(nSigmaTPCkaonMinus) > %.1f && ", nSigmaK[j], nSigmaK[j]);
      }else{
         c+= Form("abs(nSigmaTPCkaonPlus) > %.1f && abs(nSigmaTPCkaonMinus) > %.1f && ", nSigmaK[0], nSigmaK[0]);
      }
   }
      

   if( !(var.Contains("embedding") || var.Contains("Embedding") || var.Contains("EMBEDDING")) ){  //embedding is not able to reconstruct vertex 
      
      //c += Form("chiSquarepion > %.d && chiSquareproton > %.d && chiSquarekaon > %.d && ", chiSquarePi[0], chiSquareP[0], chiSquareK[0]);
      
      if(var.Contains("dcaZInCm") || var.Contains("dcazincm") || var.Contains("DCAZINCM")){
         c += Form("abs(dcaZInCm0) < %.1f && abs(dcaZInCm1) < %.1f && ", dcaZInCm[j], dcaZInCm[j]);
      }else{
         c += Form("abs(dcaZInCm0) < %.1f && abs(dcaZInCm1) < %.1f && ", dcaZInCm[0], dcaZInCm[0]);
      }

      if(var.Contains("dcaXYInCm") || var.Contains("dcaxyincm") || var.Contains("DCAXYINCM")){
         c += Form("dcaXYInCm0 < %.1f && dcaXYInCm1 < %.1f && ", dcaXYInCm[j], dcaXYInCm[j]);
      }else{
         c += Form("dcaXYInCm0 < %.1f && dcaXYInCm1 < %.1f && ", dcaXYInCm[0], dcaXYInCm[0]);
      }

      if(var.Contains("vertexZInCm") || var.Contains("vertexzincm") || var.Contains("VERTEXZINCM") || var.Contains("Vz")){
         c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[j]);
      }else{
         c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[0]);
      }
   }

   //c += Form("isBemcHit0 == 1 && isBemcHit1 == 1 && "); // require both electrons to have matched BEMC hit

   //remove the last " && "
   if(c.EndsWith("&& ") || c.EndsWith(" &&")){
      c.Remove(c.Length() - 3, 3);
   }

   // cout << "Condition for variable " << var << ": " << c << " and j = " << j << endl;

   return c;

}


void Plot::loadStarlightTree( TString filename, TString treename ){


   if(filename == ""){
      filename = MCFileName;
      treename = MCTreeName;
   }

   if(DEBUGMODE) cout << "Connecting starlight file..." << endl;


   starlightFile = shared_ptr<TFile>(new TFile(filename, "READ"));
   if(!starlightFile || starlightFile->IsZombie() || !starlightFile->IsOpen()){
      cerr << "ERROR in Plot::loadStarlightTree(): Could not get file with starlight data. Returning." << endl;
      return;
   }
   
   
   if(DEBUGMODE) cout << "Checking starlight tree..." << endl;
   starlightTree = (TTree*)starlightFile->Get(treename);
   if(!starlightTree || starlightTree->IsZombie() || starlightTree->GetEntries() <= 0){
      cerr << "ERROR in Plot::loadStarlightTree(): Could not get starlight tree. Returning." << endl;
      return;
   }

   if(DEBUGMODE) cout << "Successfully connected starlight tree" << endl;

}
