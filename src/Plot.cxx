#include "Plot.h"

Plot::Plot(const string mInputList, const char* filePath){
   outputPosition = filePath;
   inputPosition = mInputList;
   outFile = unique_ptr<TFile>(new TFile(outputPosition, "RECREATE"));
}

Plot::Plot(const string mInputList, unique_ptr<TFile> &file) : inputPosition(mInputList), outFile(file.get()) {
   outputPosition = outFile->GetName();
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
   textSTAR -> SetFillColorAlpha(0,0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText("STAR Internal");
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl - 0.05, xr, yr - 0.05,"brNDC");
   textpp510 -> SetTextSize(textSize);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillColor(0);
   textpp510 -> SetTextFont(62);
   textpp510->AddText("p+p #sqrt{s} = 510 GeV");
   textpp510 -> Draw("same");
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

   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize+0.01);
   textSTAR -> SetFillColorAlpha(kWhite, 0.0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText("STAR Internal");
   textSTAR -> Draw("same");
 
   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl-0.05, xr, yr-0.05,"brNDC");
   textpp510 -> SetTextSize(textSize );
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(4000);
   textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // white background with 0% opacity
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
   textSTAR -> SetFillColorAlpha(kWhite, 0.0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText("STAR Internal");
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(xl, yl-0.05, xr, yr-0.05,"brNDC");
   textpp510 -> SetTextSize(textSize + textSizeRel);
   textpp510 -> SetTextAlign(33);
   textpp510 -> SetFillStyle(4000);
   textpp510 -> SetFillColorAlpha(kWhite, 0.0);  
   textpp510 -> SetTextFont(42);
   textpp510->AddText("RP zerobias embedding");
   textpp510->AddText("MC: p + p #rightarrow J/#psi");
   textpp510 -> Draw("same");
}

void Plot::CreateText(TString writtenText, int textF,double xl, double yl, double xr, double yr)
{
   text = new TPaveText(xl, yl, xr, yr,"brNDC");
   text -> SetTextSize(textSize);
   text -> SetFillColor(0);
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
   line->SetLineWidth(3);
}//SetLineStyle

void Plot::CreateLine(double xl, double yl, double xr, double yr)
{
   line = new TLine(xl, yl, xr, yr);
   line->SetLineStyle( lineStyle );
   line->SetLineColor( lineColor );
   line->SetLineWidth( lineWidth );
}//CreateLine

bool Plot::ConnectInputTree(const string& input, TString nameOfTree, TTree *&TREE, TTree *&BCGTREE) {

   chain = new TChain(nameOfTree);

   bcgChain = new TChain(nameOfTree + TString("_Bcg"));
   

   vector<TFile*> inputFiles;
   //cout << "Input file: " << input.c_str() << endl;
   if(input.find(".root") != string::npos){   //!!!! something is wrong
         cout << "Input from root file: "<< input << endl;
      TFile *inputFile = new TFile(input.c_str(), "read");
      if(!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()){
         cout<< "Couldn't open input root file..."<<endl;
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
         if(line.find(".root") == string::npos)
            continue;
         
         TFile *inputFile = new TFile(line.c_str(), "read");
         if(!inputFile || inputFile->IsZombie() || !inputFile->IsOpen()){
            cout << "Couldn't open: " << line.c_str() << endl;
            continue;
         }
         inputFiles.push_back(inputFile);  // Store the file so it doesn't get deleted

         currentTree = dynamic_cast<TTree*>( inputFile->Get(nameOfTree) );
         if(currentTree){
            chain->AddFile(line.c_str());
            bcgChain->AddFile(line.c_str());
         } else {
            cout << "Name of tree: " << nameOfTree << endl;
            cout << "Couldn't open .root file with name: " << line.c_str() << endl;
         }  
         lineId++;
         //all histograms should be added together from the input files, saved to .root file, and then the .root file should be opened with pointer inFile. So far i am pointing to the last inputFile, because I use only one. Fuck off i will finish this one day ;)
      }//while

      instr.close();
      
   }//else if
   
   TREE = dynamic_cast<TTree*>( chain );
   BCGTREE = dynamic_cast<TTree*>( bcgChain );

   if(!BCGTREE || BCGTREE->IsZombie() || BCGTREE->GetEntries() == 0) {
      cerr << "Did not load background tree. Continuing without it..." << endl;
   }

   cout << "Input from " << inputFiles.size() << " files..." << endl;
   return true;
}


bool Plot::handleHistograms(TString dir){

   vector<pair<TH1*,TString>> hists1D;
   vector<pair<TH2*,TString>> hists2D;


   hists1D = GetAllTH1(dir);
   if(hists1D.size() == 0){
      cerr << "Couldn't load 1D histograms from file in dir "<< dir << ". Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists1D.size(); ++i){
      if(!hists1D[i].first){
         cerr << "Couldn't load histogram " << hists1D[i].second << " in dir "<< dir << ". Leaving..." << endl;
         return false;
      }
      TH1General(hists1D[i].second, hists1D[i].first, dir);
      cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
   }

   
   // create canvases for all TH2F
   hists2D = GetAllTH2(dir);
   if(hists2D.size() == 0){
      cerr << "Couldn't load 2D histograms from file. Leaving..." << endl;
      return false;
   }
   for (unsigned int i = 0; i < hists2D.size(); ++i){
      if(!hists2D[i].first){
         cerr << "Couldn't load histogram " << hists2D[i].second << ". Leaving..." << endl;
         return false;
      }
      TH2General(hists2D[i].second, hists2D[i].first, dir);
      cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
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
            cerr << "Directory not found: " << dir << endl;
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
                  cout << "Histogram " << h2->GetName() << " already saved. Skipping." << endl;
               }
               histograms.push_back(make_pair(h2, h2->GetName()));
               cout << "Added 2D histogram: " << h2->GetName() << endl;
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
            cerr << "Directory not found: " << dir << endl;
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
                cout << "Added 1D histogram: " << h1->GetName() << endl;
            }
        }

        // Clean up if ReadObj creates new objects
        // Do not delete h2 here; ownership is typically with the directory unless cloned
    }

    return histograms;
}


void Plot::TH1General(TString nameOfHist,TH1*& hist,TString dir ) {

   if (!hist){
      cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   if (nameOfHist == "") {
      CreateCanvas(&canvas, hist->GetName() + TString("_canvas"), widthTypical, heightTypical );
   }else{
      CreateCanvas(&canvas, nameOfHist  + TString("_canvas"), widthTypical, heightTypical );
   }

   canvas->Clear();
   SetHistStyle(hist, kBlack, markerStyleTypical);
   hist->Draw("same hist");
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
   
   if(strcmp(hist->GetName(),"hNTracksRP") == 0){
      canvas->SetLogy();
   }
   
   //DrawSTARpp510JPsi(0.6,0.85,0.93,0.93, 0.01);
   DrawSTARInternal();
   canvas->Update();


   outFile->cd();
   if(dir != ""){
      outFile->cd(dir);
   }
   canvas->Write();
   //hist->Write();
}


void Plot::TH2General(TString nameOfHist , TH2*& hist, TString dir){
   TCanvas* c = nullptr;
   if (nameOfHist == "") {
      CreateCanvas(&c, TString("canvas_") + hist->GetName(), widthTypical, heightTypical );
   }else{
      CreateCanvas(&c, nameOfHist + TString("_canvas"), widthTypical, heightTypical );
   }

   c->Clear();
   SetGPad(false,0.12, 0.16,0.11, 0.06 );
   if (!hist){
      cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
      return;
   }

   SetTH2Style(hist);

   hist->Draw("COLZ");
   DrawSTARpp510JPsi(0.53,0.85,0.82,0.93, 0.01);
   if(nameOfHist.Contains("hRPcorr") ){
      DrawFiducial();
   }
   else if(nameOfHist.Contains("hNSigma")){
         gPad->SetLogz();
         hist->SetMinimum(0.98);
   }


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



TH1D* Plot::loadInvMassHist(int numBins, Double_t minRange, Double_t maxRange, TString c){ 

   // create a canvas that will hold both fits
   CreateCanvas(&canvas, "invMassJPsi", widthTypical, heightTypical);
   SetGPad();

   TH1D *hSignal = new TH1D("invMassJPsi", "invMassJPsi", numBins, minRange, maxRange);
   TH1D *hBcg = new TH1D("invMassJPsiBcg", "invMassJPsiBcg", numBins, minRange, maxRange);

   hSignal->GetXaxis()->SetTitle("m_{ee} [GeV/c^{2}]");

   TString cmd;
   tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), c);
   hSignal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

   bcgTree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), c);
   hBcg->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

   if(hSignal->GetEntries() == 0 ){
      cout << "Could not load signal to histograms" << endl;
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
      c += Form("nHitsDEdx0 > %d && nHitsDEdx1 > %d &&", nHitsDedx[0], nHitsDedx[0]);
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

   c += Form("chiSquarepion > %.d && chiSquareproton > %.d && chiSquarekaon > %.d && ", chiSquarePi[0], chiSquareP[0], chiSquareK[0]);


   if( !(var.Contains("embedding") || var.Contains("Embedding") || var.Contains("EMBEDDING")) ){  //embedding is not able to reconstruct vertex 

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

   //remove the last " && "
   if(c.EndsWith("&& ")){
      c.Remove(c.Length() - 3, 3);
   }

   return c;

}