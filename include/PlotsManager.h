#ifndef PlotsManager_h
#define PlotsManager_h


#include "PlotAnaV0.h"
#include "PlotGeneral.h"
#include "PlotAnaV0Mult.h"
#include "PlotTofEff.h"
#include "PlotTofEffMult.h"
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <cstring>

using namespace std;
//using namespace UTIL;


TFile *mOutFile, *file1, *file2;
const char* nameOfTree;
std::vector<std::pair<TString, bool>> plots;

std::vector<std::pair<TH1D*, TString>> hists1D;
std::vector<std::pair<TH2F*, TString>> hists2D;

Int_t nInputFiles;

TH2F *hDecayLPointingA, *hDecayLPointingACut, *hArmenterosPodolanski, *hInvMassEta;
TH2F *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPKcorr;
TH1D *hNPairV0, *hPosZ, *hPosZCut, *hVtxDiff;
TH1D *hPt, *hEta, *hNhitsDEdx, *hNfitHits, *hAnalysisFlow;
//topologyCuts
TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK;

vector<TH2F*> PIDplots;
vector<TH1D*> histograms;

Plot* mPlot;

TFile *CreateOutputFile(const string& out);
void Clear();
bool saveHists();
bool connectHists(int argc, char *argv[]);
vector<pair<TH2F*, TString>> GetAllTH2FHistograms(const string& fileName);
vector<pair<TH1D*, TString>> GetAllTH1DHistograms(const string& fileName);


bool connectHists(int argc, char *argv[]) {

   const string& input = argv[1];
   TFile *inFile;
   if(input.find(".root") != string::npos){
      cout << "Input from only one root file." << endl;

      hists1D = GetAllTH1DHistograms(input);
      hists2D = GetAllTH2FHistograms(input);

      if(hists1D.size() == 0 || hists2D.size() == 0){
         cerr << "Did not load histograms, couldn't find any." << endl;
         return false;
      }
   } else if(input.find(".list") != string::npos){
      ifstream instr(input.c_str());
      if (!instr.is_open()){
         cout<< "Couldn't open: "<< input.c_str() <<endl;
         return false;
      }

      int iFile = 0;
      string line;
      while(getline(instr, line)) {
         if (line.empty())
            continue;

         //define histograms based on the first file
         if(iFile == 0){
            hists1D = GetAllTH1DHistograms(line);
            hists2D = GetAllTH2FHistograms(line);
            ++iFile;
            continue;
         }

         vector<pair<TH1D*, TString>> hists1D_current = GetAllTH1DHistograms(line);
         if(hists1D.size() != hists1D_current.size()){
            cerr << "Size of hists1D is not equal to hists1D_current for file: " << line << endl;
            return false;
         }

         for (int i = 0; i < hists1D.size(); ++i){
            for (int j = 0; j < hists1D_current.size(); ++j){  
               if(hists1D_current[j].second == hists1D[i].second)
                  hists1D[i].first->Add(hists1D_current[j].first);
            }//current1D loop
         }//hists1D loop

         vector<pair<TH2F*, TString>> hists2D_current = GetAllTH2FHistograms(line);
         if(hists2D.size() != hists2D_current.size()){
            cerr << "Size of hists2D is not equal to hists2D_current for file: " << line << endl;
            return false;
         }
         for (int i = 0; i < hists2D.size(); ++i){
            for (int j = 0; j < hists2D_current.size(); ++j){  
               if(hists2D_current[j].second == hists2D[i].second)
                  hists2D[i].first->Add(hists2D_current[j].first);
            }//current2D loop
         }//hists2D loop


      }//while
   }
   else{
      cerr << "Didn't find .list or .root in input file. Returning..." << endl;
      return false;
   }
   
   if(saveHists()){
      return true;
   } else {
      return false;
   }
}



vector<pair<TH2F*, TString>> GetAllTH2FHistograms(const string& fileName) {
    // Open the ROOT file
   TString inputFileName = fileName;

   file2 = TFile::Open(inputFileName);
   if (!file2 || file2->IsZombie()) {
      std::cerr << "Error opening file: " << fileName << std::endl;
      return {};
   }

   // Vector to store pointers to TH1D histograms
   vector<pair<TH2F*,TString>> histograms;

   // Iterate over all keys in the file
   TKey *key;
   TIter next(file2->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      if (TH2F *h1 = dynamic_cast<TH2F*>(obj)) {
         // If the object is a TH2F histogram, add it to the vector
         histograms.push_back(make_pair(h1, h1->GetName()));
      } else {
         // If not a TH1D, delete the object to free memory
         delete obj;
      }
   }

   return histograms;
}


vector<pair<TH1D*, TString>> GetAllTH1DHistograms(const string& fileName) {
    // Open the ROOT file
   TString inputFileName = fileName;
   file1 = TFile::Open(inputFileName);
   if (!file1 || file1->IsZombie()) {
      std::cerr << "Error opening file: " << fileName << std::endl;
      return {};
   }

   // Vector to store pointers to TH1D histograms
   vector<pair<TH1D*,TString>> histograms;

   // Iterate over all keys in the file
   TKey *key;
   TIter next(file1->GetListOfKeys());
   while ((key = (TKey*)next())) {
      // Retrieve the object pointed by the key. Use ReadObj() to avoid memory leaks caused by Clone()
      TObject *obj = key->ReadObj();
      if (TH1D *h1 = dynamic_cast<TH1D*>(obj)) {
         // If the object is a TH2F histogram, add it to the vector
         histograms.push_back(make_pair(h1, h1->GetName()));
      } else {
         // If not a TH1D, delete the object to free memory
         delete obj;
      }

   }

   return histograms;
}

bool saveHists() {

   if (!mOutFile || mOutFile->IsZombie()) {
    std::cerr << "Error: Output file is not open or is corrupted." << std::endl;
    return false; // or handle the error appropriately
   } else{
      cout << "outfile connected successfully." << endl;
   }

   mOutFile->cd();

   
   for (int i = 0; i < hists1D.size(); ++i){
      //cout << "histogram with name " << hists1D[i].first->GetName() << endl;

      if(hists1D[i].first){
         hists1D[i].first->Write();
      }else{
         return false;
      }
   }

   for (int i = 0; i < hists2D.size(); ++i){
      hists2D[i].first->Write();
   }

   
   return true;
}


//_____________________________________________________________________________
TFile *CreateOutputFile(const string& out) {

   TFile *outputFile = TFile::Open(out.c_str(), "RECREATE");
   if(!outputFile) 
      return 0x0;

   return outputFile;
}//CreateOutputFile


void Clear() {
   if(mOutFile){
      mOutFile->Close();
      delete mOutFile;
   }

   if(file1){
      file1->Close();
      delete file1;
   }
   if(file2){
      file2->Close();
      delete file2;
   }
}



#endif


