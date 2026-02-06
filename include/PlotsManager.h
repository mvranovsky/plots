#ifndef PlotsManager_h
#define PlotsManager_h


#include "Config.h"
#include "PlotAnaV0.h"
//#include "PlotGeneral.h"
#include "PlotAnaV0Mult.h"
#include "PlotTofEff.h"
#include "PlotTofEffMult.h"
#include "PlotAnaJPsi.h"
#include "PlotGoodRun.h"
#include "PlotEmbeddingJPsi.h"
#include "CrossSectionMaker.h"
#include "PlotZeroBias.h"
#include "Libraries.h"
#include "PlotBemcEfficiency.h"
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <cstring>

using namespace std;
//using namespace UTIL;


const char* nameOfTree;
std::vector<std::pair<TString, bool>> plots;


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


vector<string> getFilenames(const string& input){

   vector<string> filenames;
   //cout << "Input file: " << input.c_str() << endl;
   if(input.find(".root") != string::npos){
      filenames.push_back(input);
   } 
   else if(input.find(".list") != string::npos ){
      ifstream instr(input.c_str());
      if (!instr.is_open()){
         cout<< "Couldn't open: "<<input.c_str()<<endl;
         return {};
      }
      string line;
      while(getline(instr, line)) {
         if(line.empty())
            continue;
         shared_ptr<TFile> inputFile(new TFile(line.c_str(), "read"));
         if(!inputFile || inputFile->IsZombie()){
            cout << "Couldn't open: " << line.c_str() << endl;
            return {};
         }
         filenames.push_back(line);
      }//while
   }//else if

   return filenames;
}

vector<TDirectory*> getDirectories(shared_ptr<TFile>& file) {
   vector<TDirectory*> directories;
   if (!file || file->IsZombie()) {
      cerr << "Error: File is not open or is a zombie." << endl;
      return directories;
   }

   TIter nextKey(file->GetListOfKeys());
   TKey* key;
   while ((key = (TKey*)nextKey())) {
      TClass* cls = TClass::GetClass(key->GetClassName());
      if (!cls || !cls->InheritsFrom("TDirectory")) continue;
      TDirectory* dir = static_cast<TDirectory*>(key->ReadObj());
      directories.push_back(dir);
   }
   return directories;

}


bool connectHists(const char* inputPath, const char* outputFilename, bool keepExistingHists = false) {
   map<string, TH1*> mergedTH1Hists;
   map<string, TH2*> mergedTH2Hists;


   std::vector<string> filenames = getFilenames(inputPath);

   if(filenames.size() == 0 ){
      cout << "Empty vector of filenames" << endl;
      return false;
   }

   for (const auto& filename : filenames) {
      std::shared_ptr<TFile> file( new TFile(filename.c_str(), "READ") );
      if (!file || file->IsZombie()) {
         std::cerr << "Error opening file: " << filename << std::endl;
         continue;
      }
      TIter nextKey(file->GetListOfKeys());
      TKey* key;

      for(const auto& dir : getDirectories(file)) {
         if (!dir) continue;
         // Iterate over all keys in the directory
         TIter nextKey(dir->GetListOfKeys());
         while ((key = (TKey*)nextKey())) {
            TObject* obj = key->ReadObj();
            if (!obj) continue;

            // Check if it is a histogram (TH1 base class)
            TClass* cls = TClass::GetClass(key->GetClassName());
            if (!cls) continue;
            if (cls->InheritsFrom("TH1D") || cls->InheritsFrom("TH1F") || cls->InheritsFrom("TH1I")) {
               TH1* hist = static_cast<TH1*>(obj);

               string histName;
               if(TString(hist->GetName()).Contains("hVtxZFillNum") ){
                  histName = ( TString("hVtxZFillNum/") + TString(hist->GetName()) ).Data();
               }else{
                  histName = ( TString(dir->GetName()) + TString("/") + TString(hist->GetName()) ).Data();
               }


               // Merge histograms with the same name
               if (mergedTH1Hists.find(histName) == mergedTH1Hists.end()) {
                  // First instance, clone the histogram
                  mergedTH1Hists[histName] = static_cast<TH1*>(hist->Clone());
                  mergedTH1Hists[histName]->SetDirectory(nullptr); // Prevent ownership transfer
               } else {
                  // Add to existing histogram
                  mergedTH1Hists[histName]->Add(hist);
               }
            }
            // Check if it is a TH2 class
            if (cls->InheritsFrom("TH2D") || cls->InheritsFrom("TH2F") || cls->InheritsFrom("TH2I")) {
               TH2* hist = static_cast<TH2*>(obj);
               string histName = (TString(dir->GetName()) + TString("/") + TString(hist->GetName()) ).Data();

               // Merge histograms with the same name
               if (mergedTH2Hists.find(histName) == mergedTH2Hists.end()) {
                  // First instance, clone the histogram
                  mergedTH2Hists[histName] = static_cast<TH2*>(hist->Clone());
                  mergedTH2Hists[histName]->SetDirectory(nullptr); // Prevent ownership transfer
               } else {
                  // Add to existing histogram
                  mergedTH2Hists[histName]->Add(hist);
               }
            }
         }
      }
   }

   TString opt = "RECREATE";
   if(keepExistingHists) opt = "UPDATE";
   // Save merged histograms to new ROOT file
   std::shared_ptr<TFile> histFile(new TFile(outputFilename, opt));
   if (!histFile || histFile->IsZombie()) {
      std::cerr << "Error creating output file: " << outputFilename << std::endl;
      return false;
   }

   map<TString,int> saved;
   histFile->cd();
   for (const auto& entry : mergedTH1Hists) {
      // get the directory name from the histogram name
      size_t pos = entry.first.find_last_of('/');
      if (pos != string::npos) {
         string dirName = entry.first.substr(0, pos);
         if (!histFile->GetDirectory(dirName.c_str())) {
            histFile->mkdir(dirName.c_str());
         } 
         histFile->cd(dirName.c_str());
      } else {
         histFile->cd();
      }
      if(saved.find(entry.first) == saved.end()){
         saved[entry.first] = 1;
         entry.second->Write();
      }
      histFile->cd();
   }
   for (const auto& entry : mergedTH2Hists) {
      // get the directory name from the histogram name

      size_t pos = entry.first.find_last_of('/');
      if (pos != string::npos) {
         string dirName = entry.first.substr(0, pos);
         if (!histFile->GetDirectory(dirName.c_str())) {
            histFile->mkdir(dirName.c_str());
         }
         histFile->cd(dirName.c_str());
      } else {
         histFile->cd();
      }
      if(saved.find(entry.first) == saved.end()){
         saved[entry.first] = 1;
         entry.second->Write();
      }else{   
         cout << "Histogram " << entry.first << " already saved. Skipping." << endl;
      }
      histFile->cd();
   }
   histFile->Close();
   //histFile.reset();

   cout << "Merged histograms saved in " << outputFilename << endl;
   return true;
}

void runCrossSection(const char* outputFile) {
   if(DEBUGMODE) cout << "Running CrossSectionMaker..." << endl;

   //TString anaDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaJPsi_withRP_1.12.25/merged/StRP_production.list";
   TString anaDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaJPsi_noRP_sysStudy_2.12.25/merged/StRP_production.list";
   TString embedDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/EmbeddingJPsi_sysStudy/merged/StRP_production.list";
   TString goodRunDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaGoodRun_2.12.25/merged/StRP_production.list";
   TString zeroBiasDir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaZeroBias_26.8.25/merged/StRP_production.list";

   if(!connectHists(anaDir, "histFile.root", false))  return;

   if(!connectHists(embedDir, "histFile.root", true))  return;

   if(!connectHists(goodRunDir, "histFile.root", true))  return;

   if(!connectHists(zeroBiasDir, "histFile.root", true))  return;

   CrossSectionMaker *crossSectionMaker = new CrossSectionMaker(anaDir, embedDir, goodRunDir, zeroBiasDir, TString(outputFile));
   crossSectionMaker->Make();
   delete crossSectionMaker;
}


#endif

