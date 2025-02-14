#ifndef PlotsManager_h
#define PlotsManager_h


#include "PlotAnaV0.h"
//#include "PlotGeneral.h"
#include "PlotAnaV0Mult.h"
#include "PlotTofEff.h"
#include "PlotTofEffMult.h"
#include "PlotAnaJPsi.h"
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <cstring>

using namespace std;
//using namespace UTIL;


TFile *file1, *file2;
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
bool connectHists(const char* inputPath, const char* outputFilename);
vector<string> getFilenames(const char* inputPath);


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
         return false;
      }
      string line;
      while(getline(instr, line)) {
         if(line.empty())
            continue;
         inputFile = TFile::Open(line.c_str(), "read");
         if(!inputFile){
            cout << "Couldn't open: " << line.c_str() << endl;
            return false;
         } 
         filenames.push_back(line);
      }//while
   }//else if

   return filenames;
}



bool connectHists(const char* inputPath, const char* outputFilename) {
   map<string, unique_ptr<TH1>> mergedTH1Hists;
   map<string, unique_ptr<TH2>> mergedTH2Hists;


   std::vector<string> filenames = getFilenames(inputPath.c_str());

   if(filenames.size() == 0 ){
      cout << "Empty vector of filenames" << endl;
      return false;
   }

   for (const auto& filename : filenames) {
      std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
      if (!file || file->IsZombie()) {
         std::cerr << "Error opening file: " << filename << std::endl;
         continue;
      }
      TIter nextKey(file->GetListOfKeys());
      TKey* key;
      while ((key = (TKey*)nextKey())) {
         // Get object class type
         TClass* cls = TClass::GetClass(key->GetClassName());
         if (!cls) continue;
            // Check if it is a histogram (TH1 base class)
            if (cls->InheritsFrom("TH1")) {
                TH1* hist = static_cast<TH1*>(key->ReadObj());
                string histName = hist->GetName();

                // Merge histograms with the same name
                if (mergedTH1Hists.find(histName) == mergedTH1Hists.end()) {
                    // First instance, clone the histogram
                    mergedTH1Hists[histName] = std::unique_ptr<TH1>(static_cast<TH1*>(hist->Clone()));
                    mergedTH1Hists[histName]->SetDirectory(nullptr); // Prevent ownership transfer
                } else {
                    // Add to existing histogram
                    mergedTH1Hists[histName]->Add(hist);
                }
            }
            // Check if it is a TH2 class
            if (cls->InheritsFrom("TH2")) {
                TH2* hist = static_cast<TH2*>(key->ReadObj());
                string histName = hist->GetName();

                // Merge histograms with the same name
                if (mergedTH2Hists.find(histName) == mergedTH2Hists.end()) {
                    // First instance, clone the histogram
                    mergedTH2Hists[histName] = unique_ptr<TH2>(static_cast<TH2*>(hist->Clone()));
                    mergedTH2Hists[histName]->SetDirectory(nullptr); // Prevent ownership transfer
                } else {
                    // Add to existing histogram
                    mergedTH2Hists[histName]->Add(hist);
                }
            }

            // possible to add other objects to merge and save to histfile
        }
    }

    // Save merged histograms to new ROOT file
    std::unique_ptr<TFile> histFile(TFile::Open(outputFilename.c_str(), "RECREATE"));
    if (!histFile || histFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFilename << std::endl;
        return;
    }

    histFile->cd();
    for (const auto& [name, hist] : mergedTH1Hists) {
        hist->Write();
    }

    for (const auto& [name, hist] : mergedTH2Hists) {
        hist->Write();
    }

    std::cout << "Merged histograms saved in " << outputFilename << std::endl;
}

#endif


