#ifndef PlotsManager_h
#define PlotsManager_h


#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"
#include "PlotAnaV0.h"
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>

using namespace std;
//using namespace UTIL;


TFile *mOutFile;
const char* nameOfTree;
std::vector<std::pair<TString, bool>> plots;
//------------------------------------------------
//only temporary until RunDefPlots.h works
/*
const char* nameOfTree = "recTreeV0";
std::vector<std::pair<TString, bool>> plots = {
    {"invMassK0s", true},
    {"invMassLambda", true},
    {"hNSigmaPiPcorr", false},
    {"hNSigmaPiKcorr",false},
    {"hNSigmaPKcorr",false},
    {"hNPairV0", true},
    {"hPosZ", true},
    {"hPosZCut", true},
    {"hVtxDiff", true},
    {"hDcaDaughters", true},
    {"hDcaDaughtersCut", true},
    {"hDcaBeamline", true},
    {"hDcaBeamlineCut", true},
    {"hPointingAngle", true},
    {"hPointingAngleCut", true},
    {"hDecayLength", true},
    {"hDecayLengthCut", true},
    {"hEta",true},
    {"hPt",true},
    {"hNfitHits",true},
    {"hNhitsDEdx", true},
    {"hAnalysisFlow", true},
    {"hDecayLPointingA", true},
    {"hDecayLPointingACut", true}

};
bool runAnaV0 = true;
const char* YAxisDescription = "counts";
///-----------------------------------------
*/
Int_t nInputFiles;

TH2F *hDecayLPointingA, *hDecayLPointingACut;
TH2F *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPKcorr;
TH1D *hNPairV0, *hPosZ, *hPosZCut, *hVtxDiff;
TH1D *hPt, *hEta, *hNhitsDEdx, *hNfitHits, *hAnalysisFlow;
//topologyCuts
TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;

vector<TH2F*> PIDplots;
vector<TH1D*> histograms;

PlotAnaV0* mPlot;

TFile *CreateOutputFile(const string& out);
bool ConnectInputHist(int argc, char** argv);
//bool defineAnalysis();
bool DefineHists(const string& input); 
bool fillHist(const string& input);
void Clear();
bool saveHists();



bool ConnectInputHist(int argc, char** argv){
	//add histograms together if there are multiple input .root files

   const string& input = argv[1];
   TFile *inFile;
   if(input.find(".root") != string::npos){
    	cout << "Input from only one root file." << endl;

      DefineHists(input);

      if(!fillHist(input.c_str())){
         cerr << "Couldn't properly fill histograms from the only file." << endl;
         return false;
      }

    	return true;
   } 
   else if(input.find(".list") != string::npos){
    	

    	ifstream instr(input.c_str());
      if (!instr.is_open()){
        	cout<< "Couldn't open: "<< input.c_str() <<endl;
         return false;
      }
      //cout << "Input from: " << input.c_str() << endl;

      int iFile = 0;
      string line;
    	while(getline(instr, line)) {
         if (line.empty())
            continue;

         //define histograms based on the first file
         if(iFile == 0){
            DefineHists(line);
            ++iFile;
            continue;
         }

	      inFile = TFile::Open(line.c_str(), "read");
	      if(!inFile){
	         cout << "Couldn't open: " << line.c_str() << endl;
	         return false;
	      } 

         if(!fillHist(line)){
            cerr << "Couldn't properly fill histograms from file: " << line.c_str() << endl;
            return false;
         }
         inFile->Close();

	   }//while

   }
   else{
   	cerr << "Didn't find .list or .root in input file. Returning..." << endl;
   	return false;
   }

	saveHists();

	return true;
}

bool fillHist(const string& input){

   TH1D *hist;
   TH2F *PIDplot;

   TFile* inputFile = TFile::Open(input.c_str(), "read");

   if(!inputFile){
      cout << "Couldn't open: " << input.c_str() << endl;
      return false;
   } 

   for (unsigned int iHist = 2; iHist < plots.size(); ++iHist){
   	if (!plots[iHist].second)
   		continue;

	  	if (iHist == 2 || iHist == 3 || iHist == 4 || iHist == 20 || iHist == 21 || iHist == 22){ //considering only 3 PID plots
         inputFile->cd("/PID");
         PIDplot = (TH2F*)gDirectory->Get( plots[iHist].first );
         hist = (TH1D*)gDirectory->Get( plots[iHist].first );
	  		if (PIDplot){
	  		   PIDplots[iHist-2]->Add(PIDplot);
	  		} else if(hist){
            histograms[iHist-17]->Add(hist);
         }else{
	  			cerr << "Couldn't get PID histogram " << plots[iHist].first << " from file."<< endl;
	  			return false;
         }

	  	} else if(iHist == 22 || iHist == 23){
         inputFile->cd("/topologyCuts");
         PIDplot = (TH2F*)gDirectory->Get( plots[iHist].first);
         if (!PIDplot){
            cerr << "Couldn't get PID histogram " << plots[iHist].first << " from file."<< endl;
            return false;
         }
         PIDplots[iHist-19]->Add(PIDplot);

      }else {
         inputFile->cd("/cuts");
	  		hist = (TH1D*)gDirectory->Get( plots[iHist].first);
	  		if(!hist){
            inputFile->cd("/topologyCuts");
            hist = (TH1D*)gDirectory->Get(plots[iHist].first);
         }
         if(!hist){
            inputFile->cd();
            hist = (TH1D*)gDirectory->Get( plots[iHist].first);
         }
         if(!hist){   
	  			cerr << "Couldn't get histogram " << plots[iHist].first << " from file."<< endl;
	  			return false;
	  		}
	  		histograms[iHist-5]->Add(hist);
	  	}//else
      //cout << "Just added histogram: " << plots[iHist].first << " from: " << input.c_str() << endl;
	}//for
   inputFile->Close();
   return true;
}//fillHist

/*
bool histDefinition(const string& input) {
   vector<TString> hists;
   vector<TH1D*> hists1D;
   vector<TH2F*> hists2D;

   // open input file
   TFile* inputFile = TFile::Open(input.c_str(), "read");
   if(!inputFile){
      cout << "Couldn't open: " << input.c_str() << endl;
      return false;
   }  


   TList* list = file->GetListOfKeys();
   TIter next(list);
   TKey* key;
   // Iterate over all keys in the file and if its TH1D*, save them 
   while ((key = dynamic_cast<TKey*>(next()))) {
      TObject* obj = key->ReadObj();

      // Check if the object is a TH1D
      if (dynamic_cast<TH1D*>(obj)) {
         TH1D* h1 = dynamic_cast<TH1D*>(obj);
         TString name = h1->GetName();
         hists.push_back(name);
         hists1D.push_back(h1);
      }
   }

   while ((key = dynamic_cast<TKey*>(next()))) {
      TObject* obj = key->ReadObj();

      // Check if the object is a TH2F
      else if (dynamic_cast<TH2F*>(obj)) {
         TH2F* h2 = dynamic_cast<TH2F*>(obj);
         TString name = h2->GetName();
         h2File->cd();
         h2->Write();
         h2Names.push_back(name);
      }
   }


}
*/

bool DefineHists(const string& input){

	//mOutFile->mkdir("PID")->cd();

   TFile* inputFile = TFile::Open(input.c_str(), "read");
   if(!inputFile){
      cout << "Couldn't open: " << input.c_str() << endl;
      return false;
   }  

	inputFile->cd("/PID");

   hNSigmaPiPcorr = (TH2F*)gDirectory->Get(plots[2].first);
   PIDplots.push_back(hNSigmaPiPcorr);

   hNSigmaPiKcorr = (TH2F*)gDirectory->Get(plots[3].first);
   PIDplots.push_back(hNSigmaPiKcorr);

   hNSigmaPKcorr = (TH2F*)gDirectory->Get(plots[4].first);
   PIDplots.push_back(hNSigmaPKcorr);

   for (int i = 0; i < PIDplots.size(); ++i){
      if(!PIDplots[i])
         cout << "Did not get PID " << i << endl; 
   }

   inputFile->cd("/cuts");

   hNPairV0 = (TH1D*)gDirectory->Get(plots[5].first);
   histograms.push_back(hNPairV0);

   hPosZ = (TH1D*)gDirectory->Get(plots[6].first);
   histograms.push_back(hPosZ);

   hPosZCut = (TH1D*)gDirectory->Get(plots[7].first);
   histograms.push_back(hPosZCut);

   hVtxDiff = (TH1D*)gDirectory->Get(plots[8].first);
   histograms.push_back(hVtxDiff);

   inputFile->cd("/topologyCuts");

   hDcaDaughters = (TH1D*)gDirectory->Get(plots[9].first);
   histograms.push_back(hDcaDaughters);

   hDcaDaughtersCut = (TH1D*)gDirectory->Get(plots[10].first);
   histograms.push_back(hDcaDaughtersCut);

   hDcaBeamline = (TH1D*)gDirectory->Get(plots[11].first);
   histograms.push_back(hDcaBeamline);

   hDcaBeamlineCut = (TH1D*)gDirectory->Get(plots[12].first);
   histograms.push_back(hDcaBeamlineCut);

   hPointingAngle = (TH1D*)gDirectory->Get(plots[13].first);
   histograms.push_back(hPointingAngle);

   hPointingAngleCut = (TH1D*)gDirectory->Get(plots[14].first);
   histograms.push_back(hPointingAngleCut);

   hDecayLength = (TH1D*)gDirectory->Get(plots[15].first);
   histograms.push_back(hDecayLength);

   hDecayLengthCut = (TH1D*)gDirectory->Get(plots[16].first);
   histograms.push_back(hDecayLengthCut);

   inputFile->cd("/cuts");

   hEta = (TH1D*)gDirectory->Get(plots[17].first);
   histograms.push_back(hEta);

   hPt = (TH1D*)gDirectory->Get(plots[18].first);
   histograms.push_back(hPt);

   hNfitHits = (TH1D*)gDirectory->Get(plots[19].first);
   histograms.push_back(hNfitHits);

   hNhitsDEdx = (TH1D*)gDirectory->Get(plots[20].first);
   histograms.push_back(hNhitsDEdx);

   inputFile->cd();

   hAnalysisFlow = (TH1D*)gDirectory->Get(plots[21].first);
   histograms.push_back(hAnalysisFlow);

   inputFile->cd("/topologyCuts");

   hDecayLPointingA = (TH2F*)gDirectory->Get(plots[22].first);
   PIDplots.push_back(hDecayLPointingA);

   hDecayLPointingACut = (TH2F*)gDirectory->Get(plots[23].first);
   PIDplots.push_back(hDecayLPointingACut);

   return true;
}


bool defineAnalysis(){
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
    } else {
    	cout << "All analyses set to false. Returning..." << endl;
    	return false;
    }
}//defineAnalysis


bool saveHists() {
   mOutFile->cd();

   for (int i = 0; i < PIDplots.size(); ++i){
      PIDplots[i]->Write();
   }

   for (int i = 0; i < histograms.size(); ++i){
      histograms[i]->Write();
   }
}


//_____________________________________________________________________________
TFile *CreateOutputFile(const string& out) {

   TFile *outputFile = TFile::Open(out.c_str(), "recreate");
   if(!outputFile) 
      return 0x0;

   return outputFile;
}//CreateOutputFile


void Clear() {
   if(mOutFile){
      mOutFile->Close();
      delete mOutFile;
   }

}



#endif

