#ifndef Plot_h
#define Plot_h

// include headers
#include "Libraries.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"

using namespace std;

class Plot{
	public: 
		Plot(TFile* mOutFile, const string mInputList, const char* filePath/*, const char* treeName, std::vector<std::pair<TString, bool>> plotsFromManager*/);
		virtual ~Plot();




		virtual void Make(){cout<<"Hi my name is make"<<endl;};
        virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
   
        //TFile* outFile;
   	protected:
   		//can be accessed from within this class or from derived classes
		
		bool ConnectInputTree(const string& input);
		void CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth, int canvasHeight);
		void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle);
		void SetTH2Style(TH2* hist);
		void CreateLegend(TLegend **legend,double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
		void DrawSTARInternal(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
		void CreateText(TString writtenText, double xl, double yl, double xr, double yr);
		void SetGPad(bool isLogY = false, double left = 0.14, double right = 0.07, double bottom = 0.11, double top = 0.06 ); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
		void DrawFiducial();
		void SetLineStyle(TLine* line);
		void DrawSTARpp510(double xl = 0.2, double yl = 0.9, double xr = 0.4, double yr = 0.9);
		int fitGaussPol2(TH1D **histToFit, Int_t binWidth, Double_t minRange, Double_t maxRange, Double_t pol0 = 0.01, Double_t pol1 = 0.01, Double_t pol2 = 0.01, Double_t amp = 10000,Double_t mean = 0.497, Double_t sigma = 0.001);
		bool defineAnalysis();

		//bool defineAnalysis();
		
		TFile *outFile, *inFile;
		//Util *mUtil;
		const char* nameOfTree;
		std::vector<std::pair<TString, bool>> plots;

		//TTree *mTree[nStudies];
		//RecTree *mRecTree[nStudies];
		vector<TH2F*> PIDplots;
		//vector<TH1D*> histograms;
		const char* outputPosition;
		string inputPosition;

		/*
		//------------------------------------------------
		//only temporary until RunDefPlots.h works
		const char* nameOfTree = "recTreeV0";
		std::vector<std::pair<TString, bool>> plots = {
		    {"invMassK0s", true},
		    {"invMassLambda", true},
		    {"hNSigmaPiPcorr", true},
		    {"hNSigmaPiKcorr",true},
		    {"hNSigmaPKcorr",true},
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
		///-----------------------------------------
		
		TH2F *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPKcorr;
		TH1D *hNPairV0, *hPosZ, *hPosZCut, *hVtxDiff;
		TH1D *hPt, *hEta, *hNhitsDEdx, *hNfitHits;
		//topologyCuts
		TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
		TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
		*/

		TCanvas *canvas;
		TLegend *legend;
		TPaveText *text;
		TTree *tree;
		TChain *chain;
		TF1 *func, *bcg;
		TH1D *histToFit;

		//Util *mUtil;


    	//some preset constants for histograms
		const int widthTypical = 1000;
		const int heightTypical = 800;
		const int markerStyleTypical = 20;
		TString ppSTAR = "p+p #sqrt{s} = 510 GeV";

		const TString yAxisTitle = "counts"; //entries/counts
		const int textAlign = 12;
		const int textFont = 42; 
		const Double_t textSize = 0.03;
		const Double_t labelSize = 0.05;
		const Double_t ticksSize = 0.05;
		const int markerStyle = 20;		


/*
// general stuff
bool runPlots(TString inputFile, TString embedFile);
bool Init(TString inputFile, TString embedFile);
bool ConnectInput( TString inputFile, TFile **file); 
TFile *CreateOutputFile(const string& out);

*/


};

#endif
