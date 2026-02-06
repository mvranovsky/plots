#ifndef Plot_h
#define Plot_h

// include headers
#include "Libraries.h"
#include <memory>
#include "Config.h"
#include "TRandom3.h"

using namespace std;



class Plot{
	public: 
		// can be accessed from outside the class
		Plot(const string mInputList, const char* filePath);
		Plot(const string mInputList, shared_ptr<TFile> file);

		virtual ~Plot();
		
		virtual void Make(){cout<<"Hi my name is make"<<endl;};
        virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
		virtual void Finish() {cout << "Hi I should not be there but I am"<<endl;};
		
		//can be accessed from within this class or from derived classes
		bool ConnectInputTree(const string& input, TString nameOfTree,TTree *&TREE, TTree *&BCGTREE);
		void CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth = 1000 , int canvasHeight = 800);
		void CreateLegend(TLegend **legend,double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
		void DrawSTARInternal(double xl = 0.8, double yl = 0.89, double xr = 0.93, double yr = 0.93); // 0.8, 0.89, 0.93,0.93
		void DrawSTARInternalZB(double xl = 0.8, double yl = 0.89, double xr = 0.93, double yr = 0.93);
		void SetGPad(bool isLogY = false, double left = 0.14, double right = 0.07, double bottom = 0.11, double top = 0.06 ); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
		void CreateText(TString writtenText, int textF = 42,double xl = 0.65, double yl = 0.68, double xr = 0.88, double yr = 0.88);
		void DrawFiducial();
		void SetLineStyle(TLine*& line);
		void DrawSTARpp510(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		void DrawSTARpp510JPsi(double xl = 0.73, double yl = 0.86, double xr = 0.93, double yr = 0.93, double textSizeRel = 0.01);
		void DrawEmbeddingpp510JPsi(double xl = 0.73, double yl = 0.90, double xr = 0.93, double yr = 0.95, double textSizeRel = 0.0);
		void CreateLine(double xl, double yl, double xr, double yr);
		void CreateDashedLine(double xl, double yl, double xr, double yr);
		void DrawElipse(double centerX, double centerY, double radius1, double radius2);
		void DrawEtaVtxZFiducial();
		TString getCondition(TString var = "", int j = 0);
		TH1D* loadInvMassHist(int numBins, Double_t minRange, Double_t maxRange, TString c, bool removeBcg = true); // function to load invariant mass histogram from tree with given cuts. If removeBcg is true, background histogram will be subtracted from signal histogram
		
		void SetTGraphStyle(TGraph*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphAsymmErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void TH1General(TString nameOfHist,TH1*& hist, TString dir = "", TString addOn = "");
		void SetHistStyle(TH1*& hist, Int_t color, Int_t markStyle);
		void TH2General(TString nameOfHist , TH2*& hist, TString dir = "", TString addOn = "");
		void SetTH2Style(TH2*& hist);
		void loadStarlightTree( TString filename = "", TString treename = "" );
		
		TTree *tree, *bcgTree;
		TTree *starlightTree;

	protected:
		shared_ptr<TFile> outFile;  // Shared pointer
		shared_ptr<TFile> histFile;  // This one can stay shared_ptr if Plot owns it
		shared_ptr<TFile> starlightFile;
		const char* outputPosition;
		string inputPosition;
		
	
		TString convertToString(double val);
		int makeInt(double val);
		void nSigmaCorrPlot(int particles, TString condition = "");  // 1 = electron, 2 = pion, 3 = kaon, 4 = proton


		// methods to handle TH1 and TH2 histograms
		bool handleHistograms(TString dir = "", TString addOn = "");
		vector<pair<TH1*, TString>> GetAllTH1(TString dir = "");
		vector<pair<TH2*, TString>> GetAllTH2(TString dir = "");
		void TH1Specifics(TString nameOfHist, TH1*& hist, TString addOn);
		void TH2Specifics(TString nameOfHist, TH2*& hist, TString addOn);
		
		
		vector<TH2F*> PIDplots;
		//vector<TH1D*> histograms;
		
		const int widthTypical = 1000;
		const int heightTypical = 800;

		TCanvas *canvas;
		TLegend *legend, *legend2;
		TLine *line;

		TLine *left01, *left02, *left03, *left04, *left05, *left06;
		TLine *l, *l1, *l2, *l3, *l4, *l5, *l6;
		TEllipse *circle;

		TPaveText *text;
		TChain *chain, *bcgChain;
		TF1 *func, *bcg;
		TH1D *histToFit;

    	//some preset constants for histograms
		const int markerStyleTypical = 20;
		TString ppSTAR = "p+p #sqrt{s} = 510 GeV";

		TString ppSTARJPsi[4] = {"p+p #rightarrow p #oplus J/#psi #oplus", "#sqrt{s} = 510 GeV","p+p #rightarrow #oplus J/#psi #oplus", "#sqrt{s} = 510 GeV"};

		const TString yAxisTitle = "counts"; //entries/counts
		const int textAlign = 12;
		const int textFont = 42; 
		const Double_t textSize = 0.03;
		const Double_t labelSize = 0.05;
		const Double_t ticksSize = 0.05;
		const int markerStyle = 20;	
		const int lineStyle = 1;
		const int lineColor = kBlue;
		const int lineWidth = 2;

		vector<Int_t> colorSet = { 1, 2, 4, 6, 8};
    	vector<Int_t> markerSet = { 20, 24, 25, 26, 28};


		bool noRomanPots;
	

		// 0. value is the default one, 1. is the tight condition, 2. is the loose condition
		Int_t nHitsDedx[3] = {15, 17, 12};
		TString nHitsDedxDescription[3] = {Form("Nominal (N^{dEdx}_{hits} < %d)", nHitsDedx[0]), Form("Tight (N^{dEdx}_{hits} < %d)", nHitsDedx[1]), Form("Loose (N^{dEdx}_{hits} < %d)", nHitsDedx[2])};
		Int_t nHitsFit[3] = {15, 17, 12};
		TString nHitsFitDescription[3] = {Form("Nominal (N^{fit}_{hits} < %d)", nHitsFit[0]), Form("Tight (N^{fit}_{hits} < %d)", nHitsFit[1]), Form("Loose (N^{fit}_{hits} < %d)", nHitsFit[2])};
		Double_t dcaZInCm[3] = {1, 0.7, 1.2};
		TString dcaZInCmDescription[3] = {Form("Nominal (DCA_{Z} < %.1f)", dcaZInCm[0]), Form("Tight (DCA_{Z} < %.1f)", dcaZInCm[1]), Form("Loose (DCA_{Z} < %.1f)", dcaZInCm[2])};
		Double_t dcaXYInCm[3] = {1.5, 1.2, 1.8};
		TString dcaXYInCmDescription[3] = {Form("Nominal (DCA_{XY} < %.1f)", dcaXYInCm[0]), Form("Tight (DCA_{XY} < %.1f)", dcaXYInCm[1]), Form("Loose (DCA_{XY} < %.1f)", dcaXYInCm[2])};

		Int_t chiSquareE[3] = {9, 7, 11};
		TString chiSquareEDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareE[0]), Form("Tight (chi^{2} < %d)", chiSquareE[1]), Form("Loose (chi^{2} < %d)", chiSquareE[2])};
		Int_t chiSquarePi[3] = {10, 12, 8};
		TString chiSquarePiDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquarePi[0]), Form("Tight (chi^{2} < %d)", chiSquarePi[1]), Form("Loose (chi^{2} < %d)", chiSquarePi[2])};
		Int_t chiSquareK[3] = {10, 12, 8};
		TString chiSquareKDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareK[0]), Form("Tight (chi^{2} < %d)", chiSquareK[1]), Form("Loose (chi^{2} < %d)", chiSquareK[2])};
		Int_t chiSquareP[3] = {10,12,8};
		TString chiSquarePDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareP[0]), Form("Tight (chi^{2} < %d)", chiSquareP[1]), Form("Loose (chi^{2} < %d)", chiSquareP[2])};

		Int_t vertexZInCm[3] = {100, 80, 120};
		TString vertexZInCmDescription[3] = {Form("Nominal (|V_{Z}| < %d)", vertexZInCm[0]), Form("Tight (|V_{Z}| < %d)", vertexZInCm[1]), Form("Loose (|V_{Z}| < %d)", vertexZInCm[2])};
		Double_t etaHadron[3] = {0.9, 0.8, 1.0};
		TString etaHadronDescription[3] = {Form("Nominal (|#eta_{e}| < %.1f)", etaHadron[0]), Form("Tight (|#eta_{e}| < %.1f)", etaHadron[1]), Form("Loose (|#eta_{e}| < %.1f)", etaHadron[2])};

		Double_t nSigmaK[3] = {1.3, 1.5, 1.0};
		TString nSigmaKDescription[3] = {Form("Nominal (|n#sigma_{K}| < %.1f)", nSigmaK[0]), Form("Tight (|n#sigma_{K}| < %.1f)", nSigmaK[1]), Form("Loose (|n#sigma_{K}| < %.1f)", nSigmaK[2])};

		Double_t nSigmaP[3] = {1.3, 1.5, 1.0};
		TString nSigmaPDescription[3] = {Form("Nominal (|n#sigma_{p}| < %.1f)", nSigmaP[0]), Form("Tight (|n#sigma_{p}| < %.1f)", nSigmaP[1]), Form("Loose (|n#sigma_{p}| < %.1f)", nSigmaP[2])};
		

};

#endif
