#ifndef Plot_h
#define Plot_h

// include headers
#include "Libraries.h"



using namespace std;



class Plot{
	public: 
		// can be accessed from outside the class
		Plot(const string mInputList, const char* filePath);
		virtual ~Plot();

		virtual void Make(){cout<<"Hi my name is make"<<endl;};
        virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
   
        //TFile* outFile;
   	protected:
   		//can be accessed from within this class or from derived classes
		
		bool ConnectInputTree(const string& input, TString nameOfTree, bool alsoBcgTree = false);
		void CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth, int canvasHeight);
		void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle);
		void SetTH2Style(TH2* hist);
		void CreateLegend(TLegend **legend,double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
		void DrawSTARInternal(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
		void CreateText(TString writtenText, double xl, double yl, double xr, double yr);
		void SetGPad(bool isLogY = false, double left = 0.14, double right = 0.07, double bottom = 0.11, double top = 0.06 ); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
		void DrawFiducial();
		void SetLineStyle(TLine* line);
		void DrawSTARpp510(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		void DrawSTARpp510JPsi(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		int fitGaussPol2(TH1D **histToFit, Int_t binWidth, Double_t minRange, Double_t maxRange, Double_t pol0 = 0.01, Double_t pol1 = 0.01, Double_t pol2 = 0.01, Double_t amp = 10000,Double_t mean = 0.497, Double_t sigma = 0.001);
		vector<pair<TH1D*, TString>> GetAllTH1D();
		vector<pair<TH2F*, TString>> GetAllTH2F();
		void TH1DGeneral(TString nameOfHist,TH1D* hist);
		void TH2FGeneral(TString nameOfHist , TH2F* hist);
		void Clear();

		//bool defineAnalysis();
		
		std::unique_ptr<TFile> inFile, outFile;
		//Util *mUtil;
		std::vector<std::pair<TString, bool>> plots;

		//TTree *mTree[nStudies];
		//RecTree *mRecTree[nStudies];
		vector<TH2F*> PIDplots;
		//vector<TH1D*> histograms;
		const char* outputPosition;
		string inputPosition;


		TCanvas *canvas;
		TLegend *legend, *legend2;
		TPaveText *text;
		TTree *tree, *bcgTree;
		TChain *chain, *bcgChain;
		TF1 *func, *bcg;
		TH1D *histToFit;

    	//some preset constants for histograms
		const int widthTypical = 1000;
		const int heightTypical = 800;
		const int markerStyleTypical = 20;
		TString ppSTAR = "p+p #sqrt{s} = 510 GeV";
		TString ppSTARJPsi = "p+p #rightarrow p #oplus J/#Psi #oplus p, #sqrt{s} = 510 GeV";

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
