#ifndef Plot_h
#define Plot_h

// include headers
#include "Libraries.h"
#include <memory>


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
		
		//helper functions
		void CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth, int canvasHeight);
		void CreateLegend(TLegend **legend,double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
		void DrawSTARInternal(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
		void SetGPad(bool isLogY = false, double left = 0.14, double right = 0.07, double bottom = 0.11, double top = 0.06 ); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
		void CreateText(TString writtenText, double xl, double yl, double xr, double yr);
		void DrawFiducial();
		void SetLineStyle(TLine* line);
		void DrawSTARpp510(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		void DrawSTARpp510JPsi(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		
		// methods to handle TH1 and TH2 histograms
		bool handleHistograms();
		void TH1General(TString nameOfHist,TH1*& hist);
		void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle);
		void TH2General(TString nameOfHist , TH2*& hist);
		void SetTH2Style(TH2* hist);
		vector<pair<TH1*, TString>> GetAllTH1();
		vector<pair<TH2*, TString>> GetAllTH2();
		
		std::unique_ptr<TFile> 	outFile, histFile;

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

		TString ppSTARJPsi[2] = {"p+p #rightarrow p #oplus J/#psi #oplus p", "#sqrt{s} = 510 GeV"};

		const TString yAxisTitle = "counts"; //entries/counts
		const int textAlign = 12;
		const int textFont = 42; 
		const Double_t textSize = 0.03;
		const Double_t labelSize = 0.05;
		const Double_t ticksSize = 0.05;
		const int markerStyle = 20;		

};

#endif
