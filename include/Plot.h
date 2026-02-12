#ifndef Plot_h
#define Plot_h

// include headers
#include "Libraries.h"
#include "Config.h"

using namespace std;



class Plot{
	//can be accessed from within this class or from derived classes
	public: 
		Plot(const std::string& mInputList, const char* filePath);
		Plot(const std::string& mInputList,
			const std::shared_ptr<TFile>& file);

		Plot(const Plot&) = default;
		virtual ~Plot();
		
		virtual void Make(){cout<<"Hi my name is make"<<endl;};
        virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
		virtual void Finish() {cout << "Hi I should not be there but I am"<<endl;};
		
		bool ConnectInputTree(const string& input, TString nameOfTree,TTree *&TREE); // method for connecting input tree

		// helper methods for ROOT objects
		void CreateCanvas(TCanvas **canvas, TString canvasName, int canvasWidth = 1000 , int canvasHeight = 800);
		void SetGPad(bool isLogY = false, double left = 0.14, double right = 0.07, double bottom = 0.11, double top = 0.06 ); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
		void CreateLegend(TLegend **legend,double xl = 0.74, double yl = 0.74, double xr = 0.97, double yr = 0.89);
		void CreateText(TString writtenText, int textF = 42,double xl = 0.65, double yl = 0.68, double xr = 0.88, double yr = 0.88);
		void SetLineStyle(TLine*& line);
		void CreateLine(double xl, double yl, double xr, double yr);
		void CreateDashedLine(double xl, double yl, double xr, double yr);
		
		// draw a tag in top right corner (description, process, energy of collision)
		void DrawSTARTag(double xl = 0.2, double yl = 0.85, double xr = 0.4, double yr = 0.85, double textSizeRel = 0.);
		
		
		//!!! change according to Gabe's analysis !!!
		TString getCondition(TString var = "", int j = 0);
		TH1D* loadInvMassHist(int numBins, Double_t minRange, Double_t maxRange, TString c, bool removeBcg = true); // function to load invariant mass histogram from tree with given cuts. If removeBcg is true, background histogram will be subtracted from signal histogram
		
		// set style for TGraph
		void SetTGraphStyle(TGraph*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphAsymmErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);

		// create plot for a 1D histogram
		void TH1General(TH1*& hist,TString nameOfHist, TString dir = "");  
		void TH1General(TH1D*& hist,TString nameOfHist, TString dir = "");  
		void TH1General(TH1F*& hist,TString nameOfHist, TString dir = "");  
		void TH1General(TH1I*& hist,TString nameOfHist, TString dir = ""); 
		
		// set style for a 1D histogram
		
		/*
		void SetTH1Style(TH1*& hist, Int_t Color = color, Int_t markStyle = markerStyle);
		void SetTH1Style(TH1D*& hist, Int_t Color = color, Int_t markStyle = markerStyle);
		void SetTH1Style(TH1F*& hist, Int_t Color = color, Int_t markStyle = markerStyle);
		void SetTH1Style(TH1I*& hist, Int_t Color = color, Int_t markStyle = markerStyle);
		*/

		// create a plot for a 2D hist
		void TH2General( TH2*& hist,TString nameOfHist , TString dir = "");
		void TH2General( TH2D*& hist,TString nameOfHist , TString dir = "");
		void TH2General( TH2F*& hist,TString nameOfHist , TString dir = "");
		void TH2General( TH2I*& hist,TString nameOfHist , TString dir = "");

		// set style for a 2D hist
		void SetTH2Style(TH2*& hist);
		void SetTH2Style(TH2D*& hist);
		void SetTH2Style(TH2F*& hist);
		void SetTH2Style(TH2I*& hist);

		// load Monte Carlo simulation tree
		void loadMCTree( TString filename, TString treename );
		
		TTree *tree, *bcgTree;
		TTree *MCTree;

		template <typename T>
		static void SetTH1Style(T*& hist, Int_t color, Int_t markStyle)
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
			hist->SetLineColor(color);
			hist->SetLineStyle(1);
			hist->SetLineWidth(1);  
			hist->SetMarkerSize(2);
			hist->SetMarkerColor(color);
			hist->SetMarkerStyle(markerStyle);
		}//SetHistStyle

	protected:
		string outputPosition;
		string inputPosition;
		shared_ptr<TFile> outFile;  // Shared pointer
		shared_ptr<TFile> histFile;  // This one can stay shared_ptr if Plot owns it
		shared_ptr<TFile> MCFile;

		

		// methods to handle TH1 and TH2 histograms
		bool saveHistograms(TString dir = "");
		vector<pair<TH1*, TString>> GetAllTH1();
		vector<pair<TH2*, TString>> GetAllTH2();

		TCanvas *canvas;
		TLegend *legend, *legend2;
		TLine *line;

		TLine *l, *l1, *l2, *l3, *l4, *l5, *l6;

		TPaveText *text;
		TChain *chain, *bcgChain;
		TF1 *func, *bcg;
		TH1D *histToFit;








};

#endif
