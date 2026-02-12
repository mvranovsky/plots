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
		void DrawSTARTag(double xl = 0.7, double yl = 0.9, double xr = 0.9, double yr = 0.9, double textSizeRel = 0.);
		
		// set style for TGraph
		void SetTGraphStyle(TGraph*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);
		void SetTGraphStyle(TGraphAsymmErrors*& graph, Int_t color = kBlue, Int_t markStyle = 20);

		// create plot for a 1D histogram
		void TH1General(TH1*& hist,TString nameOfHist, TString dir = "");  

		// create a plot for a 2D histogram
		void TH2General( TH2*& hist,TString nameOfHist , TString dir = "");

		// load Monte Carlo simulation tree
		void loadMCTree( TString filename, TString treename );
		
		TTree *tree, *bcgTree;
		TTree *MCTree;

		template<typename T1>
		static void SetTH1Style(T1* hist, Int_t color, Int_t markStyle)
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
		}//SetTH1Style


		template<typename T2>
		static void SetTH2Style(T2* hist)
		{
			// Define color gradient (blue → green → yellow)
			const Int_t NRGBs = 3;
			const Int_t NCont = 255;
			Double_t stops[NRGBs]  = { 0.00, 0.50, 1.00 };
			Double_t red[NRGBs]    = { 0.00, 0.00, 1.00 };
			Double_t green[NRGBs]  = { 0.00, 1.00, 1.00 };
			Double_t blue[NRGBs]   = { 1.00, 1.00, 0.00 };

			TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
			gStyle->SetNumberContours(NCont);

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
