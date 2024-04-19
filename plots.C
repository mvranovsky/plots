 // c++ headers
#include <iostream>
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <unistd.h>
#include <cstdlib>
#include <string.h>
#include <string>
#include <utility> 


// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TDirectory.h"


using namespace std;

//TString namesOfPlots[] = {TString("invMassK0s"), TString("invMassLambda"), TString("PIDPlots"), TString("vertexZ"), TString("vertexZCut"), TString("vertexZDiff"), TString("numberOfPairs")};
//vector<bool> plotGraphs = {true, true, true, true, true, true, true};


//define all histograms here which you want to run
vector<pair<TString, bool>> plotsToCreate = {
    {"invMassK0s", true},
    {"invMassLambda", true},
    {"PIDPlots", true},
    {"vertexZ", true},
    {"vertexZCut", true},
    {"vertexZDiff", true},
    {"numberOfPairs", true}
};



TString input = "/home/michal/Desktop/STAR/data/dataV02/output.root";
TString outputDirectory = "/home/michal/Desktop/STAR/figuresV0/";
TString nameOfTree = "recTreeV0";

TString selection[] = {TString("CPT"), TString("N_{hits}^{fit} > 25"), TString("N^{dEdx}_{hits} > 15"), TString("#eta < 0.7"), TString("V0 Selection"), TString("|Vertex_{z}| < 80cm"), TString("totQ = 0")};
TString PIDplots[] = {TString("hNSigmaPiPcorr"), TString("hNSigmaPiKcorr"),TString("hNSigmaPKcorr") };
enum NUMBEREDSTATES {STATE0, STATE1, STATE2, STATE3, STATE4, nStates};

/*
//initializing variables which will be used later
TCanvas* c = nullptr;
TPaveText *text1 = nullptr, *text2 = nullptr, *text3 = nullptr, *textSTAR = nullptr;
TLine *line1 = nullptr, *line2 = nullptr, *line3 = nullptr;
TLegend* legend = nullptr;
*/


//settings for plots
int widthTypical = 800;
int heightTypical = 800;
TString yAxisDescription = "entries"; //entries/counts
int textAlign = 12;
int textFont = 42; 
Double_t textSize = 0.06;
Double_t labelSize = 0.05;
int markerStyle = 20;

//settings for margins
Double_t marginLeft = 0.15;
Double_t marginRight = 0.05;
Double_t marginBottom = 0.10;
Double_t marginTop = 0.05;


//settings for lines and arrows
int lineStyle = 9;
int lineWidth = 2; 
int lineColor = kRed;

//settings for STAR text
int ppSTARfont = 22;
double ppSTARsize = 0.04;
TString ppSTAR = "p+p #sqrt{s} = 510 GeV";


//functions to simplify all processes
void createCanvas(TCanvas*& c,int index , int width, int height);
void addText(TPaveText*& text,TString content, Double_t height1, Double_t height2, Double_t width1, Double_t width2);
void STARpp510(TPaveText*& text, Double_t height1, Double_t height2, Double_t width1, Double_t width2);
void drawLegend(TLegend*& legend, TH1D* hist, TH1D* bcg ,Double_t width1, Double_t width2,Double_t height1, Double_t height2 );
void writeCuts(TPaveText*& text, Double_t width1, Double_t width2,Double_t height1, Double_t height2 );
void drawLine(TLine*& line, Double_t height1, Double_t height2, Double_t width1, Double_t width2);
void drawArrow(TArrow*& line, Double_t height1, Double_t height2, Double_t width1, Double_t width2);
void histStyle(TH1D*& hist, TString xAxisDescription, Double_t xOffset, Double_t yOffset, int color, bool fillColor ,bool logy);
//void clear();


//functions for specific histograms
void invMassK0s(TTree *tree, Double_t minRange, Double_t maxRange, int nBins);
void invMassLambda(TTree *tree, Double_t minRange, Double_t maxRange, int nBins);
void vertexZ(TFile *data, TString nameOfHist ,Double_t minRange, Double_t maxRange, int nBins);
void PIDplot(TFile *data, TString nameOfHist ,Double_t minRange[], Double_t maxRange[], Int_t nBins[]);


void plots(){

	cout << "Initializing plots.C to create graphs from input: " << input << endl;
	cout << "Name of tree: " << nameOfTree << endl;
	cout << "Saving the created plots to directory: " << outputDirectory << endl;

	
	TFile *data = TFile::Open(input);
	if (!data){
		cout<<"Error: cannot open datafile: "<<input<<endl;
		return;
	}

    TTree *tree = dynamic_cast<TTree*>(data->Get(nameOfTree)); //tree name
    if (!tree){
    	cout<<"Error: cannot open tree with name: "<< nameOfTree <<endl;
		return;
    }


	cout <<"Beginning to plot."<< endl; 
	if (plotsToCreate[0].second) {//invMassK0s
		cout << "Plotting invariant mass distribution of pion pairs. " << endl;
		invMassK0s(tree, 0.2, 2, 180); 
	}

    if (plotsToCreate[1].second) { //invMassLambda
        cout << "Plotting invariant mass distribution of pion proton pairs. " << endl;
        invMassLambda(tree, 1,3,200);
    }

    if (plotsToCreate[2].second){ //PID plots

        cout << "Plotting PID correlation plots." << endl;
        for (int i = 0; i < sizeof(PIDplots); ++i){
            //PIDplot(data, PIDplots[i], minRange, maxRange, nBins);
        }
    }
    if (plotsToCreate[3].second){
        cout << "Plotting z position of vertex." << endl;
        vertexZ(data, "cuts/hPosZ", -150, 150, 60);
    }



}

void invMassK0s(TTree *tree, Double_t minRange, Double_t maxRange, int nBins){
	
    TCanvas* c;
	createCanvas(c,0, widthTypical, heightTypical);

    TH1D *signalFinal = new TH1D("signalK0s", "signal", nBins, minRange, maxRange);
    TH1D *bcgFinal = new TH1D("bcgK0s", "bcg", nBins, minRange, maxRange);

    TString cmd, condition;
    for (int i = 0; i < nStates; ++i){

    	cmd = TString("invMass") + i + TString(">>hist(") + nBins + TString(",") + minRange + TString(",") + maxRange + TString(")");
    	condition = TString("pairID") + i + TString(" == 0") + TString(" && totQ") + i + TString(" == 0");
    	cout << cmd << endl;
    	cout << condition << endl;
    	tree->Draw(cmd, condition);
    	signalFinal->Add((TH1D*)gPad->GetPrimitive("hist"));

	    cmd = TString("invMass") + i + TString(">>bcg(") + nBins + TString(",") + minRange + TString(",") + maxRange + TString(")");
    	condition = TString("pairID") + i + TString(" == 0") + TString(" && totQ") + i + TString(" != 0");
    	cout << cmd << endl;
    	cout << condition << endl;    	
    	tree->Draw(cmd, condition);
    	cmd = TString("bcg");
    	bcgFinal->Add((TH1D*)gPad->GetPrimitive(cmd));
    }
	c->Clear();



    histStyle(signalFinal, "m_{#pi #pi} [GeV/c^{2}]",0.9,1.4, kBlack, true, false);
    histStyle(bcgFinal, "m_{#pi #pi} [GeV/c^{2}]",0.9,1.4, kRed, true, false);
    TPaveText *text1, *text2;
    TLegend *leg;
    STARpp510(text1, 0.65,0.9, 0.85,0.9);
    drawLegend(leg, signalFinal, bcgFinal, 0.45, 0.8, 0.7, 0.65);
    //writeCuts(text2, 0.7, 0.9, 0.82, 0.3);
    signalFinal->Draw("same");
    bcgFinal->Draw("same");
    text1->Draw("same");
    //text2->Draw("same");
    leg->Draw("same");

    //c->Update();
    c->SaveAs(outputDirectory + plotsToCreate[0].first + TString(".pdf"));  
    c->Close();

    //clear();
	
}


void invMassLambda(TTree *tree, Double_t minRange, Double_t maxRange, int nBins){
    
    TCanvas* c;
    createCanvas(c,0, widthTypical, heightTypical);

    TH1D *signalFinal = new TH1D("signalK0s", "signal", nBins, minRange, maxRange);
    TH1D *bcgFinal = new TH1D("bcgK0s", "bcg", nBins, minRange, maxRange);

    TString cmd, condition;
    for (int i = 0; i < nStates; ++i){

        cmd = TString("invMass") + i + TString(">>hist(") + nBins + TString(",") + minRange + TString(",") + maxRange + TString(")");
        condition = TString("pairID") + i + TString(" == 1") + TString(" && PairID") + i + TString(" == 2") + TString(" && totQ") + i + TString(" == 0");
        cout << cmd << endl;
        cout << condition << endl;
        tree->Draw(cmd, condition);
        signalFinal->Add((TH1D*)gPad->GetPrimitive("hist"));

        cmd = TString("invMass") + i + TString(">>bcg(") + nBins + TString(",") + minRange + TString(",") + maxRange + TString(")");
        condition = TString("pairID") + i + TString(" == 1")+ TString(" && PairID") + i + TString(" == 2") + TString(" && totQ") + i + TString(" != 0");
        cout << cmd << endl;
        cout << condition << endl;      
        tree->Draw(cmd, condition);
        cmd = TString("bcg");
        bcgFinal->Add((TH1D*)gPad->GetPrimitive(cmd));
    }
    c->Clear();



    signalFinal->Draw("same");
    histStyle(signalFinal, "m_{p #pi} [GeV/c^{2}]",0.9,1.4, kBlack, false, false);
    histStyle(bcgFinal, "m_{p #pi} [GeV/c^{2}]",0.9,1.4, kRed, false, false);
    TPaveText *text1, *text2;
    TLegend *leg;
    STARpp510(text1, 0.65,0.9, 0.85,0.9);
    drawLegend(leg, signalFinal, bcgFinal, 0.45, 0.8, 0.7, 0.65);
    //writeCuts(text2, 0.7, 0.9, 0.82, 0.3);
    bcgFinal->Draw("same");
    text1->Draw("same");
    //text2->Draw("same");
    leg->Draw("same");

    c->SaveAs(outputDirectory + plotsToCreate[1].first + TString(".pdf"));  
    c->Close();

}

void PIDplot(TFile *data, TString nameOfHist ,Double_t minRange[], Double_t maxRange[], Int_t nBins[]){
    
    TCanvas* c;
    createCanvas(c,0, widthTypical, heightTypical);

    TH2D* hist = (TH2D*)data->Get(nameOfHist);
    if (!hist){
        cerr << "Cannot find histogram with the name " << nameOfHist << ". Continuing with other plots." << endl;
        return;
    }

    TPaveText *text1;
    STARpp510(text1, 0.65,0.85, 0.85,0.9);
    //addText(text2, "|Vertex_{z}| < 80 cm" 0.65,0.85, 0.8,0.7);
    //histStyle(hist, "Vertex_{z} [cm]", 0.9,1.4,kBlack, false, true);

    hist->Draw("colz");
    text1->Draw("same");



    c->SaveAs(outputDirectory + nameOfHist + TString(".pdf"));  


}

void vertexZ(TFile *data, TString nameOfHist, Double_t minRange, Double_t maxRange, int nBins){

    TCanvas* c;
    createCanvas(c,0, widthTypical, heightTypical);


    TH1D* hist = (TH1D*)data->Get(nameOfHist);
    if (!hist){
        cerr << "Cannot find histogram with the name " << nameOfHist << ". Continuing with other plots." << endl;
        return;
    }


    TPaveText *text1, *text2;
    STARpp510(text1, 0.65,0.85, 0.85,0.9);
    //addText(text2, "|Vertex_{z}| < 80 cm" 0.65,0.85, 0.8,0.7);
    histStyle(hist, "Vertex_{z} [cm]", 0.9,1.4,kBlack, false, true);
    hist->Draw("same");
    text1->Draw("same");

    c->SaveAs(outputDirectory + plotsToCreate[3].first + TString(".pdf"));  
    c->Close();

}




void histStyle(TH1D*& hist, TString xAxisDescription, Double_t xOffset, Double_t yOffset, int color, bool fillColor ,bool logy){

	hist->GetXaxis()->SetTitleFont(textFont);
    hist->GetYaxis()->SetTitleFont(textFont);
    hist->GetXaxis()->SetLabelFont(textFont);
    hist->GetYaxis()->SetLabelFont(textFont);
    hist->GetXaxis()->SetLabelSize(labelSize);
    hist->GetYaxis()->SetLabelSize(labelSize);
    hist->GetXaxis()->SetTitleSize(labelSize);
    hist->GetYaxis()->SetTitleSize(labelSize);
    hist->GetXaxis()->SetTitleOffset(xOffset);
    hist->GetYaxis()->SetTitleOffset(yOffset);
    hist->GetYaxis()->SetTitle(yAxisDescription);
    hist->GetXaxis()->SetTitle(xAxisDescription);


    gStyle->SetOptStat("");


    hist->SetMarkerColor(color);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(markerStyle);
    hist->SetLineColor(color);
    hist->SetLineStyle(1);
    hist->SetLineWidth(1);
    if (fillColor){
		//hist->SetFillColorAlpha(color, 0.3);
		hist->SetMaximum(3500000);
	}

}

void createCanvas(TCanvas*& c, int index , int width, int height) {
	c = new TCanvas("c", plotsToCreate[index].first, width, height);
    c->SetBottomMargin(marginBottom);
    c->SetLeftMargin(marginLeft);
    c->SetRightMargin(marginRight);
    c->SetTopMargin(marginTop);
}

void addText(TPaveText*& text,TString content, Double_t width1, Double_t width2, Double_t height1, Double_t height2){
    text = new TPaveText(width1,height1,width2,height2, "NDC");
    text->SetTextSize(textSize);
    text->SetTextAlign(textAlign);
    text->SetFillColor(0);
    text->SetTextFont(textFont);
    text->AddText(content);
    //text->Draw("same");
}

void STARpp510(TPaveText*& text, Double_t width1, Double_t width2,Double_t height1, Double_t height2){
	text = new TPaveText(width1,height1,width2,height2, "NDC");
    text->SetTextSize(ppSTARsize);
    text->SetTextAlign(textAlign);
    text->SetFillColor(0);
    text->SetTextFont(ppSTARfont);
    text->AddText(ppSTAR);
    //text->Draw("same");
}

void drawLegend(TLegend*& legend, TH1D* hist, TH1D* bcg ,Double_t width1, Double_t width2,Double_t height1, Double_t height2 ){
	legend = new TLegend(width1, height1, width2, height2);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(textSize-0.02);
    legend->SetTextFont(textFont);
    legend->SetMargin(0.1);
    legend->AddEntry(hist, "Data (unlike sign pairs)","ple");
    legend->AddEntry(bcg, "Data (like sign pairs)","ple");
    //legend->Draw("same");
}
void writeCuts(TPaveText*& text, Double_t width1, Double_t width2,Double_t height1, Double_t height2 ){
    text = new TPaveText(width1,height1,width2,height2, "NDC");
    text->SetTextSize(textSize);
    text->SetTextAlign(textAlign);
    text->SetFillColor(0);
    text->SetTextFont(textFont);
    for (int i = 0; i < sizeof(selection); ++i){
    	text->AddText(selection[i]);
    }
    //text->Draw("same");
}

void drawLine(TLine*& line, Double_t width1, Double_t width2, Double_t height1, Double_t height2){
	line = new TLine(width1, height1, width2, height2 );
    line->SetLineStyle(lineStyle);
    line->SetLineColor(lineColor);
    line->SetLineWidth(lineWidth);
    line->Draw("same");     


}
void drawArrow(TArrow*& line, Double_t width1, Double_t width2, Double_t height1, Double_t height2){
	/*
	//TArrow does not have the option to enter NDC coordinates, this is conversion
    double x1_user = gPad->GetX1() + width1 * (gPad->GetX2() - gPad->GetX1());
    double y1_user = gPad->GetY1() + height1 * (gPad->GetY2() - gPad->GetY1());
    double x2_user = gPad->GetX1() + width2 * (gPad->GetX2() - gPad->GetX1());
    double y2_user = gPad->GetY1() + height2 * (gPad->GetY2() - gPad->GetY1());
    */

	line = new TArrow(width1, height1, width2, height2);
    line->SetLineStyle(lineStyle);
    line->SetLineColor(lineColor);
    line->SetLineWidth(lineWidth);
    //line->Draw("same");

}
/*
void clear(){
    if (c) delete c;
    if (text1) delete text1;
    if (text2) delete text2;
    if (text3) delete text3;
    if (textSTAR) delete textSTAR;
    if (line1) delete line1;
    if (line2) delete line2;
    if (line3) delete line3;
    if (legend) delete legend;

    c = nullptr;
    text1 = nullptr;
    text2 = nullptr;
    text3 = nullptr;
    textSTAR = nullptr;
    line1 = nullptr;
    line2 = nullptr;
    line3 = nullptr;
    legend = nullptr;
}
*/