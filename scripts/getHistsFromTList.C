#include <TFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TIterator.h>
#include <iostream>
#include <sstream> // For std::stringstream
#include <cstring>  // Include this for std::strcmp


using namespace std;


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

double xl = 0.6;
double yl = 0.85;
double xr = 0.9;
double yr = 0.85;

void getHistsFromTList(const char* filename) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the TList object from the file
    TList* list = dynamic_cast<TList*>(file->Get("HistList")); // Change "YourListName" to the actual name of the TList object
    if (!list) {
        std::cerr << "Error: TList not found in the file." << std::endl;
        file->Close();
        return;
    }

    // Create an iterator for the TList
    TIterator* it = list->MakeIterator();
    TObject* obj;
    int count = 0;

    while ((obj = it->Next())) {
        if (obj->IsA() == TH1D::Class() || obj->IsA() == TH2D::Class()) {
            // Create a new canvas for each histogram
            TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
            TLine *line;
            TArrow *arrow;


            canvas->cd(); // Make the canvas current
            gPad->SetMargin(0.13,0.03,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
            gPad->SetTickx();
            gPad->SetTicky(); 

             // setup the drawing style
            gStyle->SetFrameBorderMode(0);
            gStyle->SetFrameFillColor(0);
            gStyle->SetCanvasBorderMode(0);
            gStyle->SetCanvasColor(0);
            gStyle->SetPadBorderMode(0);
            gStyle->SetPadColor(0);
            gStyle->SetStatColor(0);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetNumberContours(99);
            gStyle->SetPalette(55);

            if (obj->IsA() == TH1D::Class()) {
                double xRel = 0.;
                double yRel = 0.;
                TH1D *hist;
                TH1D* histOriginal = dynamic_cast<TH1D*>(obj);
                if ((strcmp(histOriginal->GetName(), "hEta") == 0 || strcmp(histOriginal->GetName(), "hEtaAfter") == 0)){
                    hist = (TH1D*)histOriginal->Rebin(1, "");                    
                } else if( strcmp(histOriginal->GetName(), "hDcaDaughters") == 0 || strcmp(histOriginal->GetName(), "hDcaBeamline") == 0) {
                    hist = (TH1D*)histOriginal->Rebin(2, "");
                }

                hist->GetYaxis()->SetTitle("counts");
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
                //hist->GetYaxis()->SetRangeUser(0, 0.4); 
                hist->SetLineColor(kBlack);
                hist->SetLineStyle(1);
                hist->SetLineWidth(1);  
                hist->SetMarkerSize(2);
                hist->SetMarkerColor(kBlack);
                hist->SetMarkerStyle(markerStyle);

                hist->Draw("same E hist");

                //lebo nespravna podmienka
                if (strcmp(hist->GetName(), "hEtaAfter") == 0) {
                    for (int i = 1; i <= 4; ++i) {
                        hist->SetBinContent(i, 0);
                        hist->SetBinError(i, 0);
                    }
                }
                if (strcmp(hist->GetName(), "hEta") == 0 || strcmp(hist->GetName(), "hDecayLength") == 0) {
                    xRel = 0.06;
                    yRel = 0.03;
                }

                if(strcmp(hist->GetName(), "hDecayLength") == 0){

                    line = new TLine(3, 0, 3, 840e+6 );
                    line -> SetLineColor(kRed);
                    line -> SetLineStyle(kDashed);
                    line -> SetLineWidth(2);
                    line -> Draw("same");
                }


                if(strcmp(hist->GetName(), "hPointingAngle") == 0){
                    canvas->SetLogy();
                    xRel = 0.06;
                    yRel = 0.03;

                    line = new TLine(0.9, 0, 0.9, 5e+10 );
                    line -> SetLineColor(kRed);
                    line -> SetLineStyle(kDashed);
                    line -> SetLineWidth(4);
                    line -> Draw("same");

                    arrow = new TArrow(3, 4e+9, 1, 4e+9, 0.04, "|>");
                    arrow -> SetLineColor(kRed);
                    arrow -> SetLineStyle(kDashed);
                    arrow -> SetFillColor(kRed);
                    arrow -> SetAngle(40);
                    arrow -> SetLineWidth(4);
                    arrow -> Draw();

                }

                if(strcmp(hist->GetName(), "hDcaDaughters") == 0 ){
                    line = new TLine(3, 1600e+6, 3, 4e+9 );
                    line -> SetLineColor(kRed);
                    line -> SetLineStyle(kDashed);
                    line -> SetLineWidth(4);
                    line -> Draw("same");

                    arrow = new TArrow(3, 4e+9, 1, 4e+9, 0.04, "|>");
                    arrow -> SetLineColor(kRed);
                    arrow -> SetLineStyle(kDashed);
                    arrow -> SetFillColor(kRed);
                    arrow -> SetAngle(40);
                    arrow -> SetLineWidth(4);
                    arrow -> Draw();

                }
                if( strcmp(hist->GetName(), "hDcaBeamline") == 0){
                    line = new TLine(3, 1.8e+9, 3, 10e+9 );
                    line -> SetLineColor(kRed);
                    line -> SetLineStyle(kDashed);
                    line -> SetLineWidth(4);
                    line -> Draw("same");

                    arrow = new TArrow(3, 10e+9, 1.5, 10e+9, 0.04, "|>");
                    arrow -> SetLineColor(kRed);
                    arrow -> SetLineStyle(kDashed);
                    arrow -> SetFillColor(kRed);
                    arrow -> SetLineWidth(4);
                    arrow -> SetAngle(40);
                    arrow -> Draw();
                }

                TPaveText *textpp510;
                textpp510 = new TPaveText(xl+xRel, yl + yRel, xr+xRel, yr + yRel,"brNDC");
                textpp510 -> SetTextSize(textSize + 0.02);
                textpp510 -> SetTextAlign(11);
                textpp510 -> SetFillColor(0);
                textpp510 -> SetTextFont(62);
                textpp510->AddText(ppSTAR);
                textpp510 -> Draw("same");

                canvas->SaveAs(TString("preselection/") + hist->GetName() + TString("_presel.pdf"));

            } else if (obj->IsA() == TH2D::Class()) {
                TH2D* hist = dynamic_cast<TH2D*>(obj);
                gPad->SetMargin(0.11,0.12,0.11,0.06); // (Float_t left, Float_t right, Float_t bottom, Float_t top)

                if(strcmp(hist->GetName(), "hDecayLPointingA") == 0){
                    
                    hist->GetYaxis()->SetRangeUser(1,0);
                    //canvas->SetLogy();
                }else{
                    canvas->SetLogy();
                    canvas->SetLogx();
                }

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



                if (hist) {
                    hist->Draw("COLZ"); // Use "COLZ" to draw 2D histograms with color
                }else{
                    cout << "Couldn't get 2D histogram." << endl;
                }

                TPaveText *textpp510;
                textpp510 = new TPaveText(xl-0.05, yl, xr-0.05, yr,"brNDC");
                textpp510 -> SetTextSize(textSize + 0.02);
                textpp510 -> SetTextAlign(11);
                textpp510 -> SetFillColor(0);
                textpp510 -> SetTextFont(62);
                textpp510 -> SetTextColor(kWhite);
                textpp510->AddText(ppSTAR);
                textpp510 -> Draw("same");
                textpp510 -> SetFillStyle(3000);


                canvas->SaveAs(TString("preselection/") + hist->GetName() + TString("_presel.pdf"));

            }


            // Save the canvas as a PDF

            // Cleanup
            delete canvas;
        }
    }

    // Cleanup
    file->Close();
    delete file;
    delete it;
}
