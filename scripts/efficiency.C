#include <iostream>
#include <cmath>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <RooFitResult.h>


//using namespace RooFit;
using namespace std;


Double_t etaVal1[9] = {12427, 16518, 20031,18410,1, 20507, 19177, 23420, 14435};
Double_t etaErr1[9] = {3348, 509, 1109,827,0, 3399, 820, 730, 575};
Double_t etaVal2[9] = {6003, 10194, 11810,13533,0, 14734, 14418, 10901, 6963};
Double_t etaErr2[9] = {194, 1222,230,263,0, 390, 388, 319, 139};

Double_t etaAltVal1[9] = {10893,17494, 16042, 16222,12796,15897, 15929,19459,11294};
Double_t etaAltErr1[9] = {1810, 4908, 557, 280,415 ,562, 521, 2729, 425};
Double_t etaAltVal2[9] = {4346, 7895, 9317, 11109,11625, 11511, 10988, 8274, 5002};
Double_t etaAltErr2[9] = {871, 383, 160, 195, 1084,196, 202, 201, 152};


Double_t phiVal1[6] = {23120,20783, 20867, 24829,24533, 25811 };
Double_t phiErr1[6] = {728, 910, 5959, 1358, 1056, 1279};
Double_t phiVal2[6] = {14241, 14675, 14737, 15821, 17489, 15953};
Double_t phiErr2[6] = {278, 880, 1118, 492, 1383, 261};

Double_t phiAltVal1[6] = {23518,23003,22622, 23563, 23883, 22387};
Double_t phiAltErr1[6] = {3302, 912, 2633, 701, 4026, 3357};
Double_t phiAltVal2[6] = {11885, 13142, 14323, 14919, 16873, 14921};
Double_t phiAltErr2[6] = {190, 244, 827, 532, 326, 485};


Double_t vzVal1[10] = {5065, 7823, 11708,16708, 21898, 24952, 19974, 13317,8753,5375 };
Double_t vzErr1[10] = {364, 445, 1898, 705, 818, 4084, 827, 648, 436, 333};
Double_t vzVal2[10] = {2752, 4644, 7901, 11617, 16288, 16367, 12326, 8950, 5473, 2962};
Double_t vzErr2[10] = {105, 197, 281, 242, 312, 1176, 329, 802, 156, 110};

Double_t vzAltVal1[10] = {5317, 9226,11913, 15552,20077, 22087,19454, 13558, 9287,5426 };
Double_t vzAltErr1[10] = {264, 2935, 3832, 471, 769, 797, 2779, 2908, 1328, 231};
Double_t vzAltVal2[10] = {2863, 4394, 7498, 9915, 13108, 13638, 11281, 8179, 5726, 3449};
Double_t vzAltErr2[10] = {94, 136, 413, 216, 211, 245, 181, 254, 184, 100};


Double_t ptVal1[6] = {83295, 32477, 12631, 4982,2204, 1145};
Double_t ptErr1[6] = {1681, 505, 1084, 2352,206, 180};
Double_t ptVal2[6] = { 57505, 21988, 7366, 2881,1085,408 };
Double_t ptErr2[6] = {1688, 1984, 166, 127,383, 43};

Double_t ptVal1Bins[5] = {83295, 32477, 12631, 4982, 3105};
Double_t ptErr1Bins[5] = {1681, 505, 1084, 2352, 226};
Double_t ptVal2Bins[5] = { 57505, 21988, 7366, 2881, 1563};
Double_t ptErr2Bins[5] = {1688, 1984, 166, 127, 112};


Double_t ptAltVal1[6] = {82820, 32087,15261, 5893, 2426, 1153};
Double_t ptAltErr1[6] = {1053, 778, 2835, 1586, 722, 146};
Double_t ptAltVal2[6] = {50033, 19369, 6411, 2530, 967, 400};
Double_t ptAltErr2[6] = {390,292, 176, 821, 72, 39};

Double_t ptAltVal1Bins[5] = {82820, 32087,15261, 5893, 3082};
Double_t ptAltErr1Bins[5] = {1053, 778, 2835, 1586, 183};
Double_t ptAltVal2Bins[5] = { 50033, 19369, 6411, 2530, 1414};
Double_t ptAltErr2Bins[5] = {390,292, 176, 821, 80};


//TString discriminator = "Tag";
TString discriminator = "Probe";
TString tag;

void plot(int runSwitch, int nBins, Double_t lowRange, Double_t topRange, TString xAxisDescription, bool useMean);
void plot2D(TH2F* hist1, TH2F* hist2, TString xAxisDescription, TString yAxisDescription);
vector<pair<double,double>> fillVector1(int runSwitch);
vector<pair<double,double>> fillVector2(int runSwitch);
TString getName(int runSwitch);
void PlotPtSpecialBin();


void efficiency(const string inputTag) {

	cout << "Beginning effectivity plots." << endl;

    tag = inputTag;

    TString inputFile = inputTag + "/AnalysisOutput.root";

    /*
	//file1 means 1 ToF track, file2 means 2 ToF tracks
   	TFile* file = new TFile(inputFile);
   	if(!file){
   		cerr << "Couldn't open the input file. Leaving." << endl;
   		return;
   	}
   	cout << "Opened input file." << endl;
    */


   	//pseudorapidity plot
   	plot(1, 9, -0.9, 0.9 , "#eta [-]", false);

   	//pT plot
   	plot(3, 6, 0.2, 1.4 , "p_{T} [GeV/c]", false);
    PlotPtSpecialBin();

    //phi plot
    plot(2, 6, -180, 180, "#phi [-]", false);

    //Vz plot
    plot(4, 10, -100, 100, "V_{Z} [cm]", false);

    //eta-phi plot
    //plot2D((TH2F*)file->Get("hEtaPhi1" + discriminator), (TH2F*)file->Get("hEtaPhi2"+discriminator), "#eta [-]", "#phi [-]");

   	cout << "All plots created, Goodbye."<< endl;
    //file->Close(); // Close the ROOT file

}

vector<pair<double,double>> fillVector1(int isAlternate ,int runSwitch){
    vector<pair<double,double>> result;

    if(!isAlternate){
        if (runSwitch == 1){ //eta
            for (int i = 0; i < 9; ++i){
                result.push_back( make_pair( etaVal1[i], etaErr1[i] ) );        
            }
        } else if(runSwitch == 2){ //phi
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( phiVal1[i], phiErr1[i] ) );
            }
        } else if(runSwitch == 3){ //pT
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( ptVal1[i], ptErr1[i] ) );   
            }
        } else if(runSwitch == 4){ //Vz
            for (int i = 0; i < 10; ++i){
                result.push_back( make_pair( vzVal1[i], vzErr1[i] ) );   
            }
        } else {
            cerr << "Incorrect runSwitch. Returning 0 vector..." << endl;
        }
    } else {
        if (runSwitch == 1){ //eta
            for (int i = 0; i < 9; ++i){
                result.push_back( make_pair( etaAltVal1[i], etaAltErr1[i] ) );        
            }
        } else if(runSwitch == 2){ //phi
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( phiAltVal1[i], phiAltErr1[i] ) );
            }
        } else if(runSwitch == 3){ //pT
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( ptAltVal1[i], ptAltErr1[i] ) );   
            }
        } else if(runSwitch == 4){ //Vz
            for (int i = 0; i < 10; ++i){
                result.push_back( make_pair( vzAltVal1[i], vzAltErr1[i] ) );   
            }
        } else {
            cerr << "Incorrect runSwitch. Returning 0 vector..." << endl;
        }
    }

    return result;
}


vector<pair<double,double>> fillVector2(int isAlternate,int runSwitch){
    vector<pair<double,double>> result;

    if(!isAlternate){
        if (runSwitch == 1){ //eta
            for (int i = 0; i < 9; ++i){
                result.push_back( make_pair( etaVal2[i], etaErr2[i] ) );        
            }
        } else if(runSwitch == 2){ //phi
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( phiVal2[i], phiErr2[i] ) );
            }
        } else if(runSwitch == 3){ //pT
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( ptVal2[i], ptErr2[i] ) );   
            }
        } else if(runSwitch == 4){ //Vz
            for (int i = 0; i < 10; ++i){
                result.push_back( make_pair( vzVal2[i], vzErr2[i] ) );   
            }
        } else {
            cerr << "Incorrect runSwitch. Returning 0 vector..." << endl;
        }
    } else {
        if (runSwitch == 1){ //eta
            for (int i = 0; i < 9; ++i){
                result.push_back( make_pair( etaAltVal2[i], etaAltErr2[i] ) );        
            }
        } else if(runSwitch == 2){ //phi
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( phiAltVal2[i], phiAltErr2[i] ) );
            }
        } else if(runSwitch == 3){ //pT
            for (int i = 0; i < 6; ++i){
                result.push_back( make_pair( ptAltVal2[i], ptAltErr2[i] ) );   
            }
        } else if(runSwitch == 4){ //Vz
            for (int i = 0; i < 10; ++i){
                result.push_back( make_pair( vzAltVal2[i], vzAltErr2[i] ) );   
            }
        } else {
            cerr << "Incorrect runSwitch. Returning 0 vector..." << endl;
        }
    }

    return result;
}


void plot(int runSwitch, int nBins, Double_t lowRange, Double_t topRange, TString xAxisDescription, bool useMean){


    vector<pair<double, double>> vector1, vector2, vector1Alt, vector2Alt;

    vector1 = fillVector1(false,runSwitch);
    vector1Alt = fillVector1(true, runSwitch);
    vector2 = fillVector2(false,runSwitch);
    vector2Alt = fillVector2(true,runSwitch);

    if(vector1.size() != vector2.size() || vector1Alt.size() != vector2Alt.size()){
        cerr << "Vector sizes are not the same. Leaving." << endl;
        return;
    }

	TCanvas *c = new TCanvas("c", "Fit Result", 800, 600);

    TH1D *hist1 = new TH1D("hist1","hist1", nBins, lowRange, topRange);
    TH1D *hist2 = new TH1D("hist2","hist2", nBins, lowRange, topRange);

    TH1D *hist1Alt = new TH1D("hist1Alt","hist1", nBins, lowRange, topRange);
    TH1D *hist2Alt = new TH1D("hist2Alt","hist2", nBins, lowRange, topRange);

    hist1->Sumw2();
    hist2->Sumw2();
    hist1Alt->Sumw2();
    hist2Alt->Sumw2();

    for (int iBin = 1; iBin <= nBins ; ++iBin){

        hist2Alt->SetBinContent(iBin, vector2Alt[iBin-1].first);
        hist2Alt->SetBinError(iBin, vector2Alt[iBin-1].second);
        hist1Alt->SetBinContent(iBin, vector1Alt[iBin-1].first);
        hist1Alt->SetBinError(iBin, vector1Alt[iBin-1].second);

        if(iBin == 4 && runSwitch == 1)
            continue;
        hist2->SetBinContent(iBin, vector2[iBin-1].first);
        hist2->SetBinError(iBin, vector2[iBin-1].second);
        hist1->SetBinContent(iBin, vector1[iBin-1].first);
        hist1->SetBinError(iBin, vector1[iBin-1].second);

    }



    TEfficiency* pEff = new TEfficiency(*hist2,*hist1);
    TGraphAsymmErrors *gr = pEff->CreateGraph();
    gr->GetYaxis()->SetTitle("efficiency");
    gr->GetXaxis()->SetTitle(xAxisDescription);
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed); 
    gr->GetYaxis()->SetRangeUser(0,1);
    gr->Draw("AP");

    TEfficiency* pEffAlt = new TEfficiency(*hist2Alt,*hist1Alt);
    TGraphAsymmErrors *grAlt = pEffAlt->CreateGraph();
    grAlt->GetYaxis()->SetTitle("efficiency");
    grAlt->GetXaxis()->SetTitle(xAxisDescription);
    grAlt->SetMarkerStyle(21);
    grAlt->SetMarkerColor(kBlue);
    grAlt->SetLineColor(kBlue); 
    grAlt->GetYaxis()->SetRangeUser(0,1);
    grAlt->Draw("same P");

	TPaveText* text2 = new TPaveText(0.15, 0.85 ,0.4, 0.85 ,"brNDC");
    text2 -> SetTextSize(0.04);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");



    TLegend *leg = new TLegend(0.55, 0.9, 0.75, 0.8);
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(0.03);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(grAlt," efficiency K^{0}_{S} (old upcDst)", "lep");
    leg -> AddEntry(gr," efficiency K^{0}_{S} (extended upcDst)", "lep");
    leg -> Draw("same");

    TString name = getName(runSwitch);

    c->SaveAs("eff_" + tag + "/eff" + name + ".pdf");  // Save the canvas to a file

}

TString getName(int runSwitch){
    TString name;
    if (runSwitch == 1){
        name = "eta";
    } else if(runSwitch == 2){
        name = "phi";
    } else if(runSwitch == 3){
        name = "pT";
    } else if(runSwitch == 4){
        name = "Vz";
    } else {
        cerr << "Incorrect runSwitch, leaving..." << endl;
    }
    return name;
}


void PlotPtSpecialBin(){

    double binEdges[] = {0.2, 0.4, 0.6, 0.8, 1., 1.4};
    int numBins = sizeof(binEdges)/sizeof(double) - 1;  // number of bins is one less than number of edges


    TCanvas *c = new TCanvas("c", "Fit Result", 800, 600);

    TH1D *hist1 = new TH1D("hist1","hist1", numBins,binEdges);
    TH1D *hist2 = new TH1D("hist2","hist2", numBins, binEdges);

    TH1D *hist1Alt = new TH1D("hist1Alt","hist1", numBins,binEdges);
    TH1D *hist2Alt = new TH1D("hist2Alt","hist2", numBins, binEdges);

    hist1->Sumw2();
    hist2->Sumw2();

    hist1Alt->Sumw2();
    hist2Alt->Sumw2();

    for (int iBin = 1; iBin <= numBins ; ++iBin){
        hist1Alt->SetBinContent(iBin, ptAltVal1Bins[iBin-1]);
        hist1Alt->SetBinError(iBin, ptAltErr1Bins[iBin-1]);
        hist2Alt->SetBinContent(iBin, ptAltVal2Bins[iBin-1]);
        hist2Alt->SetBinError(iBin, ptAltErr2Bins[iBin-1]);

        hist2->SetBinContent(iBin, ptVal2Bins[iBin-1]);
        hist2->SetBinError(iBin, ptErr2Bins[iBin-1]);
        hist1->SetBinContent(iBin, ptVal1Bins[iBin-1]);
        hist1->SetBinError(iBin, ptErr1Bins[iBin-1]);
    }


    TEfficiency* pEff = new TEfficiency(*hist2,*hist1);
    TGraphAsymmErrors *gr = pEff->CreateGraph();
    gr->GetYaxis()->SetTitle("efficiency");
    gr->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed); 
    gr->GetYaxis()->SetRangeUser(0,1);
    gr->Draw("AP");

    TEfficiency* pEffAlt = new TEfficiency(*hist2Alt,*hist1Alt);
    TGraphAsymmErrors *grAlt = pEffAlt->CreateGraph();
    grAlt->GetYaxis()->SetTitle("efficiency");
    grAlt->GetXaxis()->SetTitle("p_{T} [GeV]");
    grAlt->SetMarkerStyle(21);
    grAlt->SetMarkerColor(kBlue);
    grAlt->SetLineColor(kBlue); 
    grAlt->GetYaxis()->SetRangeUser(0,1);
    grAlt->Draw("same P");

    TPaveText* text2 = new TPaveText(0.15, 0.85 ,0.4, 0.85 ,"brNDC");
    text2 -> SetTextSize(0.04);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");


    TLegend *leg = new TLegend(0.55, 0.9, 0.75, 0.8);
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(0.03);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(grAlt," efficiency K^{0}_{S} (old upcDst)", "lep");
    leg -> AddEntry(gr," efficiency K^{0}_{S} (extended upcDst)", "lep");
    leg -> Draw("same");


    c->SaveAs("eff_" + tag + "/effPtBins.pdf");  // Save the canvas to a file

}





/*
void plot2D(TH2F* hist1, TH2F* hist2, TString xAxisDescription, TString yAxisDescription){


    TCanvas *c = new TCanvas("c", "Fit Result", 800, 600);

    //dividing histograms? GetBinCenter() alebo mame specificku funkciu na to?



    hist2->GetZaxis()->SetTitle("efficiency");
    hist2->GetXaxis()->SetTitle(xAxisDescription);
    hist2->GetYaxis()->SetTitle(yAxisDescription);

    hist2->Divide(hist1);
    hist2->Draw("same E");


    TPaveText* text2 = new TPaveText(0.6, 0.85, 0.85, 0.85,"brNDC");
    text2 -> SetTextSize(0.04);
    text2 -> SetFillColor(0);
    text2 -> SetTextFont(62);
    text2 -> AddText("pp #sqrt{s} = 510 GeV");
    text2 -> Draw("same");



    TLegend *leg = new TLegend(0.7, 0.8 ,0.9, 0.55 );
    leg -> SetTextAlign(11);
    leg -> SetFillStyle(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(0.05);
    leg -> SetTextFont(62);
    leg -> SetMargin(0.1);   
    leg -> AddEntry(hist2, discriminator + " effiency", "lep");
    leg -> Draw("same");





    c->SaveAs(tag + "/efficiency" + hist2->GetName() + ".pdf");  // Save the canvas to a file




}
*/