#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>  // for std::pair
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include <iterator>
#include <string>

#include "TROOT.h"
#include <fstream>
#include <sstream>

// ----------------------------------------------------
// usage:
// root -l -b
// .L ProbRetainEvent.C+
// ProbRetainEvent("ProbRetainEvent.root")
// ----------------------------------------------------

using namespace std;
TFile *file;

struct runInfo {
    int runNumber;
    double luminosity;
    double probRetainEvent;
};

vector<runInfo> getData(TString inputList);
TGraph* fillGraph(vector<runInfo> data);
double fitRes(double luminosity, double A, double B) ;
vector<double> fitGaussian(TH1D *h, vector<runInfo> data);
void writeToFile(const string& filename, const vector<runInfo>& data, double mean, double sigma, double A, double B) ;

void ProbRetainEvent(TString filename ) {  // argument is the name of the output file
    // This function reads a ROOT file containing a TTree with event data
    // and creates a histogram of the probability of retaining an event
    
    vector<runInfo> data = getData("triggEff2.list");
    if (data.empty()) {
        cout << "No data found in the input list." << endl;
        return;
    }
    cout << "Obtained data from input list. Size: " << data.size() << endl;
    
    TGraph *graph = fillGraph(data);
    if (!graph) {
        cout << "Error filling graph." << endl;
        return;
    }
    
    cout << "Obtained graph. Size: " << graph->GetN() << endl;
    

    file = new TFile(filename, "recreate");
    if (!file->IsOpen()) {
        cout << "Error opening file: " << filename << endl;
        return;
    }
    file->cd();

    TCanvas *c = new TCanvas("c", "ProbRetainEvent", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    graph->SetTitle("");

    // define TGraph
    graph->SetName("ProbRetainEvent");
    graph->SetTitle("ProbRetainEvent");
    graph->GetXaxis()->SetTitle("L [#mu b^{-1}]");
    graph->GetYaxis()->SetTitle("probability of retaining an event");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    graph->SetLineStyle(1);

    // define an exponential function to fit the data
    TF1 *fit = new TF1("fit", "[0]*exp([1]*x)", 60, 160);
    fit->SetParameters(1, -0.01);
    fit->SetParNames("A", "B");
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);
    //fit->SetLineStyle(2);
    graph->Fit(fit, "R");
    graph->GetYaxis()->SetRangeUser(0, 1);

    graph->Draw("AEP");
    fit->Draw("same");



    TLegend *legend = new TLegend(0.12, 0.7, 0.3, 0.88);
    legend->AddEntry(graph, "Data", "lp");
    legend->AddEntry(fit, "Exponential fit", "lp");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(0.04);
    legend->Draw("same");

    TPaveText *textpp510 = new TPaveText(0.77, 0.77, 0.94, 0.96,"brNDC");
    textpp510 -> SetTextSize(0.04);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillStyle(4000);
    textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #rightarrow p #oplus J/#psi #oplus p");
    textpp510->AddText( "#sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");

    TPaveText *text = new TPaveText(0.77, 0.6, 0.94, 0.75,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetTextAlign(33);
    text -> SetFillStyle(4000);
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetTextFont(82);
    // print chi square and NDF
    double chi2 = fit->GetChisquare();
    double ndf = fit->GetNDF();
    text -> AddText(Form("#chi^{2}/ndf = %.2f/%.0f", chi2, ndf));
    text -> AddText("f(L) = A*exp(B*L)");
    text -> AddText(Form("A = %.3f #pm %.3f", fit->GetParameter(0), fit->GetParError(0)));
    text -> AddText(Form("B = %.3f #pm %.3f", fit->GetParameter(1), fit->GetParError(1)));
    text -> Draw("same");


    TH1D *h = new TH1D("h", "", 100, -0.5,0.5 );
    // loop over all points in the graph 
    for (int i = 0; i < graph->GetN(); i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        double funcVal = fitRes(x, fit->GetParameter(0), fit->GetParameter(1));
        // difference between y and function value
        h->Fill(y-funcVal);   
    }

    vector<double> fitPars= fitGaussian(h,data);
    if(fitPars.size() == 0){
        cout << "Error fitting Gaussian." << endl;
        return;
    }

    writeToFile("goodRun17.list", data, fitPars[0], fitPars[1], fit->GetParameter(0), fit->GetParameter(1));

    cout << "Data written to file: goodRun17.list" << endl;


    cout << "back in main function" << endl;
    c->Write("ProbRetainEvent");
    cout << "ProbRetainEvent written to file" << endl;
    c->SaveAs("ProbRetainEvent.pdf");
    c->SaveAs("ProbRetainEvent.png");
    cout << "ProbRetainEvent saved as pdf" << endl;
    file->Close();
    cout << "ProbRetainEvent file closed" << endl;
    
}

double fitRes(double luminosity, double A, double B) {
    // Calculate the probability of retaining an event based on luminosity
    return A * exp(B * luminosity);
}

vector<double> fitGaussian(TH1D *h, vector<runInfo> data) {
    // Fit a Gaussian to the histogram
    TCanvas *c1 = new TCanvas("c1", "ProbRetainEvent", 1500, 1200);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetMargin(0.105,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    c1->SetLogy();

    if(!h || h->GetEntries() == 0) {
        cout << "Error: histogram is null." << endl;
        vector<double> empty;
        return empty;
    }

    TF1 *gaus = new TF1("gaus", "gaus", -0.5, 0.5);
    gaus->SetParameters(500, 0, 0.1); // initial parameters for the Gaussian
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    h->SetTitle("");
    h->Fit(gaus, "R");

    h->SetLineColor(kBlue);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1);
    h->SetMarkerColor(kBlue);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("data - fit");
    h->GetYaxis()->SetTitle("counts");
    h->Draw("EP");

    gaus->Draw("same");
    
    double mean = gaus->GetParameter(1);
    double sigma = gaus->GetParameter(2);

    // create 2 TLine objects plus and minus 9 sigma
    TLine *linePlus = new TLine(mean + 9 * sigma, 0, mean + 9 * sigma, 25);
    linePlus->SetLineColor(kBlack);
    linePlus->SetLineStyle(2);
    linePlus->SetLineWidth(2);
    linePlus->Draw("same");
    TLine *lineMinus = new TLine(mean - 9 * sigma, 0, mean - 9 * sigma, 25);
    lineMinus->SetLineColor(kBlack);
    lineMinus->SetLineStyle(2);
    lineMinus->SetLineWidth(2);
    lineMinus->Draw("same");

    TLegend *legend = new TLegend(0.18 , 0.7, 0.35, 0.88);
    legend->AddEntry(h, "Data - Fit", "lp");
    legend->AddEntry(gaus, "Gaussian fit", "lp");
    legend->AddEntry(linePlus, "Mean #pm 9#sigma", "l");
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(0.04);
    legend->Draw("same");

    TPaveText *textpp510 = new TPaveText(0.77, 0.77, 0.94, 0.96,"brNDC");
    textpp510 -> SetTextSize(0.04);
    textpp510 -> SetTextAlign(33);
    textpp510 -> SetFillStyle(4000);
    textpp510 -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    textpp510 -> SetTextFont(62);
    textpp510->AddText("p+p #rightarrow p #oplus J/#psi #oplus p");
    textpp510->AddText( "#sqrt{s} = 510 GeV");
    textpp510 -> Draw("same");

    TPaveText *text = new TPaveText(0.77, 0.6, 0.94, 0.75,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetTextAlign(33);
    text -> SetFillStyle(4000);
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetFillColorAlpha(kWhite, 0.0);  // Blue background with 30% opacity
    text -> SetTextFont(82);
    // print chi square and NDF
    double chi2 = gaus->GetChisquare();
    double ndf = gaus->GetNDF();
    text -> AddText(Form("#chi^{2}/ndf = %.2f/%.0f", chi2, ndf));
    text -> AddText(Form("#mu = %.3f #pm %.3f", gaus->GetParameter(1), gaus->GetParError(1)));
    text -> AddText(Form("#sigma = %.3f #pm %.3f", gaus->GetParameter(2), gaus->GetParError(2)));
    text -> Draw("same");

    file->cd();
    c1->Write("ProbRetainEventGaussianFit");
    c1->SaveAs("ProbRetainEventGaussianFit.pdf");
    c1->SaveAs("ProbRetainEventGaussianFit.png");
    c1->Close();

    cout << "Gaussian fit done." << endl;
    vector<double> fitParams;
    fitParams.push_back(mean);
    fitParams.push_back(sigma);
    return fitParams;

}

void writeToFile(const string& filename, const vector<runInfo>& data, double mean, double sigma, double A, double B) {
    ofstream file(filename);
    if (!file.is_open()) {
        cout << "Error opening file: " << filename << endl;
        return;
    }
    cout << "Size of all data: " << data.size() << endl;
    int i = 0;
    for (const auto& entry : data) {

        double theoreticalProb = fitRes(entry.luminosity, A, B);
        double diff = entry.probRetainEvent - theoreticalProb;
        if(diff > mean + 9 * sigma || diff < mean - 9 * sigma){
            //cout << "Error: probRetainEvent out of range in file: " << filename << endl;
            continue;
        }

        file << entry.runNumber << " " << entry.luminosity << " " << entry.probRetainEvent << endl;
        i++;
    }
    file.close();

    cout << "Number of entries written to file: " << i << endl;
    cout << "Data written to file: " << filename << endl;
    return;
}



vector<runInfo> getData(TString inputList){

    //open .list file with name inputList
    ifstream file;
    file.open(inputList);
    if (!file.is_open()) {
        cout << "Error opening file: " << inputList << endl;
        return {};
    }

    cout << "Reading input list: " << inputList << endl;


    vector<runInfo> data;
    string line;
    while(getline(file, line)) {
        if(line.empty()) continue; // skip empty lines

        TFile *File = new TFile(line.c_str(), "read");
        if (!File->IsOpen()) {
            cout << "Error opening file: " << line << endl;
            continue;
        }
        // find histogram which starts with "hRunNumber_" in all directories


        TDirectory *dir = (TDirectory*)File->Get("TofQA");
        if (!dir) { 
            cout << "Error: TofQA directory not found in file: " << line << endl;
            File->Close();
            continue;
         }
        
        TList *list = dir->GetListOfKeys();  
        if (!list || list->IsEmpty()) {
            //cout << "Warning: TofQA directory is empty in file: " << line << endl;
            continue;
        }
        int runNumber = 0;

        // load the only histogram in the directory
        TKey *key = (TKey*)list->First();
        if (!key) {
            //cout << "Error: no keys found in TofQA directory in file: " << line << endl;
            File->Close();
            continue;
        }
        TObject *obj = key->ReadObj();
        if (!obj) {
            //cout << "Error: object not found in TofQA directory in file: " << line << endl;
            File->Close();
            continue;
        }
        TString name = obj->GetName();
        TObjArray *tokens = name.Tokenize("_");
        if (tokens->GetEntries() == 2) {
            TString runStr = ((TObjString*)tokens->At(1))->GetString();
            runNumber = runStr.Atoi();
        }
        delete tokens;
        //cout << "Run number: " << runNumber << endl;
        

        
        double luminosity, probRetainEvent;
        TGraph *graph = (TGraph*)File->Get("myGraph");
        if(graph){
            probRetainEvent = graph->GetY()[0];
            luminosity = graph->GetX()[0];
            runInfo info;
            info.runNumber = runNumber;
            info.luminosity = luminosity;
            info.probRetainEvent = probRetainEvent;
            if(probRetainEvent <= 0){
                //cout << "Error: probRetainEvent <= 0 in file: " << line << endl;
                File->Close();
                continue;
            }

            data.push_back(info);
        }else{
            cout << "Error: graph not found in file: " << line << endl;
        }
        File->Close();
    }
    
    file.close();
    return data;

}

TGraph* fillGraph(vector<runInfo> data){
    
    TGraph *graph = new TGraph();

    std::ofstream txtFile;
    txtFile.open("ProbRetainEvent.txt");
    if (!txtFile.is_open()) {
        cout << "Error opening file: ProbRetainEvent.txt" << endl;
        return nullptr;
    }

    for (const auto& entry : data) {
        int runNumber = entry.runNumber;
        double luminosity = entry.luminosity;
        double probRetainEvent = entry.probRetainEvent;
        graph->SetPoint(graph->GetN(), luminosity, probRetainEvent);

        txtFile << runNumber << " " << luminosity << " " << probRetainEvent << endl;

    }
    txtFile.close();
    return graph;
}

