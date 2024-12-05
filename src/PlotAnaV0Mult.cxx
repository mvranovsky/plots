#include "../include/PlotAnaV0Mult.h"


PlotAnaV0Mult::PlotAnaV0Mult(TFile *mOutFile, const string mInputList, const char* filePath): Plot(mOutFile, mInputList, filePath){}

void PlotAnaV0Mult::Make(){


}

void PlotAnaV0Mult::Init(){
	//define the output file which will store all the canvases
	outFile = new TFile(outputPosition, "recreate");


	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaV0MultTree);

    if(!tree){
    	cerr << "Couldn't open tree with data. Returning." << endl;
    	return;
    }
}

void PlotAnaV0Mult::general2f(TString nameOfHist , TString cutDescription){
	//TCanvas* canvas;
	CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );
	SetGPad(false,0.12, 0.16,0.11, 0.06 );
	TH2F* hist = (TH2F*)inFile->Get(nameOfHist);
	if (!hist){
		cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
		return;
	}
	canvas->Clear();
	SetTH2Style(hist);
	hist->Draw("COLZ");
	DrawSTARpp510();
	CreateText(cutDescription, 0.1, 0.8,0.4,0.8);

	outFile->cd();
	canvas->Write();
	hist->Write();
}


void PlotAnaV0Mult::nSigmaCorr(){

	for (int i = 0; i < 3; ++i){
		CreateCanvas(&canvas,plots[2+i].first, widthTypical, heightTypical);
		SetGPad();
		TH2F* hist = (TH2F*)inFile->Get(plots[2+i].first);
		if(!hist){
			cerr << "Could not open histogram "<< plots[2+i].first << " from inFile."<< endl;
			return;
		}

		canvas->Clear();
		SetTH2Style(hist);
		hist->Draw("COLZ");
		DrawSTARpp510();


    	outFile->cd();
    	canvas->Write();
    	hist->Write();
	}
}

void PlotAnaV0Mult::histGeneral(TString nameOfHist, TString cutDescription) {

	//TCanvas* canvas;
	//cout << "Creating " << nameOfHist << " in histGeneral." << endl;
	CreateCanvas(&canvas, nameOfHist, widthTypical, heightTypical );

	TH1D* hist = (TH1D*)inFile->Get(nameOfHist);
	if (!hist){
		cerr << "Could not open histogram "<< nameOfHist << " from inFile."<< endl;
		return;
	}
	canvas->Clear();
	SetHistStyle(hist, kBlack, markerStyleTypical);


	hist->Draw("same");
	if(nameOfHist == "hAnalysisFlow"){
		SetGPad(true, 0.14, 0.05,0.11,0.06);
		DrawSTARpp510(0.6,0.9,0.9,0.9);
	}else{
		SetGPad(false, 0.14, 0.05,0.11,0.06);
		DrawSTARpp510();
	}

	//CreateText(cutDescription, 0.1, 0.8,0.4,0.8);

	outFile->cd();
	canvas->Write();
	hist->Write();
}


void PlotAnaV0Mult::invMassLambda(Double_t minRange, Double_t maxRange, int nBins){

   	//TCanvas* canvas;
	CreateCanvas(&canvas,"invMassLambda", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassLambdaSignal", "invMassLambdaSignal", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{p #pi} [GeV/c^{2}]");
    TH1D *bcgFinal = new TH1D("invMassLambdaBcg", "invMassLambdaBcg", nBins, minRange, maxRange);
    bcgFinal->GetXaxis()->SetTitle("m_{p #pi} [GeV/c^{2}]");

    
    TString cmd, condition;
    for (int i = 0; i < nStates; ++i){
    	//draw the pairs invMasses which satisfy condition
    	cmd = TString::Format("invMass%d>>hist(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("(pairID%d == 1 || pairID%d ==2) && totQ%d == 0", i, i, i);

    	tree->Draw(cmd, condition);
    	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

	    cmd = TString::Format("invMass%d>>bcg(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("(pairID%d == 1 || pairID%d == 2) && totQ%d != 0", i, i,i);
    	tree->Draw(cmd, condition);
    	bcgFinal->Add((TH1D*)gPad->GetPrimitive(TString("bcg")));
    }




	canvas->Clear();

	//save histograms to the outfile so you can use roofit library for playing


    SetHistStyle(signalFinal, kBlack, markerStyleTypical);
    SetHistStyle(bcgFinal, kRed, markerStyleTypical);
    signalFinal->Draw();
    bcgFinal->Draw("same");
    DrawSTARpp510();

    CreateLegend(&legend, 0.45, 0.9, 0.7, 0.75);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    legend->AddEntry(bcgFinal, "Data (like-sign pairs)","ple");
    legend->Draw("same");

    outFile->cd();
    canvas->Write();
    signalFinal->Write();
    bcgFinal->Write();

}



void PlotAnaV0Mult::invMassK0s(Double_t minRange, Double_t maxRange, int nBins){
	//TCanvas* canvas;
	CreateCanvas(&canvas,"invMassK0s", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassK0sSignal", "invMassK0sSignal", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{#pi #pi} [GeV/c^{2}]");
    TH1D *bcgFinal = new TH1D("invMassK0sBcg", "invMassK0sBcg", nBins, minRange, maxRange);
    bcgFinal->GetXaxis()->SetTitle("m_{#pi #pi} [GeV/c^{2}]");

    cout << "Creating histograms for K0S" << endl;
    TString cmd, condition;
    for (int i = 0; i < nStates; ++i){
    	//draw the pairs invMasses which satisfy condition
    	cmd = TString::Format("invMass%d>>hist(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("pairID%d == 0 && totQ%d == 0", i, i);
    	cout << "Command: " << cmd << endl;
    	cout << "Condition: " << condition << endl;
    	if (!tree){
    		cerr << "Couldn't open tree in PlotAnaV0Mult." << endl;
    		return;
    	}
    	tree->Draw(cmd, condition);
    	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
    	cmd = TString::Format("invMass%d>>bcg(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("pairID%d == 0 && totQ%d != 0", i, i);
	    tree->Draw(cmd, condition);
    	bcgFinal->Add((TH1D*)gPad->GetPrimitive(TString("bcg")));
    }
    cout << "Histograms for K0S created" << endl;

	canvas->Clear();


    SetHistStyle(signalFinal, kBlack, markerStyleTypical);
    SetHistStyle(bcgFinal, kRed, markerStyleTypical);
    signalFinal->Draw();
    bcgFinal->Draw("same");
    DrawSTARpp510();

    CreateLegend(&legend, 0.45, 0.9, 0.7, 0.75);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    legend->AddEntry(bcgFinal, "Data (like-sign pairs)","ple");
    legend->Draw("same");

    outFile->cd();
    canvas->Write();
    signalFinal->Write();
    bcgFinal->Write();
	
}

void PlotAnaV0Mult::fitK0s(Double_t minRange, Double_t maxRange, int nBins) {
	CreateCanvas(&canvas,"fitK0s", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("fitPeakK0sSignalFit", "fitPeakK0sSignalFit", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    
    TString cmd, condition;
    cout << "fitting peak for K0S" << endl;
    for (int i = 0; i < nStates; ++i){
    	//draw the pairs invMasses which satisfy condition
    	cmd = TString::Format("invMass%d>>hist(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("pairID%d == 0 && totQ%d == 0", i, i);
    	tree->Draw(cmd, condition);
	    signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
	}
    cout << "peak fitted for K0S" << endl;

	canvas->Clear();

	SetHistStyle(signalFinal, kBlack, markerStyleTypical);
	signalFinal->SetMarkerStyle(markerStyleTypical);
    signalFinal->Draw();
	int y = fitGaussPol2(&signalFinal,(maxRange-minRange)/nBins, minRange, maxRange);
    CreateText("Yield: " + to_string(y), 0.7, 0.85,0.85,0.85);

    DrawSTARpp510();
    CreateLegend(&legend,0.15, 0.85, 0.4, 0.7);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    //legend->AddEntry(func, "Fit signal + bcg","ple");
    legend->Draw("same");

    outFile->cd();
    canvas->Write();
    signalFinal->Write();
}

void PlotAnaV0Mult::fitLambda(Double_t minRange, Double_t maxRange, int nBins) {
	CreateCanvas(&canvas,"fitLambda", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassLambdaSignalFit", "invMassLambdaSignalFit", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{p #pi^{+}} [GeV/c^{2}]");
    
    TString cmd, condition;
    for (int i = 0; i < nStates; ++i){
    	//draw the pairs invMasses which satisfy condition
    	cmd = TString::Format("invMass%d>>hist(%d, %f, %f)", i, nBins, minRange,maxRange);
		condition = TString::Format("(pairID%d == 1 || pairID%d == 2) && totQ%d == 0", i,i,i);
    	tree->Draw(cmd, condition);
	    signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
	}
	canvas->Clear();

	SetHistStyle(signalFinal, kBlack, markerStyleTypical);
	signalFinal->GetYaxis()->SetTitleOffset(1.5);
    signalFinal->Draw();
	int y = fitGaussPol2(&signalFinal,(maxRange-minRange)/nBins, minRange, maxRange, 0.01, 0.01, 0.01, 5000, 1.15, 0.01);
    CreateText("Yield: " + to_string(y), 0.15, 0.8,0.4,0.8);

    DrawSTARpp510();
    CreateLegend(&legend,0.15, 0.7, 0.5, 0.6);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    //legend->AddEntry(func, "Fit signal + bcg","ple");
    legend->Draw("same");

    outFile->cd();
    canvas->Write();
    signalFinal->Write();
}
