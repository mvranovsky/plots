#include "../include/PlotAnaV0Mult.h"


PlotAnaV0Mult::PlotAnaV0Mult(const string mInputList, const char* filePath): Plot(mInputList, filePath){}

void PlotAnaV0Mult::Make(){

	handleHistograms();


	outFile->Close();
	histFile->Close();
}

void PlotAnaV0Mult::Init(){
	//define the output file which will store all the canvases
	outFile = unique_ptr<TFile>( TFile::Open(outputPosition, "recreate") );

	if(!outFile || outFile->IsZombie()){
		cerr << "Could not open outFile. Leaving..." << endl;
		return;
	}

    histFile = unique_ptr<TFile>( TFile::Open("histFile.root", "read") );

	if(!histFile || histFile->IsZombie()){
		cerr << "Could not open histFile. Leaving..." << endl;
		return;
	}

	//load the tree chain from the input file
	ConnectInputTree(inputPosition, nameOfAnaV0MultTree);

    if(!tree){
    	cerr << "Couldn't open tree with data. Returning." << endl;
    	return;
    }
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

    DrawSTARpp510();
    CreateLegend(&legend,0.15, 0.7, 0.5, 0.6);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    //legend->AddEntry(func, "Fit signal + bcg","ple");
    legend->Draw("same");

    outFile->cd();
    canvas->Write();
    signalFinal->Write();
}
