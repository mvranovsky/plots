#include "../include/PlotV0SingleState.h"


PlotV0SingleState::PlotV0SingleState(TFile *mOutFile, const string mInputList, const char* filePath): Plot(mOutFile, mInputList, filePath){}

void PlotV0SingleState::Make(){

	// create canvases for all TH1D 
	hists1D = GetAllTH1D();
	if(hists1D.size() == 0){
		cerr << "Couldn't load 1D histograms from file. Leaving..." << endl;
		return;
	}
	for (int i = 0; i < hists1D.size(); ++i){
		if(!hists1D[i].first){
			cerr << "Couldn't load histogram " << hists1D[i].second << ". Leaving..." << endl;
			return;
		}
		TH1DGeneral(hists1D[i].second, hists1D[i].first);
		cout << "Created canvas for 1D histogram " << hists1D[i].second << endl;
	}



	// create canvases for all TH2F
	hists2D = GetAllTH2F();
	if(hists2D.size() == 0){
		cerr << "Couldn't load 2D histograms from file. Leaving..." << endl;
		return;
	}
	for (int i = 0; i < hists2D.size(); ++i){
		if(!hists2D[i].first){
			cerr << "Couldn't load histogram " << hists2D[i].second << ". Leaving..." << endl;
			return;
		}
		TH2FGeneral(hists2D[i].second, hists2D[i].first);
		cout << "Created canvas for 2D histogram " << hists2D[i].second << endl;
	}


	// create invariant mass plots for K0s
	invMassK0s(0.3, 1., 70, true);
	fitK0s(0.46, 0.53, 100, true);
	invMassK0s(0.3, 1., 70, false);
	fitK0s(0.46, 0.53, 100, false);
	cout << "Successfully created invariant mass plot and fit of peak for K0S..." << endl;


	// create invariant mass plots for Lambda
	invMassLambda(0.9,1.7,80, true);
	fitLambda(1.1, 1.133, 60, true);
	invMassLambda(0.9,1.7,80, false);
	fitLambda(1.1, 1.133, 60, false);
	cout << "Successfully created invariant mass plot and fit of peak for Lambda..." << endl;
	

	Clear();
	cout << "All histograms successfully saved to canvases..." << endl;
	cout << "The output file is saved: " << outputPosition << endl;

}

void PlotV0SingleState::Init(){
	//define the output file which will store all the canvases
	outFile = new TFile(outputPosition, "recreate");

	if(!outFile){
		cerr << "Couldn't open output file with position: " << outputPosition << endl;
	}
	// define which analysis is supposed to run
    if(!defineAnalysis())
    	cout << "No analysis defined." << endl;
	//load the tree chain from the input file
	ConnectInputTree(inputPosition);

    if(!tree){
    	cerr << "Couldn't open tree with data. Returning." << endl;
    	return;
    }
}


void PlotV0SingleState::invMassLambda(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF){

   	//TCanvas* canvas;
	CreateCanvas(&canvas,"invMassLambda", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassLambdaSignal", "invMassLambdaSignal", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{p #pi} [GeV/c^{2}]");
    TH1D *bcgFinal = new TH1D("invMassLambdaBcg", "invMassLambdaBcg", nBins, minRange, maxRange);
    bcgFinal->GetXaxis()->SetTitle("m_{p #pi} [GeV/c^{2}]");

    
    TString cmd, condition;
    //draw the pairs invMasses which satisfy condition
    cmd = TString::Format("invMass>>hist(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("(pairID == 1 || pairID ==2) && totQ == 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("(pairID == 1 || pairID ==2) && totQ == 0 && tofHit0 + tofHit1 == 1");
    }
    cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
    tree->Draw(cmd, condition);
    signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

	cmd = TString::Format("invMass>>bcg(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("(pairID == 1 || pairID ==2) && totQ != 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("(pairID == 1 || pairID ==2) && totQ != 0 && tofHit0 + tofHit1 == 1");
    }
    cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
    tree->Draw(cmd, condition);
    bcgFinal->Add((TH1D*)gPad->GetPrimitive(TString("bcg")));
    

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
    if(is2TOF){
    	canvas->Write("invMassLambda" + TString(2));
    	signalFinal->Write(signalFinal->GetName() + TString(2));
    	bcgFinal->Write(bcgFinal->GetName() + TString(2));
    }else{
    	canvas->Write("invMassLambda" + TString(1));
    	signalFinal->Write(signalFinal->GetName() + TString(1));
    	bcgFinal->Write(bcgFinal->GetName() + TString(1));
    }

}





void PlotV0SingleState::invMassK0s(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF){
	//TCanvas* canvas;
	CreateCanvas(&canvas,"invMassK0s", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassK0sSignal", "invMassK0sSignal", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{#pi #pi} [GeV/c^{2}]");
    TH1D *bcgFinal = new TH1D("invMassK0sBcg", "invMassK0sBcg", nBins, minRange, maxRange);
    bcgFinal->GetXaxis()->SetTitle("m_{#pi #pi} [GeV/c^{2}]");

    TString cmd, condition;
    // one state per event
    
    //draw the pairs invMasses which satisfy condition
    cmd = TString::Format("invMass>>hist(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 1");
    }
    cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
    tree->Draw(cmd, condition);
    signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
    
    cmd = TString::Format("invMass>>bcg(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("pairID == 0 && totQ != 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("pairID == 0 && totQ != 0 && tofHit0 + tofHit1 == 1");
    }
    cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
	tree->Draw(cmd, condition);
    bcgFinal->Add((TH1D*)gPad->GetPrimitive(TString("bcg")));
    
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
    if(is2TOF){
    	canvas->Write("invMassK0s" + TString(2));
    	signalFinal->Write(signalFinal->GetName() + TString(2));
    	bcgFinal->Write(bcgFinal->GetName() + TString(2));
    }else{
    	canvas->Write("invMassK0s" + TString(1));
    	signalFinal->Write(signalFinal->GetName() + TString(1));
    	bcgFinal->Write(bcgFinal->GetName() + TString(1));
    }

	
}

void PlotV0SingleState::fitK0s(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF) {
	CreateCanvas(&canvas,"fitK0s", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("fitPeakK0sSignalFit", "fitPeakK0sSignalFit", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    
    TString cmd, condition;
    cout << "fitting peak for K0S" << endl;

    //draw the pairs invMasses which satisfy condition
    cmd = TString::Format("invMass>>hist(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 1");
    }
	cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
    tree->Draw(cmd, condition);
	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
	
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
    if(is2TOF){
    	canvas->Write("invMassK0sFit" + TString(2));
    	signalFinal->Write(signalFinal->GetName() + TString(2));
    }else{
    	canvas->Write("invMassK0sFit" + TString(1));
    	signalFinal->Write(signalFinal->GetName() + TString(1));
    }
}

void PlotV0SingleState::fitLambda(Double_t minRange, Double_t maxRange, int nBins, bool is2TOF) {
	CreateCanvas(&canvas,"fitLambda", widthTypical, heightTypical);
	SetGPad(); //might set log y scale, watchout

    TH1D *signalFinal = new TH1D("invMassLambdaSignalFit", "invMassLambdaSignalFit", nBins, minRange, maxRange);
    signalFinal->GetXaxis()->SetTitle("m_{p #pi^{+}} [GeV/c^{2}]");
    
    TString cmd, condition;
    //draw the pairs invMasses which satisfy condition
    cmd = TString::Format("invMass>>hist(%d, %f, %f)", nBins, minRange,maxRange);
    if(is2TOF){
		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 2");
    } else{
  		condition = TString("pairID == 0 && totQ == 0 && tofHit0 + tofHit1 == 1");
    }
    cout << "Command: " << cmd << endl;
    cout << "Condition: " << condition << endl;
    tree->Draw(cmd, condition);
	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
	
	canvas->Clear();

	SetHistStyle(signalFinal, kBlack, markerStyleTypical);
	signalFinal->GetYaxis()->SetTitleOffset(1.5);
    signalFinal->Draw();
	int y = fitGaussPol2(&signalFinal,(maxRange-minRange)/nBins, minRange, maxRange, 0.01, 0.01, 0.01, 5000, 1.15, 0.01);
    CreateText("Yield: " + TString(y), 0.15, 0.8,0.4,0.8);

    DrawSTARpp510();
    CreateLegend(&legend,0.15, 0.7, 0.5, 0.6);
    legend->AddEntry(signalFinal, "Data (unlike-sign pairs)","ple");
    //legend->AddEntry(func, "Fit signal + bcg","ple");
    legend->Draw("same");

    outFile->cd();
    if(is2TOF){
    	canvas->Write("invMassLambdaFit" + TString(2));
    	signalFinal->Write(signalFinal->GetName() + TString(2));
    }else{
    	canvas->Write("invMassLambdaFit" + TString(1));
    	signalFinal->Write(signalFinal->GetName() + TString(1));
    }
}
