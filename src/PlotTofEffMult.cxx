#include "../include/PlotTofEffMult.h"


PlotTofEffMult::PlotTofEffMult(const string mInputList, const char* filePath): Plot( mInputList, filePath){}

void PlotTofEffMult::Make(){

	
	TH1D* h1 = (TH1D*)inFile->Get("hInvMassTof1");
	TH1D* h2 = (TH1D*)inFile->Get("hInvMassTof2");

	if(h1 && h2){
		outFile->cd();
		h1->Write();
		h2->Write();
	}else{
		cout << "Did not get h1, h2" << endl;
	}
	cout << "input position: " << inputPosition << endl;	

	/*
		// create invariant mass plots for K0s
	if( inputPosition.find("MC") != string::npos || inputPosition.find("Control") != string::npos){
		fitK0s(0.46, 0.53, 35, true);
		fitK0s(0.46, 0.53, 35, false);
	}else{
		fitK0s(0.45, 0.54, 45, true);
		fitK0s(0.45, 0.54, 45, false);
	}

	cout << "Successfully created invariant mass plot and fit of peak for K0S..." << endl;
	*/
    
	efficiency(1,0); //eta
	efficiency(2,0); //phi
	efficiency(3,0); //pT
	efficiency(4,0); //Vz
	/*

	efficiency(5,0); //Dca daughters
	efficiency(6,0); //Dca Beamline
	efficiency(7,0); //pointing angle
	efficiency(8,0); //decay length

	efficiency(3, 1); //pT for |eta|<0.2
	efficiency(3, 2); //pT for |eta|>0.2

	// efficiencies
	efficiency(1,1); //eta
	efficiency(2,1); //phi
	efficiency(3,1); //pT
	efficiency(4,1); //Vz

	efficiency(1,2); //eta
	efficiency(2,2); //phi
	efficiency(3,2); //pT
	efficiency(4,2); //Vz
	*/
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

	for (int i = 0; i < histsParFit.size(); ++i){
		if(!histsParFit[i]){
			cerr << "Couldn't load histogram from fit parameters. Leaving..." << endl;
			return;
		}
		TH1DGeneral("", histsParFit[i]);
		cout << "Created canvas for 1D histogram " << histsParFit[i]->GetName() << endl;
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



	outFile->Close();
	Clear();
	MCAna->Close();
	cout << "All histograms successfully saved to canvases..." << endl;
	cout << "The output file is saved: " << outputPosition << endl;


}

void PlotTofEffMult::Init(){
	//define the output file which will store all the canvases
	outFile = unique_ptr<TFile>(TFile::Open(outputPosition, "recreate"));

	if(!outFile || outFile->IsZombie()){
		cerr << "Couldn't open output file with position: " << outputPosition << endl;
	}

	//load the tree chain from the input file
	if(!ConnectInputTree(inputPosition, nameOfTofEffMultTree) || !tree ){
		cout << "Could not connect input tree correctly." << endl;
		return;
	}


    if( inputPosition.find("MC") != string::npos || inputPosition.find("Control") != string::npos ){
    	cout << "This is an MC simulation or running Control!" << endl;

	    MCAna = new TFile("MCAnalysis.root", "recreate");

    	nBins = 35;
    	lowRange = 0.46;
    	topRange = 0.53;

    	innerBorder[0] = 10;
    	innerBorder[1] = 26;
    	outerBorder[0] = 2;
    	outerBorder[1] = 34;
    }else{
	    MCAna = new TFile("MCAnalysis.root", "read");    	
    }

    if(!MCAna){
    	cout << "Couldn't open MC ana file. Leaving..." << endl;
    	return;
    }

    outFile->cd();
    outFile->mkdir("PLUS");
    outFile->mkdir("MINUS");

}



void PlotTofEffMult::fitK0s(Double_t minRange, Double_t maxRange, int numBins, bool is2TOF) {

    TH1D *signalFinal = new TH1D("fitPeakK0sSignalFit", "fitPeakK0sSignalFit", numBins, minRange, maxRange);

    signalFinal->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    
    TString cmd, condition;
    cout << "fitting peak for K0S" << endl;

    //draw the pairs invMasses which satisfy condition

    for (int i = 0; i < 5; ++i){

	   	tree->Draw(TString::Format("invMass%d>>hist(%d, %f, %f)", i,numBins, minRange,maxRange), TString::Format("tofHit%d > 0 && tofHit%d > 0", 2*i, 2*i +1) );
	   	//tree->Draw(TString::Format("invMass%d>>+hist", i), TString::Format("tofHit%d > 0 && tofHit%d > 0", 2*i+1, 2*i));

		if(!is2TOF){
	   		tree->Draw(TString::Format("invMass%d>>+hist",i),TString::Format("tofHit%d > 0 && tofHit%d < 0",2*i,  2*i +1) );
		}
	   	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
    }


    outFile->cd();
    if(is2TOF){
    	signalFinal->Write("invMassK0sFit2");
    }else{
    	signalFinal->Write("invMassK0sFit1");
    }
}

void PlotTofEffMult::fitLambda(Double_t minRange, Double_t maxRange, int numBins, bool is2TOF) {

    TH1D *signalFinal = new TH1D("invMassLambdaSignalFit", "invMassLambdaSignalFit", numBins, minRange, maxRange);
    TH1D *histCheck = new TH1D("histCheck", "", numBins, minRange, maxRange);

    signalFinal->GetXaxis()->SetTitle("m_{p #pi^{+}} [GeV/c^{2}]");
    
   	tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "( pairID == 1 || pairID == 2 ) && totQ == 0 && tofHit0 > 0 && tofHit1 > 0");
	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
   	tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "( pairID == 1 || pairID == 2 ) && totQ == 0 && tofHit1 > 0 && tofHit0 > 0");
   	signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
	if(!is2TOF){
   		tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "( pairID == 1 || pairID == 2 ) && totQ == 0 && tofHit0 > 0 && tofHit1 < 0");
   		signalFinal->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );

   		tree->Draw(TString::Format("invMass>>hist(%d, %f, %f)", numBins, minRange,maxRange), "( pairID == 1 || pairID == 2 ) && totQ == 0 && tofHit1 > 0 && tofHit0 < 0");
   		histCheck->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
   		if(histCheck->GetEntries() != 0){
   			cerr << "somethings wrong with histCheck. Leaving..." << endl;
   			return;
   		}
	}


    outFile->cd();
    if(is2TOF){
    	signalFinal->Write("invMassLambdaFit2");
    }else{
    	signalFinal->Write("invMassLambdaFit1");
    }
}

double PlotTofEffMult::GoodnessOfFit(RooPlot*& frame, RooAddPdf& model, RooDataHist& data){
    double chiSquare = frame->chiSquare();  // This calculates chi-square per degree of freedom

    int numBins = frame->GetNbinsX();  // Number of bins used in the histogram
    int nParams = model.getParameters(data)->getSize();  // Number of floating parameters in the model
    int ndf = numBins - nParams;

    return chiSquare/ndf;
}


TString PlotTofEffMult::convertToString(double val) {

    ostringstream streamA;
    streamA << fixed << setprecision(3) << val;
    TString formattedA = streamA.str();

    return formattedA;
}

int PlotTofEffMult::makeInt(double val) {
	int result = val;
	return result;
}

vector<TString> PlotTofEffMult::getConditions(int Switch, int i, double Min, double Max){
	vector<TString> cond;
	if(Switch == 1){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && etaHadron%d > %f && etaHadron%d < %f", 2*i, 2*i + 1, 2*i + 1, Min, 2*i + 1, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && etaHadron%d > %f && etaHadron%d < %f", 2*i, 2*i + 1, 2*i, Min, 2*i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && etaHadron%d > %f && etaHadron%d < %f", 2*i, 2*i + 1, 2*i + 1, Min, 2*i + 1, Max) );
	}else if(Switch == 2){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && phiHadron%d > %f && phiHadron%d < %f", 2*i, 2*i + 1, 2*i + 1, Min*TMath::Pi()/180, 2*i + 1, Max*TMath::Pi()/180) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && phiHadron%d > %f && phiHadron%d < %f", 2*i, 2*i + 1, 2*i, Min*TMath::Pi()/180, 2*i, Max*TMath::Pi()/180) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && phiHadron%d > %f && phiHadron%d < %f", 2*i, 2*i + 1, 2*i + 1, Min*TMath::Pi()/180, 2*i + 1, Max*TMath::Pi()/180) );
	}else if(Switch == 3){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && pTInGev%d > %f && pTInGev%d < %f", 2*i, 2*i + 1, 2*i + 1, Min, 2*i + 1, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && pTInGev%d > %f && pTInGev%d < %f", 2*i, 2*i + 1, 2*i, Min, 2*i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && pTInGev%d > %f && pTInGev%d < %f", 2*i, 2*i + 1, 2*i + 1, Min, 2*i + 1, Max) );
	}else if(Switch == 4){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && vertexZInCm%d > %f && vertexZInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && vertexZInCm%d > %f && vertexZInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && vertexZInCm%d > %f && vertexZInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
	}else if(Switch == 5){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && dcaDaughtersInCm%d > %f && dcaDaughtersInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && dcaDaughtersInCm%d > %f && dcaDaughtersInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && dcaDaughtersInCm%d > %f && dcaDaughtersInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
	}else if(Switch == 6){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && dcaBeamlineInCm%d > %f && dcaBeamlineInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && dcaBeamlineInCm%d > %f && dcaBeamlineInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && dcaBeamlineInCm%d > %f && dcaBeamlineInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
	}else if(Switch == 7){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && PointingAngle%d > %f && PointingAngle%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && PointingAngle%d > %f && PointingAngle%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && PointingAngle%d > %f && PointingAngle%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
	}else if(Switch == 8){
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && decayLengthInCm%d > %f && decayLengthInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d > 0 && decayLengthInCm%d > %f && decayLengthInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
		cond.push_back( TString::Format("tofHit%d > 0 && tofHit%d < 0 && decayLengthInCm%d > %f && decayLengthInCm%d < %f", 2*i, 2*i + 1,i, Min,i, Max) );
	}

	return cond;

}


TH1D* PlotTofEffMult::getHist(int Switch, double Min, double Max, bool is2TOF, int runSeparatePions, bool runFirstOnly){

	// function which integrates the signal, it is a universal function for all variables
    TString variable, variableString, condition1, condition2;
    if( Switch == 1){//eta
    	variable = "eta";
    	variableString = "#eta";
    } else if( Switch == 2){
    	variable = "phi";
    	variableString = "#phi";
    } else if(Switch == 3){
    	variable = "pT";
    	variableString = "p_{T}";
	}else if(Switch == 4){
    	variable = "Vz";
    	variableString = "V_{Z}";
    }else if(Switch == 5){
    	variable = "DCAD";
    	variableString = "DCA_{Daughters}";
    }else if(Switch == 6){
    	variable = "DCAB";
    	variableString = "DCA_{Beamline}";
    }else if(Switch == 7){
    	variable = "PointAng";
    	variableString = "#phi_{pointingAngle}";
    }else if(Switch == 8){
    	variable = "DecayL";
    	variableString = "L_{decay}";
    }

	TH1D* hist;
    // create histograms
    if(is2TOF){
   		hist = new TH1D("effK0s2" + variable + convertToString( Min ) + TString("to") + convertToString( Max ), "fit peak for a specific bin;m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; counts", nBins, lowRange, topRange);
    }else{
   		hist = new TH1D("effK0s1" + variable + convertToString( Min ) + TString("to") + convertToString( Max ), "fit peak for a specific bin;m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; counts", nBins, lowRange, topRange);
   	}
	

   	for (int i = 0; i < 5; ++i){
   		vector<TString> conditions = getConditions(Switch, i, Min, Max);

   		if(runSeparatePions == 1 && Switch == 3){
   			conditions[0] += TString(" && abs(etaHadron%d) <= 0.2", 2*i +1);
   			conditions[1] += TString(" && abs(etaHadron%d) <= 0.2", 2*i);
   			conditions[2] += TString(" && abs(etaHadron%d) <= 0.2", 2*i + 1);
   		}else if(runSeparatePions == 2 && Switch == 3){
   			conditions[0] += TString(" && abs(etaHadron%d) > 0.2", 2*i +1);
   			conditions[1] += TString(" && abs(etaHadron%d) > 0.2", 2*i);
   			conditions[2] += TString(" && abs(etaHadron%d) > 0.2", 2*i + 1);   			
   		}

   		cout << conditions[0] << endl;
   		//cout << conditions[1] << endl;
   		cout << conditions[2] << endl;

   		tree->Draw( TString::Format("invMass%d>>hist(%d, %f, %f)",i, nBins, lowRange, topRange) , conditions[0] );
   		//tree->Draw( TString::Format("invMass%d>>+hist", i), conditions[1] );
   		if(!is2TOF){
   			tree->Draw( TString::Format("invMass%d>>+hist", i), conditions[2] );
   		}

		hist->Add((TH1D*)gPad->GetPrimitive( TString("hist") ) );
		if(runFirstOnly)
			break;
   	}


    //cout << "Entries in hist: " << hist->GetEntries() <<  endl;

    if(hist->GetEntries() == 0){
    	return nullptr;
    }else{
    	return hist;
    }
}

vector<pair<int,double>> PlotTofEffMult::effFit(int Switch ,double Min, double Max,Double_t signalGuess1, Double_t signalGuess2, Double_t polGuess1[2], Double_t polGuess2[2], int runSeparatePions){
    

	TH1D *hist1 = getHist(Switch, Min, Max, false, runSeparatePions, false);
	TH1D *hist2 = getHist(Switch, Min, Max, true, runSeparatePions, false);



    TString variable, variableString, condition1, condition2;
    if( Switch == 1){//eta
    	variable = "eta";
    	variableString = "#eta";
    } else if( Switch == 2){
    	variable = "phi";
    	variableString = "#phi";
    } else if(Switch == 3){
    	variable = "pT";
    	variableString = "p_{T}";
    }else if(Switch == 4){
    	variable = "Vz";
    	variableString = "V_{Z}";
    }else if(Switch == 5){
    	variable = "DCAD";
    	variableString = "DCA_{Daughters}";
    }else if(Switch == 6){
    	variable = "DCAB";
    	variableString = "DCA_{Beamline}";
    }else if(Switch == 7){
    	variable = "PointAng";
    	variableString = "#phi_{pointingAngle}";
    }else if(Switch == 8){
    	variable = "DecayL";
    	variableString = "L_{decay}";
    }

	//TH1D* hist1 = new TH1D("effK0s1" + variable + convertToString( Min ) + TString("to") + convertToString( Max ), "fit peak for a specific bin;m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; counts", nBins, lowRange, topRange);
	//TH1D* hist2 = new TH1D("effK0s2" + variable + convertToString( Min ) + TString("to") + convertToString( Max ), "fit peak for a specific bin;m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; counts", nBins, lowRange, topRange);

	double normalisation = (topRange - lowRange)/nBins;
	// considering bins from 0.45 to 0.54, signal from 0.48 to 0.52, background 0.46-0.48 and 0.52-0.54
	double sumSignal1 = 0;
	double sumBcg1 = 0;
	double sumSignal2 = 0;
	double sumBcg2 = 0;

	double sumSignal1Err = 0;
	double sumBcg1Err = 0;
	double sumSignal2Err = 0;
	double sumBcg2Err = 0;

	for (int iBin = outerBorder[0]; iBin < outerBorder[1]; ++iBin){
		if (iBin < innerBorder[0] || iBin >= innerBorder[1] ){
			sumBcg1 += hist1->GetBinContent(iBin);
			sumBcg2 += hist2->GetBinContent(iBin);
			sumBcg1Err += pow(hist1->GetBinError(iBin),2);
			sumBcg2Err += pow(hist2->GetBinError(iBin),2);
		} else if( iBin >= innerBorder[0] && iBin < innerBorder[1] ){
			sumSignal1 += hist1->GetBinContent(iBin);
			sumSignal2 += hist2->GetBinContent(iBin);
			sumSignal1Err += pow(hist1->GetBinError(iBin),2);
			sumSignal2Err += pow(hist2->GetBinError(iBin),2);

		}
	}
    // Define the observable
	/*
    RooRealVar x1("x1", "Observable1", lowRange,topRange);
    RooRealVar x2("x2", "Observable2", lowRange,topRange);

    // Convert TH1D to RooDataHist
    RooDataHist data1("data1", "dataset with x1", x1, Import(*hist1));
    RooDataHist data2("data2", "dataset with x2", x2, Import(*hist2));

    // Define the Crystal ball PDFs for the peaks
    RooRealVar mean1("mean1", "mean of Gaussian 1", 0.497611 , 0.495, 0.5); // adjust these ranges and initial values as needed
    RooRealVar sigma1("sigma1", "width of Gaussian 1", 0.004 , 0.003, 0.005);
    RooRealVar alpha1("alpha1", "Tail parameter 1", 3, 1.,8);
    RooRealVar n1("n1", "Power parameter 1", 4, 1., 8);
    RooCBShape CBfunc1("cb1", "Crystal Ball PDF 1", x1, mean1, sigma1, alpha1, n1);

    RooRealVar mean2("mean2", "mean of Gaussian 2",  0.497611 , 0.495, 0.5); // adjust these ranges and initial values as needed
    RooRealVar sigma2("sigma2", "width of Gaussian 2", 0.004 , 0.003, 0.005);
    RooRealVar alpha2("alpha2", "Tail parameter 2", 3, 1., 8);
    RooRealVar n2("n2", "Power parameter 2", 4, 1., 8);
    RooCBShape CBfunc2("cb2", "Crystal Ball PDF 2", x2, mean2, sigma2, alpha2, n2);


    // Define the polynomial PDF for the background
    RooRealVar coef1_0("coef1_0", "constant term 1", polGuess1[0], -50000, 50000);
    RooRealVar coef1_1("coef1_1", "linear term 1", polGuess1[1], -50000, 50000);
    RooPolynomial poly1("poly1", "Background Polynomial 1", x1, RooArgList(coef1_0, coef1_1)); //prehodit

    RooRealVar coef2_0("coef2_0", "constant term 2", polGuess2[0], -50000, 50000);
    RooRealVar coef2_1("coef2_1", "linear term 2", polGuess2[1], -50000, 50000);
    RooPolynomial poly2("poly2", "Background Polynomial 2", x2, RooArgList(coef2_0, coef2_1)); //prehodit poradie coef


    // Coefficients for PDFs
    RooRealVar nsig1("nsig1", "signal 1", signalGuess1, 0, 5e+5);
    RooRealVar nbkg1("nbkg1", "Background 1",300000, 100, 6e+7);
    RooAddPdf model1("model1", "Crystal ball plus polynomial 1", RooArgList(CBfunc1, poly1), RooArgList(nsig1, nbkg1)); //,poly   gaussFrac2

    RooRealVar nsig2("nsig2", "signal 2", signalGuess2, 0, 2e+5);
    RooRealVar nbkg2("nbkg2", "Background 2",1000000, 100, 1e+7);
    RooAddPdf model2("model2", "Crystal ball plus polynomial 2", RooArgList(CBfunc2, poly2), RooArgList(nsig2, nbkg2)); //,poly   gaussFrac2


    // Fit the model to the data
    model1.fitTo(data1, PrintLevel(-1));
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;

    model2.fitTo(data2, PrintLevel(-1));


    // create a canvas that will hold both fits
    TCanvas* c = new TCanvas(TString("eff") + variableString + TString("from") + convertToString(Min) + TString("to") + convertToString(Max) , "Fit Result", 1200, 800); 
    c->Divide(2,1);

    c->cd(1);
    SetGPad();
    RooPlot* frame1 = x1.frame();
    frame1->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitleOffset(1.8);
    data1.plotOn(frame1);
    model1.plotOn(frame1);
    model1.plotOn(frame1, Components(poly1), LineStyle(kDashed), LineColor(kRed));
    model1.plotOn(frame1, Components(CBfunc1), LineStyle(kDashed), LineColor(kBlue));
    frame1->Draw();



    // create a legend with results of fit
    TPaveText* text_1 = new TPaveText(0.15,0.6, 0.4, 0.4,"brNDC");
    text_1 -> SetTextSize(textSize);
    text_1 -> SetTextAlign(11);
    text_1 -> SetFillColor(0);
    text_1 -> SetTextFont(62);
    TString meanFinal1 = "#mu_{1} = " + to_string( mean1.getVal() ) + " +/- " + to_string(mean1.getError());
    TString sigmaFinal1 = "#sigma_{1} = " + to_string(sigma1.getVal()) + " +/- " + to_string(sigma1.getError());
    TString fitFinal1 = "#chi^{2}/NDF = " + to_string( GoodnessOfFit(frame1, model1, data1) );
    TString yieldText1 = "Yield = " + to_string( ceil(nsig1.getVal())) + "+/-" + to_string(ceil(nsig1.getError()));
    TString pol1_0 = "Pol_{0} = " + to_string( makeInt(coef1_0.getVal()) );
    TString pol1_1 = "Pol_{1} = " + to_string( makeInt(coef1_1 .getVal()) ); 
    TString nText1 = "n_{1} = " + to_string(n1.getVal()) + "+/-" + to_string( n1.getError() );
    TString alphaText1 = "#alpha_{1} = " + to_string(alpha1.getVal() ) + "+/-" + to_string(alpha1.getError());
    //text_1 -> AddText( pol1_0 );
    //text_1 -> AddText( pol1_1 );
    text_1 -> AddText(meanFinal1);
    text_1 -> AddText(sigmaFinal1);
    text_1 -> AddText( nText1 );
    text_1 -> AddText( alphaText1);
    text_1 -> AddText(fitFinal1);
    text_1 -> AddText(yieldText1);
    text_1 -> Draw("same");
    
    DrawSTARpp510(0.23, 0.9, 0.43, 0.9);
    CreateText("at least 1 ToF track",0.23, 0.87, 0.43, 0.87 );    
    CreateText(convertToString(Min) + TString(" < ") + variableString + TString(" < ") + convertToString(Max),0.23, 0.84, 0.43, 0.84 ); 
    TLegend *leg1;
    CreateLegend(&leg1, 0.7,0.6,0.9,0.5);
    leg1 -> AddEntry("data1", "Data", "lep");
    leg1 -> AddEntry("model1", "Total fit", "l");
    leg1 -> AddEntry("cb1", "Crystal ball", "l");
    leg1 -> Draw("same");

//--------------------------------------------------------------------------------------------------------

    c->cd(2);
    SetGPad();
    RooPlot* frame2 = x2.frame();
    frame2->GetXaxis()->SetTitle("m_{#pi^{+} #pi^{-}} [GeV/c^{2}]");
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitleOffset(1.8);
    data2.plotOn(frame2);
    model2.plotOn(frame2);
    model2.plotOn(frame2, Components(poly2), LineStyle(kDashed), LineColor(kRed));
    model2.plotOn(frame2, Components(CBfunc2), LineStyle(kDashed), LineColor(kBlue));
    //frame->GetYaxis()->SetRangeUser(yRangeBottom, yRangeTop);
    frame2->Draw();



    // create a legend with results of fit
    TPaveText* text_2 = new TPaveText(0.17,0.9, 0.3,0.75,"brNDC");
    text_2 -> SetTextSize(textSize );
    text_2 -> SetTextAlign(11);
    text_2 -> SetFillColor(0);
    text_2 -> SetTextFont(62);
    TString meanFinal2 = "#mu_{2} = " + to_string(mean2.getVal()) + " +/- " + to_string(mean2.getError());
    TString sigmaFinal2 = "#sigma_{2} = " + to_string(sigma2.getVal()) + " +/- " + to_string(sigma2.getError());
    TString fitFinal2 = "#chi^{2}/NDF = " + to_string( GoodnessOfFit(frame2, model2, data2) );
    TString yieldText2 = "Yield = " + to_string( makeInt(ceil(nsig2.getVal()))) + "+/-" + to_string(ceil(nsig2.getError()));
    TString pol2_0 = "Pol_{0} = " + to_string( makeInt(coef2_0.getVal()) );
    TString pol2_1 = "Pol_{1} = " + to_string( makeInt(coef2_1 .getVal()) ); 
    TString nText2 = "n_{2} = " + to_string(n2.getVal()) + "+/-" + to_string( n2.getError() );
    TString alphaText2 = "#alpha_{2} = " + to_string(alpha2.getVal() ) + "+/-" + to_string(alpha2.getError());
    //text_2 -> AddText( pol2_0 );
    //text_2 -> AddText( pol2_1 );
    text_2 -> AddText(meanFinal2);
    text_2 -> AddText(sigmaFinal2);
    text_2 -> AddText( nText2 );
    text_2 -> AddText( alphaText2);
    text_2 -> AddText(fitFinal2);
    text_2 -> AddText(yieldText2);
    text_2 -> Draw("same");

    DrawSTARpp510(0.65, 0.9, 0.85, 0.9);
    CreateText("2 ToF tracks",0.65, 0.87, 0.85, 0.87 );
    CreateText(convertToString(Min) + TString(" < ") + variable + TString(" < ") + convertToString(Max) ,0.65, 0.84, 0.85, 0.84 ); 
    TLegend *leg2;
    CreateLegend(&leg2, 0.7,0.72,0.9,0.62);
    leg2 -> AddEntry("data2", "Data", "lep");
    leg2 -> AddEntry("model2", "Total fit", "l");
    leg2 -> AddEntry("cb2", "Crystal ball", "l");
    leg2 -> Draw("same");

    if(runSeparatePions == 0){
    	outFile->cd();
    	c->Write(); 
    }else if(runSeparatePions == 1){
        outFile->cd("PLUS");
    	c->Write(); 
    	outFile->cd();
    }else if(runSeparatePions == 2){
    	outFile->cd("MINUS");
    	c->Write(); 
    	outFile->cd();
    }
    result.push_back( make_pair(nsig1.getVal(), nsig1.getError() ) );
    result.push_back( make_pair(nsig2.getVal(), nsig2.getError() ) );
    
    */
    vector<pair<int,double>> result;

    result.push_back( make_pair(1,0) );
    result.push_back( make_pair(0,0) );

	result.push_back( make_pair( sumSignal1 - sumBcg1, sqrt(sumSignal1Err + sumBcg1Err ) ) );
	result.push_back( make_pair( sumSignal2 - sumBcg2, sqrt(sumSignal2Err + sumBcg2Err ) ) );


    return result;

}


void PlotTofEffMult::efficiency(int switcher, int runSeparatePions = 0) {
	// runSeparatePions = 0: no
	// runSeparatePions = 1: yes, run pi +
	// runSeparatePions = 2: yes, run pi -

	// switcher
	int binCount;
	double top, bottom, difference;
    TString variable, variableString, title;
    if( switcher == 1){//eta
    	variable = "eta";
    	variableString = "#eta";
    	title = "#eta [-]";
    	binCount = 9;
    	top = 0.9;
    	bottom = -0.9;
    	difference = 0.2;
    } else if( switcher == 2){
    	variable = "phi";
    	variableString = "#phi";
       	title = "#phi [-]";
    	binCount = 6;
    	top = 180;
    	bottom = -180;
    	difference = 60;
    } else if(switcher == 3){
    	variable = "pT";
    	variableString = "p_{T}";
    	title = "p_{T} [GeV/c]";
    	binCount = 6;
    	top = 1.4;
    	bottom = 0.2;
    	difference = 0.2;
    }else if(switcher == 4){
    	variable = "Vz";
    	variableString = "V_{Z}";
      	title = "V_{Z} [cm]";
    	binCount = 10;
    	top = 100;
    	bottom = -100;
    	difference = 20;
    }else if(switcher == 5){
    	variable = "DCADaughters";
    	variableString = "DCA_{Daughters}";
      	title = "DCA_{Daughters} [cm]";
    	binCount = 5;
    	top = 1.5;
    	bottom = 0;
    	difference = 0.3;    	
    }else if(switcher == 6){
    	variable = "DCABeamline";
    	variableString = "DCA_{Beamline}";
      	title = "DCA_{Beamline} [cm]";
    	binCount = 7;
    	top = 1.4;
    	bottom = 0;
    	difference = 0.2; 
    }else if(switcher == 7){
    	variable = "PointingAngle";
    	variableString = "cos(#phi_{pointingAngle})";
      	title = "cos( #varphi_{PointingAngle})";
    	binCount = 5;
    	top = 1;
    	bottom = 0.95;
    	difference = 0.01;
    }else if(switcher == 8){
    	variable = "DecayLength";
    	variableString = "L_{decay}";
      	title = "L_{decay} [cm]";
    	binCount = 10;
    	top = 10;
    	bottom = 0;
    	difference = 1;
    }


	cout << "Beginning efficiency plot for " << variable << endl;
	vector<pair<int,double>> binContent1, binContent2, current_yield, current_bcgSub, binContentBcgSub1, binContentBcgSub2;
	Double_t arg1[2], arg2[2], arg3[2], arg4[2];
	for (int iBin = 0; iBin < binCount; ++iBin){ // considering 18 bins fro -0.9 to 0.9, skipping from -0.2 to 0.2 
		//run fits over all bins
		cout << "Running " << variable << " bin number: " << iBin << endl;

		// switching between different guesses for fits
		if(switcher ==1){
			arg1[0] = effEtaPol0Guess1[iBin];
			arg1[1] = effEtaPol1Guess1[iBin];
			arg2[0] = effEtaPol0Guess2[iBin];
			arg2[1] = effEtaPol1Guess2[iBin];
			arg3[0] = effEtaAmpGuess1[iBin];
			arg3[1] = effEtaAmpGuess2[iBin];
		} else if(switcher == 2){
			arg1[0] = effPhiPol0Guess1[iBin];
			arg1[1] = effPhiPol1Guess1[iBin];
			arg2[0] = effPhiPol0Guess2[iBin];
			arg2[1] = effPhiPol1Guess2[iBin];
			arg3[0] = effPhiAmpGuess1[iBin];
			arg3[1] = effPhiAmpGuess2[iBin];
		} else if(switcher == 3){
			arg1[0] = effPtPol0Guess1[iBin];
			arg1[1] = effPtPol1Guess1[iBin];
			arg2[0] = effPtPol0Guess2[iBin];
			arg2[1] = effPtPol1Guess2[iBin];
			arg3[0] = effPtAmpGuess1[iBin];
			arg3[1] = effPtAmpGuess2[iBin];
		} else if(switcher == 4){
			arg1[0] = effVzPol0Guess1[iBin];
			arg1[1] = effVzPol1Guess1[iBin];
			arg2[0] = effVzPol0Guess2[iBin];
			arg2[1] = effVzPol1Guess2[iBin];
			arg3[0] = effVzAmpGuess1[iBin];
			arg3[1] = effVzAmpGuess2[iBin];
		}




		if(switcher == 3 && iBin == 4 ){
			current_yield = effFit(switcher,bottom + iBin*difference, bottom + difference*(2 + iBin),arg3[0], arg3[1], arg1, arg2, runSeparatePions); 
		}else if(switcher == 6 && iBin == 6){
			//data
			current_yield = effFit(switcher,bottom + iBin*difference, 1.5,arg3[0], arg3[1], arg1, arg2, runSeparatePions); 
		}
		else{
			current_yield = effFit(switcher,bottom + iBin*difference, bottom + difference*(1 + iBin),arg3[0], arg3[1], arg1, arg2, runSeparatePions); 
		}

		binContent1.push_back(current_yield[0]);
		binContent2.push_back(current_yield[1]);
		binContentBcgSub1.push_back(current_yield[2]);
		binContentBcgSub2.push_back(current_yield[3]);
	}

	TH1D* Yield1 = new TH1D("Yield1" + variableString, "Yield1",binCount, bottom,top );
	TH1D* Yield2 = new TH1D("Yield2" + variableString, "Yield2",binCount, bottom,top );

	TH1D *YieldBcgSub1, *YieldBcgSub2;
	if(switcher == 3){
	    double binEdges[] = {0.2, 0.4, 0.6, 0.8, 1.};
	    binCount = 4;
		YieldBcgSub2 = new TH1D("YieldBcgSub2" + variableString, "Yield2",binCount,binEdges );
		YieldBcgSub1 = new TH1D("YieldBcgSub1" + variableString, "Yield1", binCount,binEdges );
	}else if(switcher == 6){
		//data
	    double binEdges[] = {0.0, 0.2,0.4, 0.6, 0.8, 1.,1.2, 1.5};
		//MC
	    //double binEdges[] = {0.0, 0.2,0.4, 0.6, 0.8, 1.};
		YieldBcgSub2 = new TH1D("YieldBcgSub2" + variableString, "Yield2",binCount,binEdges );
		YieldBcgSub1 = new TH1D("YieldBcgSub1" + variableString, "Yield1", binCount,binEdges );

	}else{
		YieldBcgSub2 = new TH1D("YieldBcgSub2" + variableString, "Yield2",binCount, bottom,top );
		YieldBcgSub1 = new TH1D("YieldBcgSub1" + variableString, "Yield1",binCount, bottom,top );
	}


	// check for inconsistencies
	if(binContent1.size() != binContent2.size() || binContent1.size() != binContentBcgSub1.size() || binContentBcgSub1.size() != binContentBcgSub2.size() ){
		cerr << "Number of bins in " << variable << " efficiency is not correct. Leaving..." << endl;
		return;
	}

	//load values to histogram
	for (int i = 0; i < binCount; ++i){

		Yield1->SetBinContent(i+1,binContent1[i].first);
		Yield1->SetBinError(i+1, binContent1[i].second);
		Yield2->SetBinContent(i+1,binContent2[i].first);
		Yield2->SetBinError(i+1, binContent2[i].second);


		YieldBcgSub1->SetBinContent(i+1,binContentBcgSub1[i].first);
		YieldBcgSub1->SetBinError(i+1, binContentBcgSub1[i].second);
		YieldBcgSub2->SetBinContent(i+1,binContentBcgSub2[i].first);
		YieldBcgSub2->SetBinError(i+1, binContentBcgSub2[i].second);
	}

	cout << "-------------------------------------------------------------------------" << endl;
	cout << "RESULTS FOR INTEGRATION:" << endl;
	for (int i = 0; i < binContentBcgSub1.size(); ++i){
		cout << variable << " bin number " << i << ",contentBcgSub 1: value: " << binContentBcgSub1[i].first << ", error: " << binContentBcgSub1[i].second << endl;
		cout << variable << " bin number " << i << ",contentBcgSub 2: value: " << binContentBcgSub2[i].first << ", error: " << binContentBcgSub2[i].second << endl;
	}
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "RESULTS FOR FITS:" << endl;
	for (int i = 0; i < binContent1.size(); ++i){
		cout << variable << " bin number " << i << ",content 1: value: " << binContent1[i].first << ", error: " << binContent1[i].second << endl;
		cout << variable << " bin number " << i << ",content 2: value: " << binContent2[i].first << ", error: " << binContent2[i].second << endl;
	}
	cout << "-------------------------------------------------------------------------" << endl;

	outFile->cd();

	
	CreateCanvas(&canvas,"eff" + variable, widthTypical + 400, heightTypical);
	//SetGPad(); //might set log y scale, watchout
	SetGPad(false,0.11, 0.05, 0.14,0.05);
    
    TEfficiency* pEff = new TEfficiency(*Yield2,*Yield1);
    TGraphAsymmErrors *gr = pEff->CreateGraph();
    gr->GetYaxis()->SetTitle("efficiency");
    gr->GetXaxis()->SetTitle(title);	
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed); 
    gr->SetMarkerSize(2);
    gr->GetYaxis()->SetRangeUser(0,1);
    //gr->Draw("same EP");

	
    TEfficiency* pEffBcgSub = new TEfficiency(*YieldBcgSub2,*YieldBcgSub1);
    TGraphAsymmErrors *grBcgSub = pEffBcgSub->CreateGraph();
    grBcgSub->GetYaxis()->SetTitle("efficiency");
    grBcgSub->GetXaxis()->SetTitle(title);
    grBcgSub->GetYaxis()->SetTitleSize(0.06);
    grBcgSub->GetXaxis()->SetTitleSize(0.06);
    grBcgSub->GetYaxis()->SetTitleOffset(0.7);
    grBcgSub->SetMarkerStyle(21);
    grBcgSub->SetMarkerColor(kBlue);
    grBcgSub->SetLineColor(kBlue); 
    grBcgSub->SetMarkerSize(3);
    //grBcgSub->GetXaxis()->SetRangeUser(-1,1);
    grBcgSub->GetYaxis()->SetRangeUser(0,1);
    grBcgSub->Draw("AEP");



    DrawSTARpp510(0.23, 0.85, 0.43, 0.85, 0.02);
    CreateLegend(&legend,0.68, 0.87, 0.9, 0.73);
    //legend->AddEntry(gr, "Efficiency K^{0}_{S} (Fitting)","ple");
    legend->AddEntry(grBcgSub, "Efficiency (data)","ple");
    legend->SetTextSize(0.05);

	TGraphAsymmErrors *grMC;
    if(inputPosition.find("MC") == string::npos){

	    if(runSeparatePions == 0){
	    	grMC = (TGraphAsymmErrors*)MCAna->Get("effInt" + variable);
	    }else if(runSeparatePions == 1){
	    	grMC = (TGraphAsymmErrors*)MCAna->Get("PLUS/effInt" + variable);
	    }else if(runSeparatePions == 2){
	    	grMC = (TGraphAsymmErrors*)MCAna->Get("MINUS/effInt" + variable);
	    }
	    
	    if(grMC){
		    grMC->SetMarkerStyle(22);
		    grMC->SetMarkerColor(kRed);
		    grMC->SetLineColor(kRed); 
	        grMC->SetMarkerSize(3);
		    grMC->GetXaxis()->SetTitle(title);
		    grMC->GetYaxis()->SetTitle("efficiency");
		    grMC->Draw("same EP");
	    	legend->AddEntry(grMC, "Efficiency (MC)", "ple");

	    	// calculate systematic errors by comparing to reconstructed MC
    		//double sysError = CalculateSystematicError(switcher,grBcgSub, grMC);
    		//cout << "Systematic error from reconstructed MC: " << sysError << endl;

	    }else{
	    	cout << "Couldn't load MC eff. Running without it." << endl;
	    }
    }else{
    	MCAna->cd();
    	gr->Write("effFit" + variable);
    	grBcgSub->Write("effInt" + variable);
    }


    TGraphAsymmErrors *graph = getTrueMC(switcher);

    // calculate sys errors from trueMC
	double sysError = CalculateSystematicError(switcher,grBcgSub, graph);
	cout << "Systematic error from true MC: " << sysError << endl;

    legend->AddEntry(graph, "true MC", "lep");

    canvas->cd();
    legend->Draw("same");
    graph->Draw("same P");
    
    outFile->cd();
    graph->Write("trueMC");
    if(runSeparatePions == 0){
    	canvas->Write("eff" + variable);
		grBcgSub->Write("effInt" + variable);
		gr->Write("effFit" + variable);
    }else if(runSeparatePions == 1){
    	outFile->cd("PLUS");
    	canvas->Write("eff" + variable);
		grBcgSub->Write("effInt" + variable);
		gr->Write("effFit" + variable); 
		outFile->cd();   	
    }else{
    	outFile->cd("MINUS");
    	canvas->Write("eff" + variable);
		grBcgSub->Write("effInt" + variable);
		gr->Write("effFit" + variable); 
		outFile->cd();  
    }


}


// Function to calculate systematic error
double PlotTofEffMult::CalculateSystematicError(int switcher,TGraphAsymmErrors* gData, TGraphAsymmErrors* gMC) {
    if (!gData || !gMC) {
        std::cerr << "Error: Input graphs are null!" << std::endl;
        return -1;
    }

    // Ensure the graphs have the same number of points
    int nDataPoints = gData->GetN();
    int nMCPoints = gMC->GetN();
    
	TString variable, variableString;
    if( switcher == 1){//eta
    	variable = "eta";
    	variableString = "#eta [-]";
    } else if( switcher == 2){
    	variable = "phi";
    	variableString = "#phi [#circ ]";
    } else if(switcher == 3){
    	variable = "pTInGev";
    	variableString = "p_{T} [GeV/c]";
    }else if(switcher == 4){
    	variable = "vertexZInCm";
    	variableString = "V_{Z} [cm]";
    }

    if(switcher == 5){// for reconstructed MC: 3
    	nDataPoints = 4;
    }else if (nDataPoints != nMCPoints) {
        std::cerr << "Error: Data and MC graphs must have the same number of points!" << std::endl;
        return -1;
    }

    // Create a TGraphAsymmErrors for the ratio
    TGraphAsymmErrors* gRatio = new TGraphAsymmErrors(nDataPoints);

    for (int i = 0; i < nDataPoints; ++i) {
        double xData, yData, exlData, exhData, eylData, eyhData;
        double xMC, yMC, exlMC, exhMC, eylMC, eyhMC;

        // Get point and errors for Data
        gData->GetPoint(i, xData, yData);
        exlData = gData->GetErrorXlow(i);
        exhData = gData->GetErrorXhigh(i);
        eylData = gData->GetErrorYlow(i);
        eyhData = gData->GetErrorYhigh(i);

        // Get point and errors for MC
        gMC->GetPoint(i, xMC, yMC);
        exlMC = gMC->GetErrorXlow(i);
        exhMC = gMC->GetErrorXhigh(i);
        eylMC = gMC->GetErrorYlow(i);
        eyhMC = gMC->GetErrorYhigh(i);

        // Skip if MC efficiency is zero to avoid division by zero
        if (yMC == 0) continue;

        // Calculate ratio and errors
        double ratio = yData / yMC;
        double errLow = ratio * sqrt(pow(eylData / yData, 2) + pow(eylMC / yMC, 2));
        double errHigh = ratio * sqrt(pow(eyhData / yData, 2) + pow(eyhMC / yMC, 2));

        // Set the point and errors in the ratio graph
        gRatio->SetPoint(i, xData, ratio);
        gRatio->SetPointError(i, exlData, exhData, errLow, errHigh);
    }

    // Fit the ratio graph with a constant function
    TF1* constantFit = new TF1("constantFit", "[0]", gRatio->GetXaxis()->GetXmin(), gRatio->GetXaxis()->GetXmax());
    gRatio->Fit(constantFit, "Q"); // "Q" suppresses fit output

    // Extract the fit parameter and its uncertainty
    double systematicError = constantFit->GetParameter(0);
    double systematicErrorUncertainty = constantFit->GetParError(0);

    // Print results
    std::cout << "Systematic error (C): " << systematicError << std::endl;
    std::cout << "Systematic uncertainty: " << systematicErrorUncertainty << std::endl;

    // Optionally draw the ratio graph and fit
    TCanvas* c1 = new TCanvas("c1", "Systematic Error", 800, 600);
    gRatio->SetMarkerStyle(21);
    gRatio->SetTitle("");
    gRatio->Draw("AP");
    gRatio->GetXaxis()->SetTitle(variableString);

    constantFit->Draw("same");

    TLegend *l;
    CreateLegend(&l,0.55, 0.9, 0.9, 0.8);
    l->AddEntry(gRatio, "data/trueMC", "lep");
    l->AddEntry(constantFit, "constant fit: " + convertToString(systematicError ) + "+/-" + convertToString( systematicErrorUncertainty ),"ple");
    l->SetTextSize(0.04);
    l->Draw("same");

    outFile->cd();

    c1->Write(TString::Format("hSysErrorFit%d", switcher)); // Save the canvas as an image (optional)
    cout << "Just saved canvas with sys error as " << TString::Format("hSysErrorFit%d", switcher) << endl;
    return systematicErrorUncertainty;
}


TGraphAsymmErrors* PlotTofEffMult::getTrueMC(int switcher){

	effFile = new TFile("trueMC.root", "r");

	if(!effFile){
		cout << "Could not get trueMC.root." << endl;
		return nullptr;
	}
	TString variable, variableString;
    if( switcher == 1){//eta
    	variable = "eta";
    	variableString = "#eta";
    } else if( switcher == 2){
    	variable = "phi";
    	variableString = "#phi";
    } else if(switcher == 3){
    	variable = "pTInGev";
    	variableString = "p_{T}";
    }else if(switcher == 4){
    	variable = "vertexZInCm";
    }

    // load TPC or TOF plots and calulcate eff hist
    auto effhist = [&](TString dir, TString val, TString suffixTotal,TString suffixPassed, int part)
    {
       // Retrieve the canvas from the file
       TH1D* hist = nullptr;

       TCanvas *canvas = (TCanvas*)effFile->Get(dir + "/" + val +suffixPassed);
       if (!canvas) {
          std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first canvas!" << std::endl;
          return hist;
       }

       TH1D *hPassed = (TH1D*)canvas->GetPrimitive("h");
       if (!hPassed) {
          std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving first histogram!" << std::endl;
          return hist;
       }

       canvas = (TCanvas*)effFile->Get(dir + "/" + val  + suffixTotal);
       if (!canvas) {
          std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving second canvas!" << std::endl;
          return hist;
       }

       TH1D *hTotal = (TH1D*)canvas->GetPrimitive("h2");
       if (!hTotal) {
          std::cerr << "Error in PlotManager::loadEfficiencies() while retrieving second histogram!" << std::endl;
          return hist;
       }

       hPassed->Divide(hTotal);
       return hPassed;
    };

    TH1D *hTrueMC = effhist("Embedding/TofEff/pionPlus", variable, "","_Tof", 0);
    TH1D *hEffMinus = effhist("Embedding/TofEff/pionMinus", variable, "","_Tof", 1);
    
    hTrueMC->Add(hEffMinus);
    hTrueMC->Scale(0.5);

    TGraphAsymmErrors* graph = new TGraphAsymmErrors(hTrueMC->GetNbinsX());

    if(switcher == 5){ //inak 2, teraz vypnute
	    for (int i = 1; i <= hTrueMC->GetNbinsX(); ++i) {
	        // Bin propertie
	        double x = hTrueMC->GetBinCenter(i);              // Center of the bin
	        double y = hTrueMC->GetBinContent(i);             // Bin content
	        double exl = hTrueMC->GetBinCenter(i) - hTrueMC->GetBinLowEdge(i);  // Left error in X
	        double exh = hTrueMC->GetBinWidth(i) - exl;       // Right error in X
	        double eyl = hTrueMC->GetBinError(i);             // Y error (assuming symmetric here)
	        double eyh = eyl;                              // Set as symmetric Y error

	        // Set the point and its errors in the graph
	        graph->SetPoint(i - 1, x*180/TMath::Pi(), y);                  // TGraph indices start at 0
	        graph->SetPointError(i - 1, exl*180/TMath::Pi(), exh*180/TMath::Pi(), 0, 0);
	    }
	}else{
	    for (int i = 1; i <= hTrueMC->GetNbinsX(); ++i) {
	        // Bin propertie
	        double x = hTrueMC->GetBinCenter(i);              // Center of the bin
	        double y = hTrueMC->GetBinContent(i);             // Bin content
	        double exl = hTrueMC->GetBinCenter(i) - hTrueMC->GetBinLowEdge(i);  // Left error in X
	        double exh = hTrueMC->GetBinWidth(i) - exl;       // Right error in X
	        double eyl = hTrueMC->GetBinError(i);             // Y error (assuming symmetric here)
	        double eyh = eyl;                              // Set as symmetric Y error

	        // Set the point and its errors in the graph
	        graph->SetPoint(i - 1, x, y);                  // TGraph indices start at 0
	        graph->SetPointError(i - 1, exl, exh, 0, 0);
	    }
	}
    graph->SetLineColor(kBlack);
    graph->SetMarkerSize(3);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(20);
    graph->GetXaxis()->SetRangeUser(-1,1);

    return graph;

}

