#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>


using namespace std;


void mergeData(){

	TFile *outputFile = new TFile("HistsForTomas.root", "recreate");


	TFile *fileMC = new TFile("AnaV0_noPV_MC/AnalysisOutput.root", "read");

	TFile *file = new TFile("AnaV0_noPV/AnalysisOutput.root", "read");


	//----------------------------------------------------------------------------------------

	outputFile->mkdir("efficiencyMC");
	outputFile->cd("efficiencyMC");

	TGraphAsymmErrors *g1 = (TGraphAsymmErrors*)fileMC->Get("effInteta");
	g1->Write("effInteta");

	g1 = (TGraphAsymmErrors*)fileMC->Get("effIntphi");
	g1->Write("effIntphi");

	g1 = (TGraphAsymmErrors*)fileMC->Get("effIntpT");
	g1->Write("effIntpT");

	g1 = (TGraphAsymmErrors*)fileMC->Get("effIntVz");
	g1->Write("effIntVz");

	//---------------------------------------------------------------------------------------

    outputFile->cd();
    outputFile->mkdir("efficiencyPLUS");
    outputFile->cd("efficiencyPLUS");

   	g1 = (TGraphAsymmErrors*)file->Get("PLUS/effInteta");
	g1->Write("effInteta");

	g1 = (TGraphAsymmErrors*)file->Get("PLUS/effIntphi");
	g1->Write("effIntphi");

	g1 = (TGraphAsymmErrors*)file->Get("PLUS/effIntpT");
	g1->Write("effIntpT");

	g1 = (TGraphAsymmErrors*)file->Get("PLUS/effIntVz");
	g1->Write("effIntVz");

	//---------------------------------------------------------------------------------------

    outputFile->cd();
    outputFile->mkdir("efficiencyMINUS");
    outputFile->cd("efficiencyMINUS");

    g1 = (TGraphAsymmErrors*)file->Get("MINUS/effInteta");
	g1->Write("effInteta");

	g1 = (TGraphAsymmErrors*)file->Get("MINUS/effIntphi");
	g1->Write("effIntphi");

	g1 = (TGraphAsymmErrors*)file->Get("MINUS/effIntpT");
	g1->Write("effIntpT");

	g1 = (TGraphAsymmErrors*)file->Get("MINUS/effIntVz");
	g1->Write("effIntVz");

	//--------------------------------------------------------------------------------------

    outputFile->cd();
    outputFile->mkdir("efficiencyFULL");
    outputFile->cd("efficiencyFULL");

    g1 = (TGraphAsymmErrors*)file->Get("effInteta");

	g1->Write("effInteta");

	g1 = (TGraphAsymmErrors*)file->Get("effIntphi");

	g1->Write("effIntphi");

	g1 = (TGraphAsymmErrors*)file->Get("effIntpT");

	g1->Write("effIntpT");

	g1 = (TGraphAsymmErrors*)file->Get("effIntVz");

	g1->Write("effIntVz");

	//--------------------------------------------------------------------------------------

	fileMC->Close();
	file->Close();
	outputFile->Close();
	cout << "All plots stored in HistsForTomas.root. Goodbye. " << endl;

}