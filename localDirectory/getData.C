#include <iostream>
#include <cmath>


using namespace std;

TFile *file;


vector<TString> ListTCanvasInFile() {
	vector<TString> canvasName;
   	//TFile file(fileName);
    TList* keys = file->GetListOfKeys();
    //for (int i = 0; i < sizeof(keys); ++i){
        //cout << "Canvas: " << keys[i].GetName() << endl;
    //}
    if (keys) {
        TIter next(keys);
        TKey *key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TCanvas")) {
                TCanvas* canvas = (TCanvas*)obj;
                canvasName.push_back( canvas->GetName() );
                cout << "good "<< canvas->GetName() <<endl;
            }
            //if(obj){
                //delete obj;
            //}
        }
    }

    return canvasName;
}

void ListTDirectoriesInFile() {
    vector<TString> canvasName;
    //TFile file(fileName);
    TList* keys = file->GetListOfKeys();
    if (keys) {
        TIter next(keys);
        TKey *key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TDirectory")) {
                TDirectory* canvas = (TDirectory*)obj;
                canvasName.push_back( canvas->GetName() );
                cout << "Directory: " << obj->GetName() << endl;

            }
            delete obj;

        }
    }

    //return canvasName;
}


void getData(const TString& inputTag) {
	
	TString inputFile = inputTag + TString("/AnalysisOutput.root");


   	file = new TFile(inputFile);
   	if(!file){
   		cerr << "Couldn't open the input file. Leaving." << endl;
   		return;
   	}
   	cout << "Opened input file." << endl;
    //ListTDirectoriesInFile();
   	vector<TString> canvasList = ListTCanvasInFile();
   	cout << "Number of canvases: " << canvasList.size() << endl;
   	if(canvasList.size() == 0){
   		cerr << "No canvases. Leaving." << endl;
   		return;
   	}

   	TCanvas *c;
   	for (int i = 0; i < canvasList.size(); ++i) {
   		
   		c = dynamic_cast<TCanvas*>(file->Get( canvasList[i] ) ); //+ TString(";1") //the ;1 is there because there are histograms saved with the same name but with ";2"
   		if(!c) {
   			cerr << "Couldn't get canvas with the name " << canvasList[i] << ". Skipping." << endl;
   			continue;
   		}

   		c->SaveAs(inputTag + TString("/") + canvasList[i] + TString(".pdf") );
   	}
   	cout << "All canvases saved. Leaving." <<endl;

//additionally, can work with histograms in RooFit
  	

}