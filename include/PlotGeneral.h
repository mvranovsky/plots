#ifndef PlotGeneral_h
#define PlotGeneral_h


#include "Util.h"
#include "Plot.h"


using namespace std;
using namespace UTIL;

class PlotGeneral : public Plot {
	public:
		//accessible from outside the class
		PlotGeneral(TFile *mOutFile, const string mInputList, const char* filePath);
      	~PlotGeneral(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

		
		vector<pair<TH1D*,TString>> hists1D;
		vector<pair<TH2F*,TString>> hists2D;

};


#endif