#ifndef PlotAnalysis_h
#define PlotAnalysis_h


#include "Plot.h"
#include "Fit.h"
#include "Config.h"


using namespace std;

class PlotAnalysis : public Plot {
	public:
		//accessible from outside the class
		PlotAnalysis(const string mInputList, const char* filePath);
		PlotAnalysis(const string mInputList, shared_ptr<TFile> file);
      	~PlotAnalysis(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;
		void Finish() override;



    private:


		
};


#endif