#ifndef SysStudy_h
#define SysStudy_h

#include "Util.h"
#include "Plot.h"
#include "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h"

using namespace RooFit;
using namespace std;
using namespace UTIL;

class SysStudy : public Plot {
	public:
		//accessible from outside the class
		SysStudy(const string mInputList, const char* filePath);
      	~SysStudy(){};
      
      	//overrides the original Make() from class Plot
      	void Make() override;
      	void Init() override;


    private:
    	//store variables, which can be accessed only from within this class

		TString convertToString(double val);
		int makeInt(double val);

};


#endif