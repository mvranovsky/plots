#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <bitset>
#include "TString.h"




const bool DEBUGMODE = true; // true = print debug messages, false = normal mode

//const TString plotDescription = "STAR Internal";
const TString plotDescription = "THIS THESIS";

const TString MCFileName = "slight_EtaCut1_100k.root";
const TString MCTreeName = "T";


const TString goodRunListWithRp = "/star/u/mvranovsk/star-upcDst/work/lists/goodRun17WithRP.list";

const TString goodRunListNoRp = "/star/u/mvranovsk/star-upcDst/work/lists/goodRun17noRP.list";



// settings for various analyses

const TString nameOfAnaBPTree = "recTreeBP";
const TString nameOfAnaBPDir = "AnaBPPlots";

const TString nameOfAnaV0Tree = "recTreeAnaV0";
const TString nameOfAnaV0Dir = "AnaV0Plots";

const TString nameOfAnaV0MultTree = "recTreeAnaV0Mult";
const TString nameOfAnaV0MultDir = "AnaV0MultPlots";

const TString nameOfTofEffTree = "recTreeTofEff";
const TString nameOfTofEffDir = "TofEffPlots";

const TString nameOfTofEffMultTree = "recTreeTofEffMult";
const TString nameOfTofEffMultDir = "TofEffMultPlots";

const TString nameOfAnaJPsiTree = "recTreeAnaJPsi";
const TString nameOfAnaJPsiDir = "AnaJPsiPlots";

const TString nameOfAnaJPSITree = "recTreeAnaJPSI";
const TString nameOfAnaJPSIDir = "AnaJPSIPlots";

const TString nameOfEmbeddingJPsiTree = "recTreeEmbeddingJPsi";
const TString nameOfEmbeddingJPsiDir = "EmbeddingJPsiPlots";

const TString nameOfAnaGoodRunTree = "recAnaGoodRun";
const TString nameOfAnaGoodRunDir = "AnaGoodRunPlots";

const TString nameOfAnaZeroBiasTree = "recTreeZeroBias";
const TString nameOfAnaZeroBiasDir = "AnaZeroBiasPlots";

const TString nameOfBemcEfficiencyTree = "recTreeBemcEfficiency";
const TString nameOfBemcEfficiencyDir = "BemcEfficiencyPlots";

const TString nameOfAnaChiCTree = "recTreeAnaChiC";
const TString nameOfAnaChiCDir = "AnaChiCPlots";




#endif //ifndef CONFIG_H