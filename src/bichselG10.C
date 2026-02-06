//
// ak chces zmenit na ine osi, musis zmenit co plotujes: whatToDraw, xmax, xmin pre funkcie, bichselZ zmenit pove a zrusit prvy Log10 v returne
//
//
#include "../include/Config.h"



//root4star [0] .L bichselG10.C 
//root4star [1] .x bichselG10.C("/gpfs01/star/pwg/truhlar/Run17_P20ic/mainAna0304/merged/StRP_production_0000.root") 
//
#if !defined(__CINT__)
// code that should be seen ONLY by the compiler
#else
#if !defined(__CINT__) || defined(__MAKECINT__)
// code that should be seen by the compiler AND rootcint
#else
// code that should always be seen
#endif
#endif
//________________________________________________________________________________
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include <stdio.h>
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TClassTable.h"
#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"

#else
class Bichsel;
#endif


// 0. value is the default one, 1. is the tight condition, 2. is the loose condition
Int_t nHitsDedx[3] = {15, 17, 12};
TString nHitsDedxDescription[3] = {Form("Nominal (N^{dEdx}_{hits} < %d)", nHitsDedx[0]), Form("Tight (N^{dEdx}_{hits} < %d)", nHitsDedx[1]), Form("Loose (N^{dEdx}_{hits} < %d)", nHitsDedx[2])};
Int_t nHitsFit[3] = {15, 17, 12};
TString nHitsFitDescription[3] = {Form("Nominal (N^{fit}_{hits} < %d)", nHitsFit[0]), Form("Tight (N^{fit}_{hits} < %d)", nHitsFit[1]), Form("Loose (N^{fit}_{hits} < %d)", nHitsFit[2])};
Double_t dcaZInCm[3] = {1, 0.7, 1.2};
TString dcaZInCmDescription[3] = {Form("Nominal (DCA_{Z} < %.1f)", dcaZInCm[0]), Form("Tight (DCA_{Z} < %.1f)", dcaZInCm[1]), Form("Loose (DCA_{Z} < %.1f)", dcaZInCm[2])};
Double_t dcaXYInCm[3] = {1.5, 1.2, 1.8};
TString dcaXYInCmDescription[3] = {Form("Nominal (DCA_{XY} < %.1f)", dcaXYInCm[0]), Form("Tight (DCA_{XY} < %.1f)", dcaXYInCm[1]), Form("Loose (DCA_{XY} < %.1f)", dcaXYInCm[2])};

Int_t chiSquareE[3] = {9, 7, 11};
TString chiSquareEDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareE[0]), Form("Tight (chi^{2} < %d)", chiSquareE[1]), Form("Loose (chi^{2} < %d)", chiSquareE[2])};
Int_t chiSquarePi[3] = {10, 12, 8};
TString chiSquarePiDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquarePi[0]), Form("Tight (chi^{2} < %d)", chiSquarePi[1]), Form("Loose (chi^{2} < %d)", chiSquarePi[2])};
Int_t chiSquareK[3] = {10, 12, 8};
TString chiSquareKDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareK[0]), Form("Tight (chi^{2} < %d)", chiSquareK[1]), Form("Loose (chi^{2} < %d)", chiSquareK[2])};
Int_t chiSquareP[3] = {10,12,8};
TString chiSquarePDescription[3] = {Form("Nominal (chi^{2} < %d)", chiSquareP[0]), Form("Tight (chi^{2} < %d)", chiSquareP[1]), Form("Loose (chi^{2} < %d)", chiSquareP[2])};

Int_t vertexZInCm[3] = {100, 80, 120};
TString vertexZInCmDescription[3] = {Form("Nominal (|V_{Z}| < %d)", vertexZInCm[0]), Form("Tight (|V_{Z}| < %d)", vertexZInCm[1]), Form("Loose (|V_{Z}| < %d)", vertexZInCm[2])};
Double_t etaHadron[3] = {0.9, 0.8, 1.0};
TString etaHadronDescription[3] = {Form("Nominal (|#eta_{e}| < %.1f)", etaHadron[0]), Form("Tight (|#eta_{e}| < %.1f)", etaHadron[1]), Form("Loose (|#eta_{e}| < %.1f)", etaHadron[2])};

double nSigmaH = 1.3; //used for protons and kaons in addition to electrons


TFile *outFile;

Bichsel *m_Bichsel = 0;
const Int_t NMasses = 10;
const Double_t Masses[NMasses] = {0.13956995,
               0.493677,
               0.93827231,
               1.87561339,
               0.51099907e-3,
               0.1056584,
               2.80925,
               2.80923, //GEANT3
               3.727417, //GEANT3
               0.13956995,
};
const Int_t   Index[NMasses] = { 2,    3,   4,   5,   0,    1,  6,    7,       8,    -2};
//const Int_t   Index[NMasses] = { 4,    3,   2,   0,   5,    1,  6,    7,       8,    -2};
const Int_t   Colors[NMasses] = { 1,    2,   4,   7,   6,    3,  9,    30,       8,    -2};
const Char_t *Names[NMasses] = {"Pion","Kaon","Proton","Deuteron","Electron","#mu","t","He3","#alpha","2#pi"};
//const Char_t *Names[NMasses] = {"Proton", "Kaon","Pion","Electron", "Deuteron","#mu","t","He3","#alpha","2#pi"};
const Int_t NF = 5;  //         0       1    2     3     4      5   6     7
const Char_t *FNames[8] = {"Girrf","Sirrf","Bz","B70","B60","B70M","dNdx","BzM"};
const Int_t Nlog2dx = 3;
const Double_t log2dx[Nlog2dx] = {0,1,2};
//________________________________________________________________________________
Double_t bichselZ(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   //Double_t pove = x[0];
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZ(TMath::Log10(pove) - TMath::Log10(mass) ,dx2)));//TMath::Exp(7.81779499999999961e-01));

   //return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}
//________________________________________________________________________________
Double_t bichselZM(Double_t *x,Double_t *par) {
   
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  (scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZM(TMath::Log10(poverm),dx2)));//TMath::Exp(7.81779499999999961e-01));
   
}
//________________________________________________________________________________
Double_t bichsel70(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   // return  TMath::Log10(scale*charge*charge*m_Bichsel->GetI70M(TMath::Log10(poverm),dx2));//TMath::Exp(7.81779499999999961e-01));
   return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}
//________________________________________________________________________________
Double_t bichsel70M(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*m_Bichsel->GetI70M(TMath::Log10(poverm),dx2));//TMath::Exp(7.81779499999999961e-01));
}
//________________________________________________________________________________
Double_t dNdx(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*StdEdxModel::instance()->dNdx(poverm,charge));//TMath::Exp(7.81779499999999961e-01));
}
#ifndef __CINT__
//________________________________________________________________________________
Double_t aleph70P(Double_t *x,Double_t *par) {
   /* 
       W.Blum, L. Rolandi "Particle Detection with Drift Chambers", page 246, eq. (9.5)
       F_g(v) = p[0]/beta**p[3]*(p[1] - beta**p[3] - log(p[2] + (beta*gamma)**-p[4]);
       F_g(v) = p[0]*(1/beta**p[3]*(p[1] - log(p[2] + 1/(beta*gamma)**p[4])) - 1) 
   */
   Double_t bg = x[0];
   Double_t b2inv = 1. + 1./(bg*bg);
   Double_t beta  = 1./TMath::Sqrt(b2inv);
   Double_t dEdx = par[0]*(-1 + TMath::Power(beta,-par[3])*(par[1] - TMath::Log(TMath::Max(1e-10,par[2] + TMath::Power(bg,-par[4])))));
   return dEdx;
};
//________________________________________________________________________________
Double_t aleph70(Double_t *x,Double_t *par) {
   static const Double_t dEdxMIP = 2.39761562607903311;
   static Double_t MIPBetaGamma = 4.;
#if 0
   struct Par_t {Int_t h, N; Double_t xmin, xmax, pars[10];};
   const Par_t Par[9] = {
      /* name          h n+1   xmin   xmax      pars[n+1] */
      /* electron */{  0,  4,   3.0,   6.0,{     0.14105,    -0.09078,     0.01901,    -0.00128,           0,           0,           0,           0,           0,           0}},
      /*     muon */{  1,  2,   0.0,   4.5,{    -0.00689,     0.00256,           0,           0,           0,           0,           0,           0,           0,           0}},
      /*     pion */{  2, -1,   0.0,   4.5,{     0.00000,           0,           0,           0,           0,           0,           0,           0,           0,           0}},
      /*     kaon */{  3,  9,  -0.1,   3.7,{     0.00869,     0.21918,    -0.88919,     1.30023,    -0.97075,     0.41298,    -0.10214,     0.01379,    -0.00079,           0}},
      /*   proton */{  4,  9,  -0.6,   3.3,{     0.03052,    -0.02423,    -0.05636,    -0.11585,     0.41292,    -0.38956,     0.16837,    -0.03494,     0.00283,           0}},
      /* deuteron */{  5,  8,  -1.0,   2.9,{     0.03523,    -0.10625,     0.04182,     0.07820,    -0.03816,    -0.02735,     0.01940,    -0.00304,           0,           0}},
      /*   triton */{  6, 10,  -1.0,   2.8,{     0.03092,    -0.07846,     0.01668,     0.00331,     0.12771,    -0.05480,    -0.13897,     0.13928,    -0.04715,     0.00555}},
      /*      He3 */{  7,  9,  -0.8,   2.9,{     0.09158,    -0.07816,     0.07039,     0.00578,    -0.16160,     0.26547,    -0.18728,     0.05988,    -0.00710,           0}},
      /*    alpha */{  8, 10,  -0.8,   2.9,{     0.09366,    -0.08276,     0.06191,     0.02631,    -0.17044,     0.30536,    -0.26867,     0.11847,    -0.02505,     0.00201}}
   };
   //  const Double_t ppar[7]     = { 0.0857988,   9.46138, 0.000206655,     2.12222,       0.974,    -1, 0.13957}; /* pion */
   static Double_t ppar[7]     = { 0.0857988,   9.46138, 0.000206655,     2.12222,       0.974,    -1, 0.13957}; /* pion */
#else
   //  static Double_t ppar[7]     = { 0.0762,  10.632, 0.134e-4, 1.863,  1.948,    -1, -1}; /* Aleph parameters from Alice TPC TDR */
   //  static Double_t ppar[7] = { 0.08942,     8.91971,     0.00024,     2.27383,     1.54174,    -1.00000, 0}; /* pion */
   static Double_t ppar[7]     = {0.12337,     6.61371,     0.00201,     2.27381,     1.54174,    -1.00000, 0 }; /* g All */
#endif
   static Double_t Norm = dEdxMIP/aleph70P(&MIPBetaGamma,ppar);
   Int_t hyp = (Int_t ) par[0];
   Int_t h = Index[hyp];
   Double_t ScaleL10 = 0;
   if (h < 0) {
      h = -h;
      ScaleL10 = TMath::Log10(2.);
   }
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t mass = Masses[hyp];
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   if (h > 6) {
      charge = 2;
      poverm *= charge;
   }
   Double_t bg = poverm;
   /* 
       W.Blum, L. Rolandi "Particle Detection with Drift Chambers", page 246, eq. (9.5)
       F_g(v) = p[0]/beta**p[3]*(p[1] - beta**p[3] - log(p[2] + (beta*gamma)**-p[4]));
       F_g(v) = p[0]*(1/beta**p[3]*(p[1] - log(p[2] + 1/(beta*gamma)**p[4])) - 1) 
   */
   Double_t dEdxL10 =  TMath::Log10(Norm*aleph70P(&bg,ppar));
#if 0
   if (Par[h].N > 0) {
      TString fName(Form("pol%i",Par[h].N-1));
      TF1 *f = (TF1 *) gROOT->GetListOfFunctions()->FindObject(fName);
      if (! f) {
         f = new TF1(fName,fName,0,1);
      }
      f->SetParameters(&Par[h].pars[0]);
      Double_t bgL10 = TMath::Log10(bg);
      dEdxL10 += f->Eval(bgL10);
   }
#endif
   return 2*TMath::Log10(charge) + dEdxL10 + ScaleL10;
}
#endif /* __CINT__ */

TString getCondition(TString var = "", int j = 0 ){

   if(var.Contains("tight") || var.Contains("TIGHT") || var.Contains("Tight") ){
      j = 1;
   }else if(var.Contains("loose") || var.Contains("LOOSE") || var.Contains("Loose")){
      j = 2;
   }

   TString c = "";
   if(var.Contains("nHitsFit") || var.Contains("NHITSFIT") || var.Contains("nhitsfit")){
      c += Form("nHitsFit0 > %d && nHitsFit1 > %d && ", nHitsFit[j], nHitsFit[j]);
   }else{
      c += Form("nHitsFit0 > %d && nHitsFit1 > %d && ", nHitsFit[0], nHitsFit[0]);
   }

   if(var.Contains("nHitsDedx") || var.Contains("NHITSDEDX") || var.Contains("nhitsdedx")){
      c += Form("nHitsDEdx0 > %d && nHitsDEdx1 > %d && ", nHitsDedx[j], nHitsDedx[j]);
   }else{
      c += Form("nHitsDEdx0 > %d && nHitsDEdx1 > %d && ", nHitsDedx[0], nHitsDedx[0]);
   }

   if(var.Contains("etaHadron") || var.Contains("etahadron") || var.Contains("ETAHADRON") || var.Contains("eta")){
      c += Form("abs(etaHadron0) < %.1f && abs(etaHadron1) < %.1f && ", etaHadron[j],etaHadron[j]);
   }else{
      c += Form("abs(etaHadron0) < %.1f && abs(etaHadron1) < %.1f && ", etaHadron[0], etaHadron[0]);
   }

   if(var.Contains("chiSquare") || var.Contains("PID") || var.Contains("chiSquareEE") || var.Contains("chiSquareE")){
      c += Form("chiSquareelectron < %.d && ", chiSquareE[j]);
   }else{
      c += Form("chiSquareelectron < %.d && ", chiSquareE[0]);
   }

   if(var.Contains("nSigmaProton") || var.Contains("NSIGMAPROTON") || var.Contains("nsigmaproton")){
      c+= Form("abs(nSigmaTPCprotonPlus) > %.1f && abs(nSigmaTPCprotonMinus) > %.1f && ", nSigmaH, nSigmaH);
   }

   if(var.Contains("nSigmaKaon") || var.Contains("NSIGMAKAON") || var.Contains("nsigmakaon")){
      c+= Form("abs(nSigmaTPCkaonPlus) > %.1f && abs(nSigmaTPCkaonMinus) > %.1f && ", nSigmaH, nSigmaH);
   }

   
   
   if( !(var.Contains("embedding") || var.Contains("Embedding") || var.Contains("EMBEDDING")) ){  //embedding is not able to reconstruct vertex 
      
      //c += Form("chiSquarepion > %.d && chiSquareproton > %.d && chiSquarekaon > %.d && ", chiSquarePi[0], chiSquareP[0], chiSquareK[0]);
      
      if(var.Contains("dcaZInCm") || var.Contains("dcazincm") || var.Contains("DCAZINCM")){
         c += Form("abs(dcaZInCm0) < %.1f && abs(dcaZInCm1) < %.1f && ", dcaZInCm[j], dcaZInCm[j]);
      }else{
         c += Form("abs(dcaZInCm0) < %.1f && abs(dcaZInCm1) < %.1f && ", dcaZInCm[0], dcaZInCm[0]);
      }

      if(var.Contains("dcaXYInCm") || var.Contains("dcaxyincm") || var.Contains("DCAXYINCM")){
         c += Form("dcaXYInCm0 < %.1f && dcaXYInCm1 < %.1f && ", dcaXYInCm[j], dcaXYInCm[j]);
      }else{
         c += Form("dcaXYInCm0 < %.1f && dcaXYInCm1 < %.1f && ", dcaXYInCm[0], dcaXYInCm[0]);
      }

      if(var.Contains("vertexZInCm") || var.Contains("vertexzincm") || var.Contains("VERTEXZINCM") || var.Contains("Vz")){
         c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[j]);
      }else{
         c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[0]);
      }
   }

   //c += Form("isBemcHit0 == 1 && isBemcHit1 == 1 && "); // require both electrons to have matched BEMC hit

   //remove the last " && "
   if(c.EndsWith("&& ") || c.EndsWith(" &&")){
      c.Remove(c.Length() - 3, 3);
   }

   // cout << "Condition for variable " << var << ": " << c << " and j = " << j << endl;

   return c;

}


TString getCondition1Track(int i,TString var = "", int j = 0 ){

   if(var.Contains("tight") || var.Contains("TIGHT") || var.Contains("Tight") ){
      j = 1;
   }else if(var.Contains("loose") || var.Contains("LOOSE") || var.Contains("Loose")){
      j = 2;
   }

   TString c = "";
   if(i == 0){
      
      if(var.Contains("nHitsFit") || var.Contains("NHITSFIT") || var.Contains("nhitsfit")){
         c += Form("nHitsFit0 > %d && ", nHitsFit[j]);
      }else{
         c += Form("nHitsFit0 > %d && ", nHitsFit[0]);
      }

      if(var.Contains("nHitsDedx") || var.Contains("NHITSDEDX") || var.Contains("nhitsdedx")){
         c += Form("nHitsDEdx0 > %d && ", nHitsDedx[j]);
      }else{
         c += Form("nHitsDEdx0 > %d && ", nHitsDedx[0]);
      }

      if(var.Contains("etaHadron") || var.Contains("etahadron") || var.Contains("ETAHADRON") || var.Contains("eta")){
         c += Form("abs(etaHadron0) < %.1f && ", etaHadron[j]);
      }else{
         c += Form("abs(etaHadron0) < %.1f && ", etaHadron[0]);
      }

      if(var.Contains("chiSquare") || var.Contains("PID") || var.Contains("chiSquareEE") || var.Contains("chiSquareE")){
         c += Form("chiSquareelectron < %.d && ", chiSquareE[j]);
      }else{
         c += Form("chiSquareelectron < %.d && ", chiSquareE[0]);
      }

      if(var.Contains("nSigmaProton") || var.Contains("NSIGMAPROTON") || var.Contains("nsigmaproton")){
         c+= Form("abs(nSigmaTPCprotonPlus) > %.1f && ", nSigmaH);
      }

      if(var.Contains("nSigmaKaon") || var.Contains("NSIGMAKAON") || var.Contains("nsigmakaon")){
         c+= Form("abs(nSigmaTPCkaonPlus) > %.1f && ", nSigmaH);
      }

      
      
      if( !(var.Contains("embedding") || var.Contains("Embedding") || var.Contains("EMBEDDING")) ){  //embedding is not able to reconstruct vertex 
         
         //c += Form("chiSquarepion > %.d && chiSquareproton > %.d && chiSquarekaon > %.d && ", chiSquarePi[0], chiSquareP[0], chiSquareK[0]);
         
         if(var.Contains("dcaZInCm") || var.Contains("dcazincm") || var.Contains("DCAZINCM")){
            c += Form("abs(dcaZInCm0) < %.1f && ", dcaZInCm[j]);
         }else{
            c += Form("abs(dcaZInCm0) < %.1f && ", dcaZInCm[0]);
         }

         if(var.Contains("dcaXYInCm") || var.Contains("dcaxyincm") || var.Contains("DCAXYINCM")){
            c += Form("dcaXYInCm0 < %.1f && ", dcaXYInCm[j]);
         }else{
            c += Form("dcaXYInCm0 < %.1f && ", dcaXYInCm[0]);
         }

         if(var.Contains("vertexZInCm") || var.Contains("vertexzincm") || var.Contains("VERTEXZINCM") || var.Contains("Vz")){
            c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[j]);
         }else{
            c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[0]);
         }
      }
   }else{
      if(var.Contains("nHitsFit") || var.Contains("NHITSFIT") || var.Contains("nhitsfit")){
         c += Form("nHitsFit1 > %d && ", nHitsFit[j]);
      }else{
         c += Form("nHitsFit1 > %d && ", nHitsFit[0]);
      }

      if(var.Contains("nHitsDedx") || var.Contains("NHITSDEDX") || var.Contains("nhitsdedx")){
         c += Form("nHitsDEdx1 > %d && ", nHitsDedx[j]);
      }else{
         c += Form("nHitsDEdx1 > %d && ", nHitsDedx[0]);
      }

      if(var.Contains("etaHadron") || var.Contains("etahadron") || var.Contains("ETAHADRON") || var.Contains("eta")){
         c += Form("abs(etaHadron1) < %.1f && ", etaHadron[j]);
      }else{
         c += Form("abs(etaHadron1) < %.1f && ", etaHadron[0]);
      }

      if(var.Contains("chiSquare") || var.Contains("PID") || var.Contains("chiSquareEE") || var.Contains("chiSquareE")){
         c += Form("chiSquareelectron < %.d && ", chiSquareE[j]);
      }else{
         c += Form("chiSquareelectron < %.d && ", chiSquareE[0]);
      }

      if(var.Contains("nSigmaProton") || var.Contains("NSIGMAPROTON") || var.Contains("nsigmaproton")){
         c+= Form("abs(nSigmaTPCprotonMinus) > %.1f && ", nSigmaH);
      }

      if(var.Contains("nSigmaKaon") || var.Contains("NSIGMAKAON") || var.Contains("nsigmakaon")){
         c+= Form("abs(nSigmaTPCkaonMinus) > %.1f && ", nSigmaH);
      }

      
      
      if( !(var.Contains("embedding") || var.Contains("Embedding") || var.Contains("EMBEDDING")) ){  //embedding is not able to reconstruct vertex 
         
         //c += Form("chiSquarepion > %.d && chiSquareproton > %.d && chiSquarekaon > %.d && ", chiSquarePi[0], chiSquareP[0], chiSquareK[0]);
         
         if(var.Contains("dcaZInCm") || var.Contains("dcazincm") || var.Contains("DCAZINCM")){
            c += Form("abs(dcaZInCm1) < %.1f && ", dcaZInCm[j]);
         }else{
            c += Form("abs(dcaZInCm1) < %.1f && ", dcaZInCm[0]);
         }

         if(var.Contains("dcaXYInCm") || var.Contains("dcaxyincm") || var.Contains("DCAXYINCM")){
            c += Form("dcaXYInCm1 < %.1f && ", dcaXYInCm[j]);
         }else{
            c += Form("dcaXYInCm1 < %.1f && ", dcaXYInCm[0]);
         }

         if(var.Contains("vertexZInCm") || var.Contains("vertexzincm") || var.Contains("VERTEXZINCM") || var.Contains("Vz")){
            c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[j]);
         }else{
            c += Form("abs(vertexZInCm) < %d && ", vertexZInCm[0]);
         }
      }
   }
   //c += Form("isBemcHit0 == 1 && isBemcHit1 == 1 && "); // require both electrons to have matched BEMC hit

   //remove the last " && "
   if(c.EndsWith("&& ") || c.EndsWith(" &&")){
      c.Remove(c.Length() - 3, 3);
   }

   // cout << "Condition for variable " << var << ": " << c << " and j = " << j << endl;

   return c;

}

//________________________________________________________________________________
TString defineAnalysis(TString inputPosition){

   if( inputPosition.Contains("AnaV0Mult") ){
      return nameOfAnaV0MultTree;
   }else if( inputPosition.Contains("AnaV0") ){
      return nameOfAnaV0Tree;
   } else if( inputPosition.Contains("TofEffMult") ){
      return nameOfTofEffMultTree;
   } else if( inputPosition.Contains("TofEff") ){
      return nameOfTofEffTree;
   } else if( inputPosition.Contains("AnaJPsi") ){
      return nameOfAnaJPsiTree;
   } else if( inputPosition.Contains("AnaJPSI") ){
      return nameOfAnaJPSITree;
   } else if( inputPosition.Contains("EmbeddingJPsi") ){
      return nameOfEmbeddingJPsiTree;
   } else if( inputPosition.Contains("AnaGoodRun") ){
      return nameOfAnaGoodRunTree;
   } else if( inputPosition.Contains("SysStudy") ){
      return nameOfAnaJPsiTree;
   } else if(inputPosition.Contains("AnaZeroBias")){
      return nameOfAnaZeroBiasTree;
   }else{
      cout << "No analysis to run. Leaving..." << endl;
      return "";
   } 
}


bool drawPIDPlot(TH2D* hdEdx, TString canvasName){

   TCanvas *cCanvas2D = new TCanvas(canvasName, canvasName, 800, 700);
   gPad->SetMargin(0.12,0.12,0.12,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);
   gStyle->SetOptDate(0);
   gStyle->SetLineWidth(2);      //axis line
   gStyle->SetFrameLineWidth(2); //frame line
   gPad->SetTickx();
   gPad->SetTicky(); 
   gStyle->SetOptStat("");
   cCanvas2D->SetGridx(0);
   cCanvas2D->SetGridy(0);
   hdEdx->GetXaxis()->SetTitleFont(22);
   hdEdx->GetYaxis()->SetTitleFont(22);
   hdEdx->GetZaxis()->SetTitleFont(22);
   hdEdx->GetXaxis()->SetLabelFont(22);
   hdEdx->GetYaxis()->SetLabelFont(22);
   hdEdx->GetZaxis()->SetLabelFont(22);

   //hdEdx->GetXaxis()->SetRangeUser(-0.8,1.0);
   //hdEdx->GetYaxis()->SetRangeUser(0.1,1.7);
   hdEdx->GetZaxis()->SetRangeUser(1.0, 50);
   hdEdx->GetXaxis()->SetLabelSize(0.05);
   hdEdx->GetYaxis()->SetLabelSize(0.05);
   hdEdx->GetZaxis()->SetLabelSize(0.05);
   hdEdx->GetXaxis()->SetTitleSize(0.05);
   hdEdx->GetYaxis()->SetTitleSize(0.05);
   hdEdx->GetZaxis()->SetTitleSize(0.05);
   hdEdx->GetXaxis()->SetTitleOffset(1.07);
   hdEdx->GetYaxis()->SetTitleOffset(1.0);
   hdEdx->GetZaxis()->SetTitleOffset(0.50);
   hdEdx->SetStats(0); 
   hdEdx->SetTitle(" ; log_{10} p [GeV/c] ;log_{10} dE/dx [keV/cm] ");
   cCanvas2D->SetLogz(1);
   //cCanvas2D->SetLogx(1);
   //cCanvas2D->SetLogy(1);
   hdEdx->Draw("COLZ");
   
   TPaveText *textSTAR;
   textSTAR = new TPaveText(0.6, 0.9, 0.85, 0.95,"brNDC");
   textSTAR -> SetTextSize(0.04);
   textSTAR -> SetFillStyle(0);
   textSTAR -> SetFillColorAlpha(kWhite,0);
   textSTAR -> SetBorderSize(0);
   textSTAR -> SetTextFont(72);
   textSTAR -> SetTextAlign(33);
   textSTAR->AddText(plotDescription);
   textSTAR -> Draw("same");

   TPaveText *textpp510;
   textpp510 = new TPaveText(0.6, 0.85, 0.85, 0.9,"brNDC");
   textpp510 -> SetTextSize(0.03);
   textpp510 -> SetFillStyle(0);
   textpp510 -> SetFillColorAlpha(kWhite,0);
   textpp510 -> SetBorderSize(0);
   textpp510 -> SetTextFont(62);
   textpp510 -> SetTextAlign(33);
   textpp510->AddText("p+p #sqrt{s} = 510 GeV");
   textpp510 -> Draw("same");

   const Char_t *type="Bz";

   TString Type(type);
   TLegend *leg = new TLegend(0.65,0.65,0.8,0.8);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetTextFont(22);
   Double_t xmax = 0.9;
   //Double_t xmax = 10;
   //  for (int h = 0; h < NMasses; h++) { // Masses
   for (int h = 0; h < NF; h++) 
   { // Masses
      Int_t f = 3;
      if      (Type.Contains("BzM",TString::kIgnoreCase))  f = 7;
      else if (Type.Contains("Bz",TString::kIgnoreCase))   f = 2;
      else if (Type.Contains("I70M",TString::kIgnoreCase)) f = 5;
      else if (Type.Contains("I70",TString::kIgnoreCase))  f = 3;
      else if (Type.Contains("I60",TString::kIgnoreCase))  f = 4;
      else if (Type.Contains("N",TString::kIgnoreCase))    f = 6;

      Int_t dx = 1;
      Char_t *FunName = Form("%s%s%i",FNames[f],Names[h],(int)log2dx[dx]);
      cout << "Make " << FunName << endl;
      //Double_t xmin = 0.1;
      Double_t xmin = -0.9;
      //    if (h == 0 || h >= 5) xmin = -0.75;
      if (h == 1) xmin = -0.80;//0.16;
      if (h == 2) xmin = -0.60;//0.25;
      if (h == 3) xmin = -0.30;//0.5;
      TF1 *func = 0;
      if      (f == 3) func = new TF1(FunName,bichsel70,xmin, xmax,2);
      else if (f == 2) func = new TF1(FunName,bichselZ ,xmin, xmax,2);
      else if (f == 5) func = new TF1(FunName,bichsel70M ,xmin, xmax,2);
      else if (f == 6) func = new TF1(FunName,dNdx ,xmin, xmax,2);
      else if (f == 7) func = new TF1(FunName,bichselZM,xmin, xmax,2);
      else { return;}

      if (h == 9) func->SetParameter(0,-Masses[h]);
      else       func->SetParameter(0,Masses[h]);

      func->SetParameter(1,1.);

      if (h >= 7 && h < 9) func->SetParameter(1,2.);

      Int_t color = Colors[h];
      #if 1
         func->SetLineColor(color);
         func->SetMarkerColor(color);
         func->SetLineStyle(1);
         func->SetLineWidth(4);
      #endif
      func->Draw("same");
      leg->AddEntry(func,Names[h], "l");
      #if !defined( __CINT__) && defined(__Aleph__)
         TF1 *fA = new TF1(Form("Aleph%s",Names[h]),aleph70,xmin,xmax, 1);
         fA->SetParameter(0,h);
         fA->SetLineColor(color);
         fA->SetMarkerColor(color);
         fA->SetLineStyle(4);
         fA->Draw("same");
         leg->AddEntry(fA,Names[h]);
      #endif
   }
   leg->Draw("same");

   cCanvas2D->Update();
   //cCanvas2D->SaveAs("dEdxAll.pdf");


   // Write histograms to output file
   outFile->cd();
   cCanvas2D->Write(canvasName, TObject::kOverwrite);
   return true;
}


void bichselG10(TString input) {  
   if (gClassTable->GetID("StBichsel") < 0 || !m_Bichsel)
   {
      gSystem->Load("libTable");
      gSystem->Load("St_base");
      gSystem->Load("StarClassLibrary");
      gSystem->Load("StBichsel");
      m_Bichsel = Bichsel::Instance();
   }

   TString nameOfTree = defineAnalysis(input);
   if(nameOfTree == ""){
      cout << "No analysis to run. Leaving..." << endl;
      return;
   }


   TString outFilePath = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + input + "/merged/AnalysisOutput.root";
   cout << "Outfile path: " << outFilePath << endl;

   outFile = TFile::Open(outFilePath, "UPDATE");
   if (!outFile || outFile->IsZombie()) {
      cout << "Error opening output file: " << outFilePath << endl;
      return;
   }

   cout << "About to run analysis whose tree is called: " << nameOfTree << endl;

   TString inputListLocation = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/" + input + "/merged/StRP_production.list";
   TChain *chain = new TChain(nameOfTree);
   TChain *chain_bcg = new TChain(nameOfTree + "_Bcg");
   ifstream instr(inputListLocation);
   if (!instr.is_open()){
      cout<< "Couldn't open: "<<inputListLocation<<endl;
      return;
   }

   int lineId = 0;
   string line;
   TFile *inputFile;
   while(getline(instr, line)) {

      if(line.empty())
         continue;

      inputFile = TFile::Open(line.c_str(), "read");
      if(!inputFile){
         cout << "Couldn't open: " << line.c_str() << endl;
         return;
      } 
      TTree *currentTree = dynamic_cast<TTree*>( inputFile->Get(nameOfTree) );
      if(currentTree){
         chain->AddFile(line.c_str());
         cout << "Successfully added file: " << line.c_str() << endl;
      } else {
         cout << "Name of tree: " << nameOfTree << endl;
         cout << "Couldn't open .root file with name: " << line.c_str() << endl;
      }  

      TTree *currentTree_bcg = dynamic_cast<TTree*>( inputFile->Get(nameOfTree + "_Bcg") );
      if(currentTree_bcg){
         chain_bcg->AddFile(line.c_str());
         cout << "Successfully added background file: " << line.c_str() << endl;
      } else {
         cout << "Name of tree: " << nameOfTree + "_Bcg" << endl;
         cout << "Couldn't open .root file with name: " << line.c_str() << endl;
      }
      lineId++;
      //nInputFiles++;
   }//while
   cout << ".root files added to chain" << endl;
   
   TTree* tree = chain;
   if(!tree){
      cout << "Couldn't get tree from chain" << endl;
      return;
   }
   cout << "Correctly loaded main tree" << endl;
   
   TTree* bcg = dynamic_cast<TTree*>( chain_bcg );
   bool isBackground = true;
   if(!bcg){
      cout << "Couldn't get background tree from chain" << endl;
      isBackground = false;
   }


   TString condition = getCondition();
   TString conditionBcg = getCondition();
   if( nameOfTree.Contains("mult") || nameOfTree.Contains("Mult") ){
      condition = "pairID == 0 && totQ == 0"; //momentumChargePlus or Minus
      conditionBcg = "pairID == 0 && totQ != 0"; //momentumChargePlus or Minus
   }

   TString whatToDraw = "TMath::Log10(dEdxInKevCm0):TMath::Log10(momentumInGev0)>>hDEdx(400,-0.8,1,400,0.1,1.7)";
   TString whatToDraw2 = "TMath::Log10(dEdxInKevCm1):TMath::Log10(momentumInGev1)>>+hDEdx";
   //TString whatToDraw = "dEdxInKevCm0:momentumInGev0>>hDEdx(400, 0.1, 10, 400, 2, 10 )";
   //TString whatToDraw2 = "dEdxInKevCm1:momentumInGev1>>+hDEdx";

   tree->Draw(whatToDraw, getCondition1Track(0)); //momentumChargePlus or Minus
   tree->Draw(whatToDraw2, getCondition1Track(1)); 
   cout << "Condition for just signal: " << getCondition1Track(0) << endl;



   TH2D* hdEdx = (TH2D*)gPad->GetPrimitive("hDEdx");
   if(!hdEdx || hdEdx->IsZombie() || hdEdx->GetEntries() == 0){
      cout << "Empty 2D histogram. Leaving..." << endl;
      return;
   }

   if(drawPIDPlot(hdEdx, "BichselCanvas_justSignal")){
      cout << "PID plot 2 drawn successfully." << endl;
   } else {
      cout << "Failed to draw PID plot." << endl;
      return;
   }

   if(isBackground){
      bcg->Draw(whatToDraw, conditionBcg); //momentumChargePlus or Minus
      bcg->Draw(whatToDraw2, conditionBcg);

      TH2D* hdEdx = (TH2D*)gPad->GetPrimitive("hDEdx");
      if(!hdEdx || hdEdx->IsZombie() || hdEdx->GetEntries() == 0){
         cout << "Empty 2D histogram. Leaving..." << endl;
         return;
      }

      if(drawPIDPlot(hdEdx, "BichselCanvas_justBackground")){
         cout << "PID plot 3 drawn successfully." << endl;
      } else {
         cout << "Failed to draw PID plot." << endl;
         return;
      }
   }

   if(nameOfTree.Contains(nameOfAnaJPsiTree.Data()) ){ // if this is AnaJPsi, create one more plot with cut out k, p
      condition = getCondition1Track(0, "nsigmaproton nsigmakaon");
      cout << "Condition for JPsi with k,p cut out: " << condition << endl;
   
      tree->Draw(whatToDraw, condition); //momentumChargePlus or Minus

      condition = getCondition1Track(1, "nsigmaproton nsigmakaon");
      
      tree->Draw(whatToDraw2, condition);

      TH2D* hdEdx2 = (TH2D*)gPad->GetPrimitive("hDEdx");
      if(!hdEdx2 || hdEdx2->IsZombie() || hdEdx2->GetEntries() == 0){
         cout << "Empty 2D histogram. Leaving..." << endl;
         return;
      }

      if(drawPIDPlot(hdEdx2, "nsigmaCut")){
         cout << "PID plot 4 drawn successfully." << endl;
      } else {
         cout << "Failed to draw PID plot." << endl;
         return;
      }

   }

   if(nameOfTree.Contains(nameOfAnaJPsiTree.Data()) ){ // if this is AnaJPsi, create one more plot with cut out k, p
      condition = "invMass > 2.88 && invMass < 3.22 && " + getCondition1Track(0, "nsigmaproton nsigmakaon");
   
      tree->Draw(whatToDraw, condition); //momentumChargePlus or Minus
      condition = "invMass > 2.88 && invMass < 3.22 && " + getCondition1Track(1, "nsigmaproton nsigmakaon");
      tree->Draw(whatToDraw2, condition);

      TH2D* hdEdx2 = (TH2D*)gPad->GetPrimitive("hDEdx");
      if(!hdEdx2 || hdEdx2->IsZombie() || hdEdx2->GetEntries() == 0){
         cout << "Empty 2D histogram. Leaving..." << endl;
         return;
      }

      if(drawPIDPlot(hdEdx2, "invMassJPsi_only")){
         cout << "PID plot 5 drawn successfully." << endl;
      } else {
         cout << "Failed to draw PID plot." << endl;
         return;
      }

   }


   // just a test of electron line
   if(nameOfTree.Contains(nameOfAnaJPsiTree.Data()) ){ // if this is AnaJPsi, create one more plot with cut out k, p
      condition = getCondition1Track(0, "") + " && abs(nSigmaTPCelectronPlus) > 1.0";
   
      tree->Draw(whatToDraw, condition); //momentumChargePlus or Minus
      condition = getCondition1Track(1, "") + " && abs(nSigmaTPCelectronMinus) > 1.0";
      tree->Draw(whatToDraw2, condition);

      TH2D* hdEdx2 = (TH2D*)gPad->GetPrimitive("hDEdx");
      if(!hdEdx2 || hdEdx2->IsZombie() || hdEdx2->GetEntries() == 0){
         cout << "Empty 2D histogram. Leaving..." << endl;
         return;
      }

      if(drawPIDPlot(hdEdx2, "noElectronsPlot")){
         cout << "PID plot 6 drawn successfully." << endl;
      } else {
         cout << "Failed to draw PID plot." << endl;
         return;
      }

   }

   outFile->Close();

   cout << "Finished all bichselG10 PID plots." << endl;



}