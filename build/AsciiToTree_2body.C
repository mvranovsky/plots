// STARlight ascii to root tree, 2 BODY FINAL STATES ONLY
// taken from STARlight util ConvertStarlightAsciiToTree.C

#include <iostream>
#include <fstream>
#include <sstream>

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;
using namespace TMath;


void AsciiToTree_2body(char ctag[])
{
  char cvar[120];

  // create output tree
  sprintf(cvar,"%s.root",ctag);
  cout << cvar << endl;
  TFile* outFile = new TFile(cvar, "RECREATE");
  if (!outFile) {
    cerr << "    error: could not create output file '" << outFileName << "'" << endl;
    return;
  }
  TTree* T = new TTree("T","T");

  if(!T){
    cerr << "Could not create TTree. Leaving..." << endl;
    return;
  }
  // tree variables
  Int_t evt;
  Double_t mVM;
  Double_t ptVM;
  Double_t rapVM;
  Double_t etaVM;
  Double_t phiVM;
  Double_t pVM;
  Double_t pxVM;
  Double_t pyVM;
  Double_t pzVM;
  Double_t ed1, ed2;
  Double_t ptd1, ptd2;
  Double_t etad1, etad2;
  Double_t zatemchd1, zatemchd2;
  Double_t phiatemchd1, phiatemchd2;
  Int_t wdgemcd1, wdgemcd2;
  // tree branches
  T->Branch("evt",&evt,"evt/I");
  T->Branch("mVM",&mVM,"mVM/D");
  T->Branch("ptVM",&ptVM,"ptVM/D");
  T->Branch("rapVM",&rapVM,"rapVM/D");
  T->Branch("etaVM",&etaVM,"etaVM/D");
  T->Branch("phiVM",&phiVM,"phiVM/D");
  T->Branch("pVM",&pVM,"pVM/D");
  T->Branch("pxVM",&pxVM,"pxVM/D");
  T->Branch("pyVM",&pyVM,"pyVM/D");
  T->Branch("pzVM",&pzVM,"pzVM/D");
  T->Branch("ed1",&ed1,"ed1/D");
  T->Branch("ed2",&ed2,"ed2/D");
  T->Branch("ptd1",&ptd1,"ptd1/D");
  T->Branch("ptd2",&ptd2,"ptd2/D");
  T->Branch("etad1",&etad1,"etad1/D");
  T->Branch("etad2",&etad2,"etad2/D");
  T->Branch("zatemchd1",&zatemchd1,"zatemchd1/D");
  T->Branch("zatemchd2",&zatemchd2,"zatemchd2/D");
  T->Branch("phiatemchd1",&phiatemchd1,"phiatemchd1/D");
  T->Branch("phiatemchd2",&phiatemchd2,"phiatemchd2/D");
  T->Branch("wdgemcd1",&wdgemcd1,"wdgemcd1/I");
  T->Branch("wdgemcd2",&wdgemcd2,"wdgemcd2/I");

  ifstream inFile;
  sprintf(cvar,"%s.out",ctag); cout << cvar << endl;
  inFile.open(cvar);
  unsigned int countLines = 0;
  while (inFile.good()) {
    string       line;
    stringstream lineStream;
		
    // read EVENT
    string label;
    int    eventNmb, nmbTracks;
    if (!getline(inFile, line))
      break;
    ++countLines;
    lineStream.str(line);
    assert(lineStream >> label >> eventNmb >> nmbTracks);
    if (!(label == "EVENT:"))
      continue;

    //if (nmbTracks != 2) {cout << "WRONG NUMBER OF DAUGHTER TRACKS" << endl; return;}

    if (eventNmb%10000==0) {cout << "event: " << eventNmb << endl;}
    // read vertex
    if (!getline(inFile, line))
      break;
    ++countLines;
    lineStream.str(line);
    assert(lineStream >> label);
    assert(label == "VERTEX:");
			
    TLorentzVector vmDaughter1; TLorentzVector vmDaughter2;
    for (int i = 0; i < nmbTracks; ++i) {
      // read tracks
      int    particleCode;
      double momentum[3];
      if (!getline(inFile, line))
	break;
      ++countLines;
      if (i>1) {continue;}
      lineStream.str(line);
      assert(lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2]);
      assert(label == "TRACK:");
      Double_t daughterMass = IDtoMass(particleCode);
      if (daughterMass < 0) {break;}
      const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
			      + momentum[2] * momentum[2] + daughterMass * daughterMass);
      if (i==0) {vmDaughter1.SetPxPyPzE(momentum[0],momentum[1],momentum[2],E);}
      else {vmDaughter2.SetPxPyPzE(momentum[0],momentum[1],momentum[2],E);}
    }
    TLorentzVector vm = vmDaughter1 + vmDaughter2;

    evt = eventNmb;

    mVM = vm.M();
    ptVM = vm.Perp();
    rapVM = vm.Rapidity();
    etaVM = vm.PseudoRapidity();
    phiVM = vm.Phi();
    pVM = vm.P();
    pxVM = vm.Px();
    pyVM = vm.Py();
    pzVM = vm.Pz();

    ed1 = vmDaughter1.E();
    ed2 = vmDaughter2.E();
    ptd1 = vmDaughter1.Perp();
    ptd2 = vmDaughter2.Perp();
    etad1 = vmDaughter1.PseudoRapidity();
    etad2 = vmDaughter2.PseudoRapidity(); 

    // extrapolate daughters to BEMC
    double pt1 = vmDaughter1.Pt(); double pz1 = vmDaughter1.Pz();
    double pt2 = vmDaughter2.Pt(); double pz2 = vmDaughter2.Pz();
    double dz1 = dzhelix(pt1,pz1); double dz2 = dzhelix(pt2,pz2);
    double dphi1 = dphihelix(pt1); double dphi2 = dphihelix(pt2);

    // take Z=0 for event, add delta-z at BEMC; -999. for Pt too low to hit BEMC
    zatemchd1 = dzhelix(pt1,pz1);
    zatemchd2 = dzhelix(pt2,pz2);

    // add/subtract dphi1,2 for +/- charge curvature, and put phi in range (0,2pi)
    if (dphi1<-900.) {phiatemchd1 = -999.;}
    else {phiatemchd1 = vmDaughter1.Phi() + dphi1;
      while (phiatemchd1<0.) {phiatemchd1 += TwoPi();}
      while (phiatemchd1>TwoPi()) {phiatemchd1 -= TwoPi();}
    }
    if (dphi2<-900.) {phiatemchd2 = -999.;}
    else {phiatemchd2 = vmDaughter2.Phi() - dphi2;
      while (phiatemchd2<0.) {phiatemchd2 += TwoPi();}
      while (phiatemchd2>TwoPi()) {phiatemchd2 -= TwoPi();}
    }

    // BEMC wedge # 1-6
    if (phiatemchd1<-900.) {wdgemcd1 = -999;}
    else {wdgemcd1 = (int)(6.*phiatemchd1/TMath::TwoPi()) + 1;}
    if (phiatemchd2<-900.) {wdgemcd2 = -999;}
    else {wdgemcd2 = (int)(6.*phiatemchd2/TMath::TwoPi()) + 1;}

    T->Fill();
  }

  T->Write();
  if (outFile) {
    outFile->Close();
    delete outFile;
  }

  cout << "Finished converting .txt file to a root tree. Data saved to " << cvar << endl;
}

double IDtoMass(int particleCode){
  double mass;
  if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
  else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
  else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
  else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
  else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
  else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
  else if (particleCode == 14)	{mass = 0.93827231;} // proton
  else {
    cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << endl;
    mass = -1.0;
    //			exit(0); 
  } 

  return mass;
}

double dzhelix(double pth, double pzh) // helix delta-z
{

  const double BSTAR = 0.5; // B-field (tesla)
  const double RCYL = 220.; // BEMC cylinder ~radius (cm)
  //const double RCYL = 213.; // TOF cylinder ~radius (cm)

  double rhelix = 100.*pth/(0.3*BSTAR); // helix radius (cm)

  if (rhelix < RCYL/2.) {return -999.;} // pt too small to hit cylinder

  double bhelix = (pzh/pth)*rhelix; // helix z-slope (cm)

  // helix evolution parameter at cylinder
  double t_cyl = ACos(1. - 0.5*(RCYL/rhelix)*(RCYL/rhelix));

  double dz = bhelix*t_cyl;
  return dz;

}

double dphihelix(double pth) // helix delta-phi
{

  const double BSTAR = 0.5; // B-field (tesla)
  const double RCYL = 220.; // BEMC cylinder ~radius (cm)
  //const double RCYL = 213.; // TOF cylinder ~radius (cm)

  double rhelix = 100.*pth/(0.3*BSTAR); // helix radius (cm)

  if (rhelix < RCYL/2.) {return -999.;} // pt too small to hit cylinder

  double dphi = ASin(RCYL/(2.*rhelix));
  return dphi;

}
