#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include<cmath>
//using namespace RooFit;
void Plot(void)
{

Int_t intb,k,j;

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

 chain->Add("PhiPi0_4s_GMC_0.root");

TCanvas* c = new TCanvas();


//  h1->SetBranchAddress("Mygenfla",&f_Mygenfla);
//  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
//  chain->SetBranchAddress("Phimass",&f_phimass);
//  chain->SetBranchAddress("Dzeromas",&f_d0mass);

//h1->Draw("Pizmassn","Mygenfla ==1 ");
//h1->Draw("Phimass","Mygenfla ==1 ");
h1->Draw("Deltam","Mygenfla ==1 ");

}//end 



