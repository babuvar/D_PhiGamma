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
using namespace RooFit;
void TajCuts(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

//  chain->Add("DtoPhiGamma_3550x2.root");
  chain->Add("DtoPhiGamma_7100.root");
//  chain->Add("test/DtoPhiGamma_7100_corr.root");
//  chain->Add("test/DtoPhiGamma_7100_err.root");
//  chain->Add("Mer_2.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, f_helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, Piflag, Mypgenfl, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_helicity2 ;

  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Mygenfla",&Mygenfla);
  h1->SetBranchAddress("Mypgenfl",&Mypgenfl);
  h1->SetBranchAddress("Dstps",&Dstp);
  h1->SetBranchAddress("Kplushel",&f_helicity);
  h1->SetBranchAddress("Gamenerg",&Gamener);
  h1->SetBranchAddress("P_piz",&P_piz);
  h1->SetBranchAddress("P_eta",&P_eta);
  h1->SetBranchAddress("Vetom1",&Vetom1);
  h1->SetBranchAddress("Vetom2",&Vetom2); 
  h1->SetBranchAddress("P_piz_2",&ProbPi);
  h1->SetBranchAddress("P_eta_2",&ProbEta);
  h1->SetBranchAddress("Phimass",&MPhi);
  h1->SetBranchAddress("Dzeromas",&MD0);
  h1->SetBranchAddress("Ga_e9e25",&E9E25);
  h1->SetBranchAddress("Piflag",&Piflag);
  h1->SetBranchAddress("Dstarcha",&Dstarcha);
  h1->SetBranchAddress("Kppid",&f_Kppid);
  h1->SetBranchAddress("Kmpid",&f_Kmpid);
  h1->SetBranchAddress("Pispid",&f_Pispid);
  h1->SetBranchAddress("Kpdr",&f_Kpdr);
  h1->SetBranchAddress("Kmdr",&f_Kmdr);
  h1->SetBranchAddress("Pisdr",&f_Pisdr);
  h1->SetBranchAddress("Kpdz",&f_Kpdz);
  h1->SetBranchAddress("Kmdz",&f_Kmdz);
  h1->SetBranchAddress("Pisdz",&f_Pisdz);

int totsig=0,totsigp=0,totsign=0;



  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);

if(Deltam > 0.1434 && Deltam < 0.1474){
if(MPhi > 1.01 && MPhi < 1.03){
if(Vetom1 < 0.119 || Vetom1 > 0.151){
if(Gamener > 0.450 && Dstp > 2.9){
if(fabs(f_helicity) < 0.4){

if(fabs(f_Kpdr) < 0.5 && fabs(f_Kmdr) < 0.5 && fabs(f_Pisdr) < 0.5){
if(fabs(f_Kpdz) < 1.5 && fabs(f_Kmdz) < 1.5 && fabs(f_Pisdz) < 1.5){
if(fabs(f_Kppid) > 0.6 && fabs(f_Kmpid) > 0.6){
if(fabs(f_Pispid) < 0.1){


if(Mygenfla==1){totsig++;
if(Dstarcha==1){totsigp++;}
else{totsign++;}
}

}}}}}}}}}

}

cout<<"Total Signal= "<<totsig<<endl;
cout<<"Total Signal P= "<<totsigp<<endl;
cout<<"Total Signal N= "<<totsign<<endl;



}




