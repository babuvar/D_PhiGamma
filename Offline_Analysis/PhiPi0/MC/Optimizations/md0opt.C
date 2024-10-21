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
void md0opt(void)
{

Int_t intb,k,j;

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

 chain->Add("PhiPi0_4s_GMC_0.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass,f_dstid,f_dstf,f_hel,f_dstm,f_d0f,f_pisf,f_phif,f_pi0f, f_kpf,f_kmf,f_multiplicity,f_pmultiplicity,f_nmultiplicity,f_d0id, NS[70][70]={0},NB[70][70]={0},sig[70][70],bestsig=1000000,besti=3, bestj=20, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;


  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
  chain->SetBranchAddress("Dstid",&f_dstid);
  chain->SetBranchAddress("D0id",&f_d0id);
  chain->SetBranchAddress("Kplushel",&f_hel);
  chain->SetBranchAddress("Dstarmas",&f_dstm);
  chain->SetBranchAddress("Dstf",&f_dstf);
  chain->SetBranchAddress("D0f",&f_d0f);
  chain->SetBranchAddress("Pisf",&f_pisf);
  chain->SetBranchAddress("Phif",&f_phif);
  chain->SetBranchAddress("Pi0f",&f_pi0f);
  chain->SetBranchAddress("Kpf",&f_kpf);
  chain->SetBranchAddress("Kmf",&f_kmf);
  chain->SetBranchAddress("Multipli",&f_multiplicity);
  chain->SetBranchAddress("Pmultipl",&f_pmultiplicity);
  chain->SetBranchAddress("Nmultipl",&f_nmultiplicity);
  chain->SetBranchAddress("Kppid",&f_Kppid);
  chain->SetBranchAddress("Kmpid",&f_Kmpid);
  chain->SetBranchAddress("Pispid",&f_Pispid);
  chain->SetBranchAddress("Kpdr",&f_Kpdr);
  chain->SetBranchAddress("Kmdr",&f_Kmdr);
  chain->SetBranchAddress("Pisdr",&f_Pisdr);
  chain->SetBranchAddress("Kpdz",&f_Kpdz);
  chain->SetBranchAddress("Kmdz",&f_Kmdz);
  chain->SetBranchAddress("Pisdz",&f_Pisdz);
  chain->SetBranchAddress("Phot1the",&f_Phot1the);
  chain->SetBranchAddress("Phot2the",&f_Phot2the);
  chain->SetBranchAddress("Photon1e",&f_Photon1e);
  chain->SetBranchAddress("Photon2e",&f_Photon2e);

int photon1cutflag=0, photon2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(Int_t i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);

//Phot1cut
if(f_Phot1the > -60 && f_Phot1the < 67 && f_Photon1e > 0.05){photon1cutflag=1;}
else if((f_Phot1the < -60 || f_Phot1the > 67) && f_Photon1e > 0.1){photon1cutflag=1;}
else{photon1cutflag=0;}

//Phot2cut
if(f_Phot2the > -60 && f_Phot2the < 67 && f_Photon2e > 0.05){photon2cutflag=1;}
else if((f_Phot2the < -60 || f_Phot2the > 67) && f_Photon2e > 0.1){photon2cutflag=1;}
else{photon2cutflag=0;}


if(fabs(f_delm-0.1454) < 0.0017) {
if(f_phimass > 1.01 && f_phimass < 1.03) {
if(f_pi0mass > 0.119 && f_pi0mass < 0.151) {
if(f_d0mass > 1.65 && f_d0mass < 2.05) {
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID
if(photon1cutflag == 1 && photon1cutflag == 1){




intb= ((f_d0mass-1.65)*100);

if(f_dstf == 1 && fabs(f_dstid)==413){
for(k=intb+1;k<=41;k++){//upper cut
for(j=intb;j>=1;j--){//lower cut
NS[k][j]++;
}}
}

//cout<<"DONE"<<endl;
else{

for(k=intb+1;k<=41;k++){//upper cut
for(j=intb;j>=1;j--){//lower cut
NB[k][j]++;
}}

}



}}}}}}}}}
    }

cout<<"done with assigning"<<endl;


for(j=1;j<=20;j++){

for(k=22;k<=41;k++){
sig[k][j]=NS[k][j]+NB[k][j];
sig[k][j]=sqrt(sig[k][j]);
sig[k][j]=sig[k][j]/NS[k][j];
if(sig[k][j]<=bestsig){bestsig=sig[k][j];besti=k;bestj=j;}
}

}

cout<<" lower cut = "<<(bestj/100)+1.65<<" GeV  upper cut = "<<(besti/100)+1.65<<" GeV"<<endl;



}//end 

