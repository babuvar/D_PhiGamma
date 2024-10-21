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
void EffLoss(void)
{

Int_t intb,k,j;

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

 chain->Add("PhiPi0_4s_GMC_0.root");


  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass,f_dstid,f_dstf,f_hel,f_dstm,f_d0f,f_pisf,f_phif,f_pi0f, f_kpf,f_kmf,f_multiplicity,f_pmultiplicity,f_nmultiplicity,f_d0id, NS[70][70]={0},NB[70][70]={0},sig[70][70],bestsig=1000000,besti=3, bestj=20, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e, f_Mygenfla;

  h1->SetBranchAddress("Mygenfla",&f_Mygenfla);
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

int totsig=0;
int survsig1=0; float effloss1;
int survsig2=0; float effloss2;
int survsig3=0; float effloss3;
int survsig4=0; float effloss4;
int survsig5=0; float effloss5;
int survsig6=0; float effloss6;
int survsig7=0; float effloss7;
int survsig8=0; float effloss8;
int survsig9=0; float effloss9;









  Int_t nevt=(int)chain->GetEntries();

  for(Int_t i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
/*
//Phot1cut
if(f_Phot1the > -60 && f_Phot1the < 67 && f_Photon1e > 0.05){photon1cutflag=1;}
else if((f_Phot1the < -60 || f_Phot1the > 67) && f_Photon1e > 0.1){photon1cutflag=1;}
else{photon1cutflag=0;}

//Phot2cut
if(f_Phot2the > -60 && f_Phot2the < 67 && f_Photon2e > 0.05){photon2cutflag=1;}
else if((f_Phot2the < -60 || f_Phot2the > 67) && f_Photon2e > 0.1){photon2cutflag=1;}
else{photon2cutflag=0;}
*/
//if(f_dstf == 1 && fabs(f_dstid)==413){
if(f_Mygenfla==1){
totsig++;

//if(f_delm > 0.14 && f_delm < 0.16) {survsig1++;}
if(f_phimass > 1.01 && f_phimass < 1.03) {survsig2++;}
//if(f_pi0mass > 0.119 && f_pi0mass < 0.151) {survsig3++;}
//if(f_pi0mom > 0.38 && f_dstps > 2.5) {survsig4++;}
//if(f_d0mass > 1.83 && f_d0mass < 1.89) {survsig5++;}
//if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){survsig6++;}
//if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){survsig7++;}
//if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){survsig8++;}
//if(fabs(f_Pispid) < 0.9){survsig9++;}
//if(photon1cutflag == 1 && photon1cutflag == 1){

}


    }

//cout<<"totsig = "<<totsig<<endl;
//cout<<"survsig = "<<survsig<<endl;
effloss1=(float)(totsig-survsig1)*100/totsig;
effloss2=(float)(totsig-survsig2)*100/totsig;
effloss3=(float)(totsig-survsig3)*100/totsig;
effloss4=(float)(totsig-survsig4)*100/totsig;
effloss5=(float)(totsig-survsig5)*100/totsig;
effloss6=(float)(totsig-survsig6)*100/totsig;
effloss7=(float)(totsig-survsig7)*100/totsig;
effloss8=(float)(totsig-survsig8)*100/totsig;
effloss9=(float)(totsig-survsig9)*100/totsig;


cout<<"Efficiency loss DeltaM = "<<effloss1<<" %"<<endl;
cout<<"Efficiency loss MPhi = "<<effloss2<<" %"<<endl;
cout<<"Efficiency loss Mpi0 = "<<effloss3<<" %"<<endl;
cout<<"Efficiency loss 2-D = "<<effloss4<<" %"<<endl;
cout<<"Efficiency loss M_D = "<<effloss5<<" %"<<endl;
cout<<"Efficiency loss |dr| = "<<effloss6<<" %"<<endl;
cout<<"Efficiency loss |dz| = "<<effloss7<<" %"<<endl;
cout<<"Efficiency loss K-PID = "<<effloss8<<" %"<<endl;
cout<<"Efficiency loss Pi-PID = "<<effloss9<<" %"<<endl;

cout<<totsig<<"\t"<<survsig2<<endl;

}//end 



