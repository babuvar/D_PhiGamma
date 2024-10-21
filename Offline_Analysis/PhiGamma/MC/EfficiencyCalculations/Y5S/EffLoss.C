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
void EffLoss(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  chain->Add("DtoPhiGamma_80760x2.root");
  
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
int survsig10=0; float effloss10;
int survsig11=11; float effloss11;

  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);

if(Mygenfla==1){

totsig++;

if(Deltam > 0.14 && Deltam < 0.16){survsig1++;}
if(MD0 > 1.68 && MD0 < 2.05){survsig2++;}
if(MPhi > 1.01 && MPhi < 1.03){survsig3++;}
if(ProbPi < 0.055){survsig4++;}
if(ProbEta < 0.39){survsig5++;}
if(E9E25 > 0.941){survsig6++;}
if(Gamener > 0.61 && Dstp > 3.1){survsig7++;}
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){survsig8++;}
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){survsig9++;}
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){survsig10++;}
if(fabs(f_Pispid) < 0.9){survsig11++;}



}

}

effloss1=(float)(totsig-survsig1)*100/totsig;
effloss2=(float)(totsig-survsig2)*100/totsig;
effloss3=(float)(totsig-survsig3)*100/totsig;
effloss4=(float)(totsig-survsig4)*100/totsig;
effloss5=(float)(totsig-survsig5)*100/totsig;
effloss6=(float)(totsig-survsig6)*100/totsig;
effloss7=(float)(totsig-survsig7)*100/totsig;
effloss8=(float)(totsig-survsig8)*100/totsig;
effloss9=(float)(totsig-survsig9)*100/totsig;
effloss10=(float)(totsig-survsig10)*100/totsig;
effloss11=(float)(totsig-survsig11)*100/totsig;



cout<<"Efficiency loss DeltaM = "<<effloss1<<" %"<<endl;
cout<<"Efficiency loss M_D = "<<effloss2<<" %"<<endl;
cout<<"Efficiency loss Mphi = "<<effloss3<<" %"<<endl;
cout<<"Efficiency loss Prob(Pi) = "<<effloss4<<" %"<<endl;
cout<<"Efficiency loss Prob(Eta) = "<<effloss5<<" %"<<endl;
cout<<"Efficiency loss e9/e25 = "<<effloss6<<" %"<<endl;
cout<<"Efficiency loss 2-D = "<<effloss7<<" %"<<endl;
cout<<"Efficiency loss |dr| = "<<effloss8<<" %"<<endl;
cout<<"Efficiency loss |dz| = "<<effloss9<<" %"<<endl;
cout<<"Efficiency loss K-PID = "<<effloss10<<" %"<<endl;
cout<<"Efficiency loss Pi-PID = "<<effloss11<<" %"<<endl;






}

