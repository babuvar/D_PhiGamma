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
  RooRealVar delmass("delmass","delmass",0.14,0.16,"GeV/c^{2}");
  RooDataSet* data_1=new RooDataSet("data_1","data_1",RooArgSet(delmass));
  RooDataSet* data_2=new RooDataSet("data_2","data_2",RooArgSet(delmass));
  RooDataSet* data_3=new RooDataSet("data_3","data_3",RooArgSet(delmass));
  RooDataSet* data_4=new RooDataSet("data_4","data_4",RooArgSet(delmass));
  RooDataSet* data_5=new RooDataSet("data_5","data_5",RooArgSet(delmass));
  RooDataSet* data_6=new RooDataSet("data_6","data_6",RooArgSet(delmass));
  RooDataSet* data_7=new RooDataSet("data_7","data_7",RooArgSet(delmass));
  RooDataSet* data_8=new RooDataSet("data_8","data_8",RooArgSet(delmass));
  RooDataSet* data_9=new RooDataSet("data_9","data_9",RooArgSet(delmass));
  RooDataSet* data_10=new RooDataSet("data_10","data_10",RooArgSet(delmass));
void fitit(int index);
void simfit_1D(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiPi0_2M.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;

  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
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
  for(int i=0;i<nevt;i++){
//   for(int i=0;i<10000;i++){ 
//    for(int i=0;i < 287036;i++){
      chain->GetEntry(i);
		  delmass.setVal(f_delm);
//Phot1cut
if(f_Phot1the > -60 && f_Phot1the < 67 && f_Photon1e > 0.05){photon1cutflag=1;}
else if((f_Phot1the < -60 || f_Phot1the > 67) && f_Photon1e > 0.1){photon1cutflag=1;}
else{photon1cutflag=0;}

//Phot2cut
if(f_Phot2the > -60 && f_Phot2the < 67 && f_Photon2e > 0.05){photon2cutflag=1;}
else if((f_Phot2the < -60 || f_Phot2the > 67) && f_Photon2e > 0.1){photon2cutflag=1;}
else{photon2cutflag=0;}

if(f_delm > 0.14 && f_delm < 0.16) {
if(f_phimass > 1.01 && f_phimass < 1.03) {
if(f_pi0mass > 0.119 && f_pi0mass < 0.151) {
if(f_pi0mom > 0.38) {//0.75
if(f_dstps > 2.5) {//2.9
if(f_d0mass > 1.83 && f_d0mass < 1.89) {
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID
if(photon1cutflag == 1 && photon1cutflag == 1){


		  nevt_ac_p++;nevt_ac_n++;
if(f_dstcts >= -1 && f_dstcts < -0.8) {data_1->add(RooArgSet(delmass));}
if(f_dstcts >= -0.8 && f_dstcts < -0.6) {data_2->add(RooArgSet(delmass));}
if(f_dstcts >= -0.6 && f_dstcts < -0.4) {data_3->add(RooArgSet(delmass));}
if(f_dstcts >= -0.4 && f_dstcts < -0.2) {data_4->add(RooArgSet(delmass));}
if(f_dstcts >= -0.2 && f_dstcts < 0) {data_5->add(RooArgSet(delmass));}
if(f_dstcts >= 0 && f_dstcts < 0.2) {data_6->add(RooArgSet(delmass));}
if(f_dstcts >= 0.2 && f_dstcts < 0.4) {data_7->add(RooArgSet(delmass));}
if(f_dstcts >= 0.4 && f_dstcts < 0.6) {data_8->add(RooArgSet(delmass));}
if(f_dstcts >= 0.6 && f_dstcts < 0.8) {data_9->add(RooArgSet(delmass));}
if(f_dstcts >= 0.8 && f_dstcts < 1.0) {data_10->add(RooArgSet(delmass));}

}}}}}}}}}}}
    }

cout<<"events read"<<endl;

for(int k=1;k<=10;k++){fitit(k);}


    exit();

}//end of simfit

void fitit(int index)
{
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(delmass));

  char name[10] = "s01f.png";

  Float_t fraction,sigma1,sigma2,sigma;
switch (index)
{
case 1: data=data_1; name[2] = '1';  break;
case 2: data=data_2; name[2] = '2';  break;
case 3: data=data_3; name[2] = '3';  break;
case 4: data=data_4; name[2] = '4';  break;
case 5: data=data_5; name[2] = '5';  break;
case 6: data=data_6; name[2] = '6';  break;
case 7: data=data_7; name[2] = '7';  break;
case 8: data=data_8; name[2] = '8';  break;
case 9: data=data_9; name[2] = '9';  break;
case 10: data=data_10; name[1] ='1'; name[2] ='0'; break;
}


  //DEFINE PDF

  RooRealVar mean("mean","mean",0.145428,0.14,0.15);
  RooRealVar sig("sig","sig",0.000397327,0,0.003);
  RooGaussian sig_pdf("sig_pdf", "signal Gaussian",delmass,mean,sig);

  RooRealVar sig1("sig1","sig1",0.00089779,0,0.003);
  RooRealVar sig2("sig2","sig2",0.00109104,0,0.003);
  RooBifurGauss sig_pdf1("sig_pdf1", "signal Gaussian",delmass,mean,sig1,sig2);
  
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("a1", "coef 1",0.002,0.0,10);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("a3", "coef 3",-0.119,-10.0,10);
  RooDstD0BG bkg("bkg","background", delmass, dm0,a1,a2,a3);

  RooRealVar f_Sig("f_Sig","f_Sig",0.6012,0,1);
  RooAddPdf sig_tot("sig_tot","sig_tot",RooArgList(sig_pdf,sig_pdf1),f_Sig);
  

  RooRealVar N_s("N_s","N_s",1000,0,1000000);
  RooRealVar N_b("N_b","N_b",1000,0,1000000);

  RooAddPdf model("model", "Model",RooArgList(sig_tot,bkg), RooArgList(N_s,N_b));



   RooFitResult* fitRes = model.fitTo(*data);


cout<<"fitting done"<<endl;


  TCanvas* c = new TCanvas() ;
  c->cd();
  RooPlot* delmass_frame = delmass.frame();
  delmass_frame->SetTitle("M_{D*} - M_{D}");
  data->plotOn(delmass_frame);

  model.plotOn(delmass_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
  model.paramOn(delmass_frame);
  delmass_frame->Draw();



  c->SaveAs(name);

  delete c;

ofstream fout;
fout.open("Parms.txt",ios::app);
fout<<"# "<<index<<endl;
fout<<"sigma ="<<sig.getVal()<<endl;
fout<<"sigma1 ="<<sig1.getVal()<<endl;
fout<<"sigma2 ="<<sig2.getVal()<<endl;
fout<<"f_Sig ="<<f_Sig.getVal()<<endl;
fout<<"----------------------------------"<<endl<<endl;
fout.close();


}






