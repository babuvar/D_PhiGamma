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
#include "RooGenericPdf.h"
#include "RooAbsPdf.h"
 #include "RooMCStudy.h"
 #include "TH2.h"
 #include "RooFitResult.h"
 #include "TStyle.h"
 #include "TDirectory.h"
using namespace RooFit;
  RooRealVar delmass("delmass","delmass",0.14,0.16,"GeV/c^{2}");
  RooRealVar helicity("helicity","helicity",-1,1,"no unit");
  RooDataSet* data_p1=new RooDataSet("data_p1","data_p1",RooArgSet(helicity,delmass));
  RooDataSet* data_n1=new RooDataSet("data_n1","data_n1",RooArgSet(helicity,delmass));
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(helicity,delmass));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(helicity,delmass));
  RooDataSet* data_p3=new RooDataSet("data_p3","data_p3",RooArgSet(helicity,delmass));
  RooDataSet* data_n3=new RooDataSet("data_n3","data_n3",RooArgSet(helicity,delmass));
  RooDataSet* data_p4=new RooDataSet("data_p4","data_p4",RooArgSet(helicity,delmass));
  RooDataSet* data_n4=new RooDataSet("data_n4","data_n4",RooArgSet(helicity,delmass));
  RooDataSet* data_p5=new RooDataSet("data_p5","data_p5",RooArgSet(helicity,delmass));
  RooDataSet* data_n5=new RooDataSet("data_n5","data_n5",RooArgSet(helicity,delmass));
  RooDataSet* data_p6=new RooDataSet("data_p6","data_p6",RooArgSet(helicity,delmass));
  RooDataSet* data_n6=new RooDataSet("data_n6","data_n6",RooArgSet(helicity,delmass));
  RooDataSet* data_p7=new RooDataSet("data_p7","data_p7",RooArgSet(helicity,delmass));
  RooDataSet* data_n7=new RooDataSet("data_n7","data_n7",RooArgSet(helicity,delmass));
  RooDataSet* data_p8=new RooDataSet("data_p8","data_p8",RooArgSet(helicity,delmass));
  RooDataSet* data_n8=new RooDataSet("data_n8","data_n8",RooArgSet(helicity,delmass));
  RooDataSet* data_p9=new RooDataSet("data_p9","data_p9",RooArgSet(helicity,delmass));
  RooDataSet* data_n9=new RooDataSet("data_n9","data_n9",RooArgSet(helicity,delmass));
  RooDataSet* data_p10=new RooDataSet("data_p10","data_p10",RooArgSet(helicity,delmass));
  RooDataSet* data_n10=new RooDataSet("data_n10","data_n10",RooArgSet(helicity,delmass));
Float_t A[11],B[11];
void fitit(int index);
void getasy();
void simfit_2D(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiPi0_4s_GMC_2.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass,f_hel, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;

  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
  chain->SetBranchAddress("Kplushel",&f_hel);
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
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
		  delmass.setVal(f_delm);
		  helicity.setVal(f_hel);
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

	      if(f_dstcharge==+1)
		{
		  nevt_ac_p++;
//data_p1->add(RooArgSet(helicity,delmass));

if(f_dstcts >= -1 && f_dstcts < -0.8) {data_p1->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.8 && f_dstcts < -0.6) {data_p2->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.6 && f_dstcts < -0.4) {data_p3->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.4 && f_dstcts < -0.2) {data_p4->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.2 && f_dstcts < 0) {data_p5->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0 && f_dstcts < 0.2) {data_p6->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.2 && f_dstcts < 0.4) {data_p7->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.4 && f_dstcts < 0.6) {data_p8->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.6 && f_dstcts < 0.8) {data_p9->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.8 && f_dstcts < 1.0) {data_p10->add(RooArgSet(helicity,delmass));}
		}
	      else
		{
		  nevt_ac_n++;
//data_n1->add(RooArgSet(helicity,delmass));

if(f_dstcts >= -1 && f_dstcts < -0.8) {data_n1->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.8 && f_dstcts < -0.6) {data_n2->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.6 && f_dstcts < -0.4) {data_n3->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.4 && f_dstcts < -0.2) {data_n4->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= -0.2 && f_dstcts < 0) {data_n5->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0 && f_dstcts < 0.2) {data_n6->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.2 && f_dstcts < 0.4) {data_n7->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.4 && f_dstcts < 0.6) {data_n8->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.6 && f_dstcts < 0.8) {data_n9->add(RooArgSet(helicity,delmass));}
if(f_dstcts >= 0.8 && f_dstcts < 1.0) {data_n10->add(RooArgSet(helicity,delmass));}
		}
}}}}}}}}}}}
    }

for(int k=1;k<=10;k++){fitit(k);}
getasy();

    exit();

}//end of simfit

void fitit(int index)
{
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(helicity,delmass));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(helicity,delmass));
  char name[10] = "s01f.png";
  Float_t fraction,sigma1,sigma2,sigma;
switch (index)
{
case 1: data_p=data_p1; data_n=data_n1; name[2] = '1';  fraction= 0.518373; sigma1=0.000951887 ; sigma2=0.00115961 ; sigma= 0.000443002;
    break;
case 2: data_p=data_p2; data_n=data_n2; name[2] = '2';  fraction= 0.620506; sigma1=0.000900976 ; sigma2=0.00111373 ; sigma= 0.000405406;
    break;
case 3: data_p=data_p3; data_n=data_n3;  name[2] = '3';   fraction= 0.650871; sigma1=0.000977793 ; sigma2=0.00113021 ; sigma= 0.000395017;
    break;
case 4: data_p=data_p4; data_n=data_n4;  name[2] = '4';   fraction= 0.634276; sigma1=0.000885001 ; sigma2=0.00105302 ; sigma= 0.000379537;
    break;
case 5: data_p=data_p5; data_n=data_n5;  name[2] = '5';   fraction= 0.609559; sigma1=0.00084263 ; sigma2=0.0010166 ; sigma= 0.000377725;
    break;
case 6: data_p=data_p6; data_n=data_n6;  name[2] = '6';   fraction= 0.664211; sigma1=0.00103672 ; sigma2= 0.00117925; sigma= 0.000408053;
    break;
case 7: data_p=data_p7; data_n=data_n7;  name[2] = '7';   fraction= 0.613573; sigma1=0.00089095 ; sigma2= 0.0010811; sigma= 0.00039559;
    break;
case 8: data_p=data_p8; data_n=data_n8;  name[2] = '8';   fraction= 0.614625; sigma1=0.000934598 ; sigma2= 0.00114655; sigma= 0.000419327;
    break;
case 9: data_p=data_p9; data_n=data_n9;  name[2] = '9';   fraction= 0.534314; sigma1=0.00094043 ; sigma2= 0.00116195; sigma= 0.000433835;
    break;
case 10: data_p=data_p10; data_n=data_n10; name[1] = '1'; name[2] = '0';  fraction= 0.664235; sigma1=0.00122308 ; sigma2=0.00174307 ; sigma=0.000667947 ;
    break;
}

Int_t maxp=data_p->numEntries(),maxn=data_n->numEntries();
//fraction=0.601244;
//sigma1=0.000897836; 
//sigma2=0.00109111;

  //DEFINE PDF
  //____________________________________________________
  //CHARGE=+1

  //Gauss
  RooRealVar mean("mean","mean",0.145428,0.14,0.15);
  RooRealVar sig("sig","sig",0.000397327,0,0.0015);//sigma);//0.000397327,0,0.003);
  RooGaussian sig_pdf_p("sig_pdf_p", "signal Gaussian",delmass,mean,sig);
 
  //Double-Gauss
  RooRealVar sig1("sig1","sig1",sigma1);//0.00089779,0,0.0020);
  RooRealVar sig2("sig2","sig2",sigma2);//0.00109104,0,0.0020);
  RooBifurGauss sig_pdf1_p("sig_pdf1_p", "signal Gaussian",delmass,mean,sig1,sig2);//***************8
  
  //Background
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("a1", "coef 1",0.002,0,0.005);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("a3", "coef 3",-0.119,-20.0,20);
  RooDstD0BG bkg_p("bkg_p","background", delmass, dm0,a1,a2,a3);

  //Full Signal
  RooRealVar f_Sig("f_Sig","f_Sig",fraction);//0.6012,0.5,1);
  RooAddPdf sig_tot_p("sig_tot_p","sig_tot_p",RooArgList(sig_pdf_p,sig_pdf1_p),f_Sig);
  
  //____________________________________________________
  //Araw Defn
  RooRealVar Araw("Araw","Araw",0,-1,1);
  RooRealVar N_t("N_t","N_t",100,0,(maxn+maxp));
  RooFormulaVar N_n("N_n","(0.5)*(1-Araw)*N_t",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+Araw)*N_t",RooArgList(Araw,N_t));
  //____________________________________________________
  //Abkg Defn
  RooRealVar Abkg("Abkg","Abkg",0,-1,1);
  RooRealVar N_tb("N_tb","N_tb",100,0,(maxn+maxp));
  RooFormulaVar N_nb("N_nb","(0.5)*(1-Abkg)*N_tb",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+Abkg)*N_tb",RooArgList(Abkg,N_tb));
  //_____________________________________________________

  RooRealVar shift1("shift1","shift1",0,-0.5,0.5);
  RooRealVar p("p","p",0);  
  RooChebychev bkgH1_p("bkgH1_p", "bkgH1_p",helicity, RooArgList(p));
  RooGenericPdf bkgH2_p("bkgH2_p","(helicity-shift1)*(helicity-shift1)",RooArgSet(helicity,shift1));
  RooRealVar f_Sig2("f_Sig2","f_Sig2",0.6012,0,1);
  RooAddPdf bkgH_p("bkgH_p","bkgH_p",RooArgList(bkgH1_p,bkgH2_p),f_Sig2);
  RooGenericPdf sigH_p("sigH_p","(helicity-shift1)*(helicity-shift1)",RooArgSet(helicity,shift1));

//Full + signal model
RooProdPdf Sig_P("Sig_P","Sig_P",RooArgSet(sig_tot_p,sigH_p));
RooProdPdf Bkg_P("Bkg_P","Bkg_P",RooArgSet(bkg_p,bkgH_p));
  RooAddPdf Model_P("Model_P", "Model_P",RooArgList(Sig_P,Bkg_P), RooArgList(N_p,N_pb));

  //_____________________________________________________
  //CHARGE=-1

  RooGaussian sig_pdf_n("sig_pdf_n", "signal Gaussian",delmass,mean,sig);
  RooBifurGauss sig_pdf1_n("sig_pdf1_n", "signal Gaussian",delmass,mean,sig1,sig2);
  RooDstD0BG bkg_n("bkg_n","background", delmass, dm0,a1,a2,a3);
  RooAddPdf sig_tot_n("sig_tot_n","sig_tot_n",RooArgList(sig_pdf_n,sig_pdf1_n),f_Sig);

  RooRealVar p1("p1","p1",0);
  RooChebychev bkgH1_n("bkgH1_n", "bkgH1_n",helicity, RooArgList(p1));
  RooGenericPdf bkgH2_n("bkgH2_n","(helicity+shift1)*(helicity+shift1)",RooArgSet(helicity,shift1));
  RooAddPdf bkgH_n("bkgH_n","bkgH_n",RooArgList(bkgH1_n,bkgH2_n),f_Sig2);
  RooGenericPdf sigH_n("sigH_n","(helicity+shift1)*(helicity+shift1)",RooArgSet(helicity,shift1));


RooProdPdf Sig_N("Sig_N","Sig_N",RooArgSet(sig_tot_n,sigH_n));
RooProdPdf Bkg_N("Bkg_N","Bkg_N",RooArgSet(bkg_n,bkgH_n));
  RooAddPdf Model_N("Model_N", "Model_N",RooArgList(Sig_N,Bkg_N), RooArgList(N_n,N_nb));



  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  RooDataSet combData("combData","combData",RooArgSet(helicity,delmass),Index(sample),Import("pos",*data_p),Import("neg",*data_n));
  
  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_P,"pos");
  simPdf.addPdf(Model_N,"neg");
  //_____________________________________________________

   RooFitResult* fitRes = simPdf.fitTo(combData,Save());
//   cout<<"var fit status = "<<fitRes->status()<<endl;

  //PLOTING
  RooPlot *xframe_1 =delmass.frame(Bins(100),Title("#pi^{+}_{s}"));
  combData.plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("bkg_p"),ProjWData(sample,combData),LineStyle(kDashed));
   Model_P.paramOn(xframe_1,data_p);
  RooPlot *xframe_2 = delmass.frame(Bins(100),Title("#pi^{-}_{s}"));
  combData.plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("bkg_n"),ProjWData(sample,combData),LineStyle(kDashed));
     Model_N.paramOn(xframe_2,data_n);  
  TCanvas* c = new TCanvas("c","c",1200,800) ;

  c->Divide(2,2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;


  //PLOTING
  RooPlot *yframe_1 =helicity.frame(Bins(100),Title("#pi^{+}_{s}"));
  combData.plotOn(yframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),ProjWData(sample,combData));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("bkgH_p"),ProjWData(sample,combData),LineStyle(kDashed));
//   Model_P.paramOn(yframe_1,data_p);

  RooPlot *yframe_2 = helicity.frame(Bins(100),Title("#pi^{-}_{s}"));
  combData.plotOn(yframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),ProjWData(sample,combData));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("bkgH_n"),ProjWData(sample,combData),LineStyle(kDashed));
//     Model_N.paramOn(yframe_2,data_n);  
//  TCanvas* c2 = new TCanvas("c2","c2",1200,800) ;

//  c2->Divide(2) ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; yframe_1->GetYaxis()->SetTitleOffset(1.4) ; yframe_1->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; yframe_2->GetYaxis()->SetTitleOffset(1.4) ; yframe_2->Draw() ;

  TCanvas* c2 = new TCanvas("c2","c2",1200,800) ;
  c2->Divide(2) ;
  c2->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  c2->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

  TCanvas* c3 = new TCanvas("c3","c3",1200,800) ;
  c3->Divide(2) ;
  c3->cd(1) ; gPad->SetLeftMargin(0.15) ; yframe_1->GetYaxis()->SetTitleOffset(1.4) ; yframe_1->Draw() ;
  c3->cd(2) ; gPad->SetLeftMargin(0.15) ; yframe_2->GetYaxis()->SetTitleOffset(1.4) ; yframe_2->Draw() ;



  c->SaveAs(name);
//  c3->SaveAs("FullFit.png");
  delete c;
  delete c2;
  delete c3;

A[index]=Araw.getVal();
B[index]=Araw.getError();
ofstream fout;
fout.open("Asym.txt",ios::app);
fout<<Araw.getVal()<<"\t"<<Araw.getError()<<"\t"<<(mean.getVal()-0.14)/0.01<<"\t"<<sig.getVal()/0.003<<"\t"<<sig1.getVal()/0.0019<<"\t"<<sig2.getVal()/0.0019<<"\t"<<f_Sig.getVal()/0.75<<"\t"<<fitRes->status()<<endl;
fout.close();


}

void getasy()
{
Int_t i=0;
Float_t M=0,S=0,D1=0;
for(i=1;i<=5;i++){
A[i]=(A[i]+A[11-i])/2;
B[i]=(1/(B[i]*B[i]))+(1/(B[11-i]*B[11-i]));
B[i]=sqrt(B[i]);
B[i]=1/B[i];
//B[i]=(B[i]*B[i])+(B[11-i]*B[11-i]);
//B[i]=sqrt(B[i]);
D1=D1+(1/(B[i]*B[i]));
M=M+(A[i]/(B[i]*B[i]));
}
S=sqrt(D1);
S=1/S;
M=M/D1;
ofstream fout;
fout.open("Asym.txt",ios::app);
//for(i=1;i<=10;i++){fout<<A[i]<<"\t"<<B[i]<<endl;}
fout<<endl<<"Combined ACP = "<<M<<" Combined Error on ACP = "<<S<<endl;
fout.close();

}





