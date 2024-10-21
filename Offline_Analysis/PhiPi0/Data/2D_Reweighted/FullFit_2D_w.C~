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
  RooRealVar w("w","w",-15,15) ;
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(helicity,delmass,w));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(helicity,delmass,w));

float x1,x2,x3,x4,x5,w1,w2,w3,w4,w5,dpw[50],dnw[50];
x1=0.0233808; x2=0.0283486; x3=0.0208331; x4=0.0174581; x5=0.0108812;
w1=(1+x1)/(1-x1); w2=(1+x2)/(1-x2); w3=(1+x3)/(1-x3); w4=(1+x4)/(1-x4); w5=(1+x5)/(1-x5);
dpw[0]=1.0;dnw[0]=w1;
dpw[1]=1.0;dnw[1]=w2;
dpw[2]=1.0;dnw[2]=w3;
dpw[3]=1.0;dnw[3]=w4;
dpw[4]=1.0;dnw[4]=w5;
dpw[5]=w5;dnw[5]=1.0;
dpw[6]=w4;dnw[6]=1.0;
dpw[7]=w3;dnw[7]=1.0;
dpw[8]=w2;dnw[8]=1.0;
dpw[9]=w1;dnw[9]=1.0;

int bin;
float n=5,m;
m=1/n;

void FullFit_2D_w(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiPi0_4s_Data.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass,f_hel,f_piscost,f_dr,f_dz,f_kid, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;

  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
  chain->SetBranchAddress("Kplushel",&f_hel);
  chain->SetBranchAddress("Piscost",&f_piscost);
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
bin=floor(f_dstcts/m)+n;

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
w.setVal(dpw[bin]);
data_p->add(RooArgSet(helicity,delmass,w));

		}
	      else
		{
		  nevt_ac_n++;
w.setVal(dnw[bin]);
data_n->add(RooArgSet(helicity,delmass,w));


		}
}}}}}}}}}}}
    }

Int_t maxp=data_p->numEntries(),maxn=data_n->numEntries();

  //DEFINE PDF
  //____________________________________________________
  //CHARGE=+1

  //Gauss
  RooRealVar mean("mean","mean",0.145428,0.14,0.15);
  RooRealVar sig("sig","sig",0.000397327,0,0.0015);//sigma);//0.000397327,0,0.003);
  RooGaussian sig_pdf_p("sig_pdf_p", "signal Gaussian",delmass,mean,sig);
 
  //Double-Gauss
  RooRealVar sig1("sigL","sigL",0.00089779,0.0008,0.0025);
  RooRealVar sig2("sigR","sigR",0.00109104,0.0008,0.0025);
  RooBifurGauss sig_pdf1_p("sig_pdf1_p", "signal Gaussian",delmass,mean,sig1,sig2);//***************8
  
  //Background
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("a1", "coef 1",0.002,0,0.01);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("a3", "coef 3",-0.119,-20.0,20);
  RooDstD0BG bkg_p("bkg_p","background", delmass, dm0,a1,a2,a3);

  //Full Signal
  RooRealVar f_Sig("f_Sig","f_Sig",0.6012,0.1,1);
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

  RooRealVar shift1("shift1","shift1",0);//,-0.5,0.5);
  RooRealVar shift2("shift2","shift2",0,-0.5,0.5);
  RooRealVar p("slope1","p",0);  
  RooChebychev bkgH1_p("bkgH1_p", "bkgH1_p",helicity, RooArgList(p));
  RooGenericPdf bkgH2_p("bkgH2_p","(helicity-shift1)*(helicity-shift1)",RooArgSet(helicity,shift1));
  RooRealVar f_Sig2("f_Sig2","f_Sig2",0.6012,0,1);
  RooAddPdf bkgH_p("bkgH_p","bkgH_p",RooArgList(bkgH1_p,bkgH2_p),f_Sig2);
  RooGenericPdf sigH_p("sigH_p","(helicity-shift2)*(helicity-shift2)",RooArgSet(helicity,shift2));

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

  RooRealVar p1("slope2","p1",0,-1,1);
  RooChebychev bkgH1_n("bkgH1_n", "bkgH1_n",helicity, RooArgList(p1));
  RooGenericPdf bkgH2_n("bkgH2_n","(helicity+shift1)*(helicity+shift1)",RooArgSet(helicity,shift1));
  RooAddPdf bkgH_n("bkgH_n","bkgH_n",RooArgList(bkgH1_n,bkgH2_n),f_Sig2);

  RooGenericPdf sigH_n("sigH_n","(helicity+shift2)*(helicity+shift2)",RooArgSet(helicity,shift2));

RooProdPdf Sig_N("Sig_N","Sig_N",RooArgSet(sig_tot_n,sigH_n));
RooProdPdf Bkg_N("Bkg_N","Bkg_N",RooArgSet(bkg_n,bkgH_n));
  RooAddPdf Model_N("Model_N", "Model_N",RooArgList(Sig_N,Bkg_N), RooArgList(N_n,N_nb));



  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(helicity,delmass,w),Index(sample),Import("pos",*data_p),Import("neg",*data_n));
  

  RooDataSet* combData_w=new RooDataSet(combData->GetName(),combData->GetTitle(),combData,*combData->get(),0,w.GetName()) ;




  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_P,"pos");
  simPdf.addPdf(Model_N,"neg");
  //_____________________________________________________

//   RooFitResult* fitRes = simPdf.fitTo(combData,Save());
   RooFitResult* fitRes = simPdf.fitTo(*combData_w,Save(),SumW2Error(kTRUE));

//   cout<<"var fit status = "<<fitRes->status()<<endl;

  //PLOTING
  RooPlot *xframe_1 =delmass.frame(Bins(100),Title("D^{0} #rightarrow #phi #pi^{0}"));
  combData_w.plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("bkg_p"),ProjWData(sample,*combData_w),LineStyle(kDashed));
   Model_P.paramOn(xframe_1,data_p);
  RooPlot *xframe_2 = delmass.frame(Bins(100),Title("D^{0} #rightarrow #phi #pi^{0}"));
  combData_w.plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("bkg_n"),ProjWData(sample,*combData_w),LineStyle(kDashed));
     Model_N.paramOn(xframe_2,data_n);  
  TCanvas* c = new TCanvas("c","c",1200,800) ;

  c->Divide(2,2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;


  //PLOTING
  RooPlot *yframe_1 =helicity.frame(Bins(100),Title("D^{0} #rightarrow #phi #pi^{0}"));
  combData_w.plotOn(yframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),ProjWData(sample,*combData_w));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("bkgH_p"),ProjWData(sample,*combData_w),LineStyle(kDashed));
//   Model_P.paramOn(yframe_1,data_p);

  RooPlot *yframe_2 = helicity.frame(Bins(100),Title(" D^{0} #rightarrow #phi #pi^{0}"));
  combData_w.plotOn(yframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),ProjWData(sample,*combData_w));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("bkgH_n"),ProjWData(sample,*combData_w),LineStyle(kDashed));
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




cout<<"Araw =\t"<<Araw.getVal()<<" +/- "<<Araw.getError()<<endl;



}





