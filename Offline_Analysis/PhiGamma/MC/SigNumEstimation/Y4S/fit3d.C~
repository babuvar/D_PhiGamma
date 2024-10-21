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
//  RooRealVar delmass("#Delta M","#Delta M",0.14,0.16,"GeV/c^{2}");
  RooRealVar delmass("delmass","delmass",0.14,0.16,"GeV/c^{2}");
//  RooRealVar helicity("cos(#theta_{hel})","helicity",-1,1,"no unit");
  RooRealVar helicity("helicity","cos(#theta_{hel})",-1,1,"no unit");
  RooRealVar dzero("M_{D^{0}}","M_{D^{0}}",1.68,2.05,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(helicity,delmass,dzero));
void fit3d(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiGamma_6xGMC.root");
//  chain->Add("PhiGamma_4s_data.root");

  
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

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  Int_t numsig=0;

  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
		  
		  delmass.setVal(Deltam);
		  dzero.setVal(MD0);

if(Deltam > 0.14 && Deltam < 0.16){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.05){
if(ProbEta < 0.215){
if(E9E25 > 0.938){
if(Gamener > 0.58){
if(Dstp > 2.52){
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID

if(Dstarcha == 1) {helicity.setVal(f_helicity);}
else if(Dstarcha == -1) {helicity.setVal(-f_helicity);}


data->add(RooArgSet(helicity,delmass,dzero));
		

}}}}}}}}}}}}
    }




  //DEFINE PDF
  //COMMON 
  RooRealVar p("p","p",0);
  RooRealVar mean("#mu","mean",0.1454327,0.14,0.15);
  RooRealVar sigma("#sigma","sigma",0.0003896,0,0.0015);
  RooRealVar sigb("#sigma_{Bkg}","sigb",0.000372,0,0.0015);
  RooRealVar sig1("sig1","sig1",0.00085);//**
  RooRealVar sig2("sig2","sig2",0.001025);//**
  RooRealVar sig1b("sig1b","sig1b",0.000758);//**
  RooRealVar sig2b("sig2b","sig2b",0.000865);//**
  RooRealVar f_Sig("f_Sig","f_Sig",0.769);//**
  RooRealVar f_Bkg2("f_Bkg2","f_Bkg2",0.562);//**
  RooRealVar F_Bkg("F_{peak/comb}","F_Bkg",0.5,0,0.8);
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("#alpha", "coef 1",0.002,0,0.007);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("#beta", "coef 3",-0.119,-20.0,20);
  RooRealVar m("m","m",1.86481);//**
  RooRealVar s("s","s",0.011571);//**
  RooRealVar a("a","a",0.734);//**
  RooRealVar n("n","n",6.3);//**
  RooRealVar p1("p1","p1",0.0);   
  RooRealVar c1("param_{expo}","c1",-9,-12,0);
  RooRealVar f1("frac_{cheby/expo}","f1",0.5,0,1);
  RooRealVar m1("m1","m1",1.84643);//**
  RooRealVar s1("s1","s1",0.01364);//**
  RooRealVar A1("A1","A1",0.163);//**
  RooRealVar n1("n1","n1",29);//**
  RooRealVar f2("frac_{cheby/1-x^{2}}","f2",0.2,0,1);
  RooRealVar shift1("shift1","shift1",0,-0.5,0.5);
  RooRealVar shift2("shift2","shift2",0,-0.5,0.5);
  //____________________________________________________

  RooRealVar N("N_{Sig}","N",100,0,1000000);
  RooRealVar N_b("N_{Bkg}","N_b",100,0,1000000);


//SIGNAL

//Signal DeltaM
  //Gauss
  RooGaussian sig_pdf("sig_pdf", "signal Gaussian",delmass,mean,sigma);
 
  //Double-Gauss
  RooBifurGauss sig_pdf1("sig_pdf1", "signal Gaussian",delmass,mean,sig1,sig2);//***************8

  //Full Signal
  RooAddPdf sig_tot("sig_tot","sig_tot",RooArgList(sig_pdf,sig_pdf1),f_Sig);

//Signal Helicity
  RooGenericPdf sig("sig","1-(helicity*helicity)",RooArgSet(helicity));

//Signal D0
  RooCBShape sig_dz("sig_dz","sig_dz", dzero, m, s, a, n); 

//Product Signal
RooProdPdf Sig("Sig","Sig",RooArgSet(sig_tot,sig,sig_dz));

//BACKGROUND

//Background1 DeltaM
  RooDstD0BG bkg("bkg","background", delmass, dm0,a1,a2,a3);

//Background1 Helicity
  RooChebychev bkgh1("bkgh1", "bkgh1",helicity, RooArgList(p));
  RooGenericPdf bkgh2("bkgh2","(helicity-shift1)*(helicity-shift1)",RooArgSet(helicity,shift1));
  RooAddPdf bkgH1("bkgH1","bkgH1",RooArgList(bkgh1,bkgh2),f2);

//Background1 D0 
  RooChebychev bkg1_D0("bkg1_D0", "bkg1_D0", dzero, RooArgList(p1));
  RooExponential bkg2_D0("bkg2_D0", "bkg2_D0", dzero, c1); 
  RooAddPdf Bkg_D0("Bkg_D0","Bkg_D0",RooArgList(bkg1_D0,bkg2_D0),RooArgList(f1));

//Product Background1
RooProdPdf Bkg1("Bkg1","Bkg1",RooArgSet(bkg,bkgH1,Bkg_D0));

//Background2 DeltaM
  //Gauss
  RooGaussian bkg_pdf("bkg_pdf", "signal Gaussian",delmass,mean,sigb);
 
  //Double-Gauss
  RooBifurGauss bkg_pdf1("bkg_pdf1", "signal Gaussian",delmass,mean,sig1b,sig2b);//***************8

  //Full Signal
  RooAddPdf bkg_tot2("bkg_tot2","bkg_tot2",RooArgList(bkg_pdf,bkg_pdf1),f_Bkg2);


//Background2 Helicity
  RooGenericPdf bkg2H("bkg2H","(helicity-shift2)*(helicity-shift2)",RooArgSet(helicity,shift2));

//Background2 D0
  RooCBShape Bkg_D0_2("Bkg_D0_2","Bkg_D0_2", dzero, m1, s1, A1, n1); 

//Product Background2
RooProdPdf Bkg2("Bkg2","Bkg2",RooArgSet(bkg_tot2,bkg2H,Bkg_D0_2));

//Full Background P
  RooAddPdf bkg_tot("bkg_tot","bkg_tot",RooArgList(Bkg1,Bkg2),F_Bkg);

//Full Model P
  RooAddPdf Model("Model","Model",RooArgList(Sig,bkg_tot),RooArgList(N,N_b));



   RooFitResult* fitRes = Model.fitTo(*data);



//Defining signal windows
 delmass.setRange("SignalforMD",0.142,0.149);	//DeltaM window 'signal2'
 dzero.setRange("SignalforDelM",1.83,1.89);	//DeltaM window 'signal2'
 delmass.setRange("RectSignal",0.142,0.149);	//DeltaM window 'signal2'
 dzero.setRange("RectSignal",1.83,1.89);	//DeltaM window 'signal2'


  //DeltaM PLOTING
RooPlot *xframe =delmass.frame(Bins(50),Title("D^{0} #rightarrow #phi #gamma"));
data.plotOn(xframe);
Model.plotOn(xframe);
Model.plotOn(xframe,Components("sig_tot"),LineStyle(kDashed));
Model.plotOn(xframe,Components("bkg"),LineStyle(kDashed),LineColor(kGreen));
Model.plotOn(xframe,Components("bkg_tot2"),LineStyle(kDashed),LineColor(kRed));
Model.paramOn(xframe,data); 

  TCanvas* C = new TCanvas("C","C",1200,800) ;
  C->cd() ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;

  //Helicity PLOTING
RooPlot *yframe =helicity.frame(Bins(50),Title("D^{0} #rightarrow #phi #gamma"));
data.plotOn(yframe);
Model.plotOn(yframe);
Model.plotOn(yframe,Components("sig_tot"),LineStyle(kDashed));
Model.plotOn(yframe,Components("bkg"),LineStyle(kDashed),LineColor(kGreen));
Model.plotOn(yframe,Components("bkg_tot2"),LineStyle(kDashed),LineColor(kRed));
//Model.paramOn(yframe,data); 

  TCanvas* C1 = new TCanvas("C1","C1",1200,800) ;
  C1->cd() ; gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;


  //MD PLOTING
RooPlot *zframe =dzero.frame(Bins(50),Title("D^{0} #rightarrow #phi #gamma"));
data.plotOn(zframe);
Model.plotOn(zframe);
Model.plotOn(zframe,Components("sig_tot"),LineStyle(kDashed));
Model.plotOn(zframe,Components("bkg"),LineStyle(kDashed),LineColor(kGreen));
Model.plotOn(zframe,Components("bkg_tot2"),LineStyle(kDashed),LineColor(kRed));
//Model.paramOn(zframe,data); 

  TCanvas* C2 = new TCanvas("C2","C2",1200,800) ;
  C2->cd() ; gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;


}

