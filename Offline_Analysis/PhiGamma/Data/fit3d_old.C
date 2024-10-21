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
  RooRealVar helicity("helicity","helicity",-1,1,"no unit");
  RooRealVar dzero("dzero","dzero",1.68,2.05,"GeV");
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(helicity,delmass,dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(helicity,delmass,dzero));
void fit3d_old(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiGamma_4s_data.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, f_helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, Piflag, Mypgenfl, ;

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

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
		  helicity.setVal(f_helicity);
		  delmass.setVal(Deltam);
		  dzero.setVal(MD0);

if(Deltam > 0.14 && Deltam < 0.16){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.055){
if(ProbEta < 0.289){
if(E9E25 > 0.944){
if(Gamener > 0.6){
if(Dstp > 2.57){

	      if(Dstarcha==+1)
		{
		  nevt_ac_p++;
data_p->add(RooArgSet(helicity,delmass,dzero));
		}
	      else if(Dstarcha==-1)
		{
		  nevt_ac_n++;
data_n->add(RooArgSet(helicity,delmass,dzero));
		}
}}}}}}}}
    }

chain->Add("PhiGamma_5s_data.root");
  Int_t nevt1=(int)chain->GetEntries();
  for(int i=nevt;i<nevt1;i++) 
    {
      chain->GetEntry(i);
		  helicity.setVal(f_helicity);
		  delmass.setVal(Deltam);
		  dzero.setVal(MD0);

if(Deltam > 0.14 && Deltam < 0.16){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.056){
if(ProbEta < 0.285){
if(E9E25 > 0.947){
if(Gamener > 0.59){
if(Dstp > 3.1){

	      if(Dstarcha==+1)
		{
		//  nevt_ac_p++;
data_p->add(RooArgSet(helicity,delmass,dzero));
		}
	      else if(Dstarcha==-1)
		{
		  //nevt_ac_n++;
data_n->add(RooArgSet(helicity,delmass,dzero));
		}
}}}}}}}}
    }






  //DEFINE PDF
  //COMMON 
  RooRealVar p("p","p",0);//,-0.5,0.5);
  RooRealVar mean("mean","mean",0.1454327,0.14,0.15);
  RooRealVar sig("sig","sig",0.0003896,0,0.0015);
  RooRealVar sigb("sigb","sigb",0.000372,0,0.0015);//sigma);//0.000397327,0,0.003);
//  RooRealVar sig1("sig1","sig1",0.000841);//2D  ,0.00089779,0,0.0050);
//  RooRealVar sig2("sig2","sig2",0.001012);//2D  ,0.00109104,0,0.0020);
//  RooRealVar sig1b("sig1b","sig1b",0.000766);//2D  ,0.00089779,0,0.0050);
//  RooRealVar sig2b("sig2b","sig2b",0.000884);//2D  ,0.00109104,0,0.0020);
  RooRealVar sig1("sig1","sig1",0.000963);//3D  ,0.00089779,0,0.0050);
  RooRealVar sig2("sig2","sig2",0.001148);//3D  ,0.00109104,0,0.0020);
  RooRealVar sig1b("sig1b","sig1b",0.001267);//3D  ,0.00089779,0,0.0050);
  RooRealVar sig2b("sig2b","sig2b",0.001540);//3D  ,0.00109104,0,0.0020);
//  RooRealVar f_Sig("f_Sig","f_Sig",0.763);//2D  ,0.6012,0.5,1);
  RooRealVar f_Sig("f_Sig","f_Sig",0.7281);//3D  ,0.6012,0.5,1);
//  RooRealVar f_Bkg2("f_Bkg2","f_Bkg2",0.574);//2D  ,0.6,0.5,1);
  RooRealVar f_Bkg2("f_Bkg2","f_Bkg2",0.423);//3D  ,0.6,0.5,1);
  RooRealVar F_Bkg("F_Bkg","F_Bkg",0.5,0,0.8);
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("a1", "coef 1",0.002,0,0.007);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("a3", "coef 3",-0.119,-20.0,20);
  RooRealVar m("m","m",1.86494);//1.86,1.83,1.88);//
  RooRealVar s("s","s",0.011488);//0.01,0,0.04);//
  RooRealVar a("a","a",0.728);//3,0,10);
  RooRealVar n("n","n",6.42);//2,0,150);
  RooRealVar p1("p1","p1",0.0);   
  RooRealVar c1("c1","c1",-9,-12,0);
  RooRealVar f1("f1","f1",0.5,0,1);
   RooRealVar m1("m1","m1",1.84686);//1.86,1.83,1.88);//
  RooRealVar s1("s1","s1",0.01336);//0.01,0,0.04);//
  RooRealVar A1("A1","A1",0.1605);//3,0,10);
  RooRealVar n1("n1","n1",93.07);//2,0,150);
  RooRealVar f2("f2","f2",0.2,0,1);
  //____________________________________________________
  RooRealVar Araw("Araw","Araw",0,-1,1);
  RooRealVar N_t("N_t","N_t",100,0,1000000);
  RooFormulaVar N_n("N_n","(0.5)*(1-Araw)*N_t",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+Araw)*N_t",RooArgList(Araw,N_t));
  //____________________________________________________
  RooRealVar Abkg("Abkg","Abkg",0,-1,1);
  RooRealVar N_tb("N_tb","N_tb",100,0,1000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-Abkg)*N_tb",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+Abkg)*N_tb",RooArgList(Abkg,N_tb));
  //_____________________________________________________


  //CHARGE=+1

//SIGNAL

//Signal DeltaM
  //Gauss
  RooGaussian sig_pdf_p("sig_pdf_p", "signal Gaussian",delmass,mean,sig);
 
  //Double-Gauss
  RooBifurGauss sig_pdf1_p("sig_pdf1_p", "signal Gaussian",delmass,mean,sig1,sig2);//***************8

  //Full Signal
  RooAddPdf sig_tot_p("sig_tot_p","sig_tot_p",RooArgList(sig_pdf_p,sig_pdf1_p),f_Sig);

//Signal Helicity
  RooGenericPdf sig_p("sig_p","1-(helicity*helicity)",RooArgSet(helicity));

//Signal D0
  RooCBShape sig_dz_p("sig_dz_p","sig_dz_p", dzero, m, s, a, n); 

//Product Signal
RooProdPdf Sig_P("Sig_P","Sig_P",RooArgSet(sig_tot_p,sig_p,sig_dz_p));

//BACKGROUND

//Background1 DeltaM
  RooDstD0BG bkg_p("bkg_p","background", delmass, dm0,a1,a2,a3);

//Background1 Helicity
  RooChebychev bkgh1_p("bkgh1_p", "bkgh1_p",helicity, RooArgList(p));
  RooGenericPdf bkgh2_p("bkgh2_p","helicity*helicity",RooArgSet(helicity));
  RooAddPdf bkgH1_p("bkgH1_p","bkgH1_p",RooArgList(bkgh1_p,bkgh2_p),f2);

//Background1 D0 
  RooChebychev bkg1_D0_p("bkg1_D0_p", "bkg1_D0_p", dzero, RooArgList(p1));
  RooExponential bkg2_D0_p("bkg2_D0_p", "bkg2_D0_p", dzero, c1); 
  RooAddPdf Bkg_D0_p("Bkg_D0_p","Bkg_D0_p",RooArgList(bkg1_D0_p,bkg2_D0_p),RooArgList(f1));

//Product Background1
RooProdPdf Bkg1_P("Bkg1_P","Bkg1_P",RooArgSet(bkg_p,bkgH1_p,Bkg_D0_p));

//Background2 DeltaM
  //Gauss
  RooGaussian bkg_pdf_p("bkg_pdf_p", "signal Gaussian",delmass,mean,sigb);
 
  //Double-Gauss
  RooBifurGauss bkg_pdf1_p("bkg_pdf1_p", "signal Gaussian",delmass,mean,sig1b,sig2b);//***************8

  //Full Signal
  RooAddPdf bkg_tot2_p("bkg_tot2_p","bkg_tot2_p",RooArgList(bkg_pdf_p,bkg_pdf1_p),f_Bkg2);


//Background2 Helicity
  RooGenericPdf bkg2H_p("bkg2H_p","helicity*helicity",RooArgSet(helicity));

//Background2 D0
  RooCBShape Bkg_D0_2p("Bkg_D0_2p","Bkg_D0_2p", dzero, m1, s1, A1, n1); 

//Product Background2
RooProdPdf Bkg2_P("Bkg2_P","Bkg2_P",RooArgSet(bkg_tot2_p,bkg2H_p,Bkg_D0_2p));

//Full Background P
  RooAddPdf bkg_tot_p("bkg_tot_p","bkg_tot_p",RooArgList(Bkg1_P,Bkg2_P),F_Bkg);

//Full Model P
  RooAddPdf Model_P("Model_P","Model_P",RooArgList(Sig_P,bkg_tot_p),RooArgList(N_p,N_pb));


  //CHARGE=-1

//Signal DeltaM
  //Gauss
  RooGaussian sig_pdf_n("sig_pdf_n", "signal Gaussian",delmass,mean,sig);
 
  //Double-Gauss
  RooBifurGauss sig_pdf1_n("sig_pdf1_n", "signal Gaussian",delmass,mean,sig1,sig2);//***************8

  //Full Signal
  RooAddPdf sig_tot_n("sig_tot_n","sig_tot_n",RooArgList(sig_pdf_n,sig_pdf1_n),f_Sig);

//Signal Helicity
  RooGenericPdf sig_n("sig_n","1-(helicity*helicity)",RooArgSet(helicity));

//Signal D0
  RooCBShape sig_dz_n("sig_dz_n","sig_dz_n", dzero, m, s, a, n); 

//Product Signal
RooProdPdf Sig_N("Sig_N","Sig_N",RooArgSet(sig_tot_n,sig_n,sig_dz_n));

//BACKGROUND

//Background1 DeltaM
  RooDstD0BG bkg_n("bkg_n","bkg_n", delmass, dm0,a1,a2,a3);

//Background1 Helicity
  RooChebychev bkgh1_n("bkgh1_n", "bkgh1_n",helicity, RooArgList(p));
  RooGenericPdf bkgh2_n("bkgh2_n","helicity*helicity",RooArgSet(helicity));
  RooAddPdf bkgH1_n("bkgH1_n","bkgH1_n",RooArgList(bkgh1_n,bkgh2_n),f2);

//Background1 D0 
  RooChebychev bkg1_D0_n("bkg1_D0_n", "bkg1_D0_n", dzero, RooArgList(p1));
  RooExponential bkg2_D0_n("bkg2_D0_n", "bkg2_D0_n", dzero, c1); 
  RooAddPdf Bkg_D0_n("Bkg_D0_n","Bkg_D0_n",RooArgList(bkg1_D0_n,bkg2_D0_n),RooArgList(f1));

//Product Background1
RooProdPdf Bkg1_N("Bkg1_N","Bkg1_N",RooArgSet(bkg_n,bkgH1_n,Bkg_D0_n));

//Background2 DeltaM
  //Gauss
  RooGaussian bkg_pdf_n("bkg_pdf_n", "signal Gaussian",delmass,mean,sigb);
 
  //Double-Gauss
  RooBifurGauss bkg_pdf1_n("bkg_pdf1_n", "signal Gaussian",delmass,mean,sig1b,sig2b);//***************8

  //Full Signal
  RooAddPdf bkg_tot2_n("bkg_tot2_n","bkg_tot2_n",RooArgList(bkg_pdf_n,bkg_pdf1_n),f_Bkg2);


//Background2 Helicity
  RooGenericPdf bkg2H_n("bkg2H_n","helicity*helicity",RooArgSet(helicity));

//Background2 D0
  RooCBShape Bkg_D0_2n("Bkg_D0_2n","Bkg_D0_2n", dzero, m1, s1, A1, n1); 

//Product Background2
RooProdPdf Bkg2_N("Bkg2_N","Bkg2_N",RooArgSet(bkg_tot2_n,bkg2H_n,Bkg_D0_2n));

//Full Background N
  RooAddPdf bkg_tot_n("bkg_tot_n","bkg_tot_n",RooArgList(Bkg1_N,Bkg2_N),F_Bkg);

//Full Model N
  RooAddPdf Model_N("Model_N","Model_N",RooArgList(Sig_N,bkg_tot_n),RooArgList(N_n,N_nb));



  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  RooDataSet combData("combData","combData",RooArgSet(helicity,delmass,dzero),Index(sample),Import("pos",*data_p),Import("neg",*data_n));
  
  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_P,"pos");
  simPdf.addPdf(Model_N,"neg");
  //_____________________________________________________

   RooFitResult* fitRes = simPdf.fitTo(combData,Save());
//   cout<<"var fit status = "<<fitRes->status()<<endl;

  //DeltaM PLOTING
  RooPlot *xframe_1 =delmass.frame(Bins(100),Title("#pi^{+}_{s}"));
  combData.plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("sig_tot_p"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("bkg_p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("bkg_tot2_p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));
   Model_P.paramOn(xframe_1,data_p);
  RooPlot *xframe_2 = delmass.frame(Bins(100),Title("#pi^{-}_{s}"));
  combData.plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("sig_tot_n"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("bkg_n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("bkg_tot2_n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));

     Model_N.paramOn(xframe_2,data_n);  
  TCanvas* c = new TCanvas("c","c",1200,800) ;

  c->Divide(2,3) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;


  //Helicity PLOTING
  RooPlot *yframe_1 =helicity.frame(Bins(100),Title("#pi^{+}_{s}"));
  combData.plotOn(yframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),ProjWData(sample,combData));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("sig_p"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("bkgH1_p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("bkg2H_p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));

//   Model_P.paramOn(yframe_1,data_p);

  RooPlot *yframe_2 = helicity.frame(Bins(100),Title("#pi^{-}_{s}"));
  combData.plotOn(yframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),ProjWData(sample,combData));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("sig_n"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("bkgH1_n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("bkg2H_n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));
//     Model_N.paramOn(yframe_2,data_n);  
//  TCanvas* c2 = new TCanvas("c2","c2",1200,800) ;

//  c2->Divide(2) ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; yframe_1->GetYaxis()->SetTitleOffset(1.4) ; yframe_1->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; yframe_2->GetYaxis()->SetTitleOffset(1.4) ; yframe_2->Draw() ;

  //D0 PLOTING
  RooPlot *zframe_1 =dzero.frame(Bins(100),Title("#pi^{+}_{s}"));
  combData.plotOn(zframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(zframe_1,Slice(sample,"pos"),ProjWData(sample,combData));
  simPdf.plotOn(zframe_1,Slice(sample,"pos"),Components("sig_dz_p"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(zframe_1,Slice(sample,"pos"),Components("Bkg_D0_p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(zframe_1,Slice(sample,"pos"),Components("Bkg_D0_2p"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));

//   Model_P.paramOn(yframe_1,data_p);

  RooPlot *zframe_2 = dzero.frame(Bins(100),Title("#pi^{-}_{s}"));
  combData.plotOn(zframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(zframe_2,Slice(sample,"neg"),ProjWData(sample,combData));
  simPdf.plotOn(zframe_2,Slice(sample,"neg"),Components("sig_dz_n"),ProjWData(sample,combData),LineStyle(kDashed));
  simPdf.plotOn(zframe_2,Slice(sample,"neg"),Components("Bkg_D0_n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGreen));
  simPdf.plotOn(zframe_2,Slice(sample,"neg"),Components("Bkg_D0_2n"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));
//     Model_N.paramOn(yframe_2,data_n);  
//  TCanvas* c2 = new TCanvas("c2","c2",1200,800) ;

//  c2->Divide(2) ;
  c->cd(5) ; gPad->SetLeftMargin(0.15) ; zframe_1->GetYaxis()->SetTitleOffset(1.4) ; zframe_1->Draw() ;
  c->cd(6) ; gPad->SetLeftMargin(0.15) ; zframe_2->GetYaxis()->SetTitleOffset(1.4) ; zframe_2->Draw() ;







  TCanvas* c2 = new TCanvas("c2","c2",1200,800) ;
  c2->Divide(2) ;
  c2->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  c2->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

  TCanvas* c3 = new TCanvas("c3","c3",1200,800) ;
  c3->Divide(2) ;
  c3->cd(1) ; gPad->SetLeftMargin(0.15) ; yframe_1->GetYaxis()->SetTitleOffset(1.4) ; yframe_1->Draw() ;
  c3->cd(2) ; gPad->SetLeftMargin(0.15) ; yframe_2->GetYaxis()->SetTitleOffset(1.4) ; yframe_2->Draw() ;

  TCanvas* c4 = new TCanvas("c4","c4",1200,800) ;
  c4->Divide(2) ;
  c4->cd(1) ; gPad->SetLeftMargin(0.15) ; zframe_1->GetYaxis()->SetTitleOffset(1.4) ; zframe_1->Draw() ;
  c4->cd(2) ; gPad->SetLeftMargin(0.15) ; zframe_2->GetYaxis()->SetTitleOffset(1.4) ; zframe_2->Draw() ;



//  c->SaveAs(name);
//  c3->SaveAs("FullFit.png");
//  delete c;
//  delete c2;
//  delete c3;



}

