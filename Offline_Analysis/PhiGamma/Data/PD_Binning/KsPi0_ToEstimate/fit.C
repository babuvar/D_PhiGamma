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
#include "TText.h"
#include "TArrow.h"
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
//#include <iostream>
using namespace RooFit;
using namespace std;
void fitit(RooRealVar delmG, RooDataSet * data_pG, RooDataSet * data_nG, int num);
TCanvas* c = new TCanvas("c","c",1200,600) ;
TCanvas* c1 = new TCanvas("c1","c1",1200,600) ;
TCanvas* c3 = new TCanvas("c3","c3",1200,600) ;
TCanvas* c4 = new TCanvas("c4","c4",1200,600) ;

  RooRealVar delm("delm","#DeltaE (GeV)", 0, 10);
RooDataSet* data_p[50];
RooDataSet* data_n[50];
float Araw[50]={0},EAraw[50]={0},Cost[50]={0},ECost[50]={0},Acp[25]={0},Afb[25]={0},EAfb[25]={0},Cost2[25]={0};
//______________________________________________________________________
float n=2,m; // number of bins is 2n
//______________________________________________________________________
m=1/n;
int i;
void fit(void)
{
for(i=0;  i<=49; i++){
data_p[i] = new RooDataSet("data_p","data_p",RooArgSet(delm));
data_n[i] = new RooDataSet("data_n","data_n",RooArgSet(delm));

} 


  //LOAD DATA FILE
  TChain* chain=new TChain("h4");

  chain->Add("data.root");
  
  Int_t nevt=(int)chain->GetEntries();
  RooRealVar delm("delm","#DeltaE (GeV)", 0.140, 0.160);
  float f_delm,f_dstcts,f_dstcharge;

  chain->SetBranchAddress("delmass",&f_delm);
  chain->SetBranchAddress("dst_ct_c",&f_dstcts);
  chain->SetBranchAddress("spi_char",&f_dstcharge);



float bin;int bin1; 
  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
    {
      chain->GetEntry(i);
          delm.setVal(f_delm);
bin=floor(f_dstcts/m)+n;//cout<<bin<<",";
bin1=bin;

if(f_dstcharge==+1){data_p[bin1]->add(RooArgSet(delm));}
else if(f_dstcharge==-1){data_n[bin1]->add(RooArgSet(delm));}

    }

//c->cd();
//  RooPlot* frame1 =delm.frame();
//  data_p[11]->plotOn(frame1); frame1->Draw() ;
//  data_n[11]->plotOn(frame1); frame1->Draw() ;

for(i=0;i<2*n;i++)
{fitit(delm,data_p[i],data_n[i],i);}

int k1;

float fb=0,cp=0,er;
ofstream fout;
fout.open("AFB.txt");

//TGraph *g1 = new TGraph();
for(i=0;i<n;i++){
EAfb[i]=0;//*
k1=(2*n)-(i+1);
Afb[i]=(Araw[i]-Araw[k1])/2;
Acp[i]=(Araw[i]+Araw[k1])/2;
//EAfb[i]=sqrt((EAraw[i]*EAraw[i])+(EAraw[k1]*EAraw[k1]));
EAfb[i]=EAfb[i]+(1/(EAraw[i]*EAraw[i]))+(1/(EAraw[k1]*EAraw[k1]));//*
EAfb[i]=sqrt(EAfb[i]);//*
EAfb[i]=1/EAfb[i];//*

fb=fb+(Afb[i]/(EAfb[i]*EAfb[i]));
cp=cp+(Acp[i]/(EAfb[i]*EAfb[i]));
er=er+(1/(EAfb[i]*EAfb[i]));
fout<<Afb[i]<<endl;
cout<<"Afb["<<i<<"] = "<<Afb[i]<<endl;
cout<<"Acp["<<i<<"] = "<<Acp[i]<<endl;
cout<<"EAfb["<<i<<"] = "<<EAfb[i]<<endl;
}
 fout.close();

er=sqrt(er);
er=1/er;
cp=cp*er*er;
fb=fb*er*er;

cout<<"AFB = ("<<fb*100<<" +/- "<<er*100<<") %"<<endl;
cout<<"ACP = ("<<cp*100<<" +/- "<<er*100<<") %"<<endl;
fout.open("Result.txt");
fout<<"AFB = ("<<fb*100<<" +/- "<<er*100<<") %"<<endl;
fout<<"ACP = ("<<cp*100<<" +/- "<<er*100<<") %"<<endl;
 fout.close();

c1->cd();
   TGraphErrors *g1 = new TGraphErrors(2*n,Cost,Araw,ECost,EAraw);
   g1->SetTitle("Measured A_{RAW} vs Cos #Theta *_{D*}");
   g1->SetMarkerColor(kRed);
   g1->SetMarkerStyle(21);
g1->Draw("AP");


c3->cd();
   TGraphErrors *g2 = new TGraphErrors(n,Cost2,Afb,ECost,EAfb);
//   TGraphErrors *g1 = new TGraphErrors(2*n,Cost,Araw,ECost,EAraw);
   g2->SetTitle("Measured A_{FB} vs |Cos #Theta *_{D*}|");
   g2->SetMarkerColor(kBlue);
   g2->SetMarkerStyle(21);
g2->Draw("AP");


c4->cd();
   TGraphErrors *g3 = new TGraphErrors(n,Cost2,Acp,ECost,EAfb);
   g3->SetTitle("Measured A_{CP} vs |Cos #Theta *_{D*}|");
   g3->SetMarkerColor(kGreen);
   g3->SetMarkerStyle(21);
g3->Draw("AP");

exit();
}









void fitit(RooRealVar delmG, RooDataSet * data_pG, RooDataSet * data_nG, int num)
  {

  //DEFINE PDF
  //____________________________________________________
  //CHARGE=+1
  RooRealVar meanG("meanG","meanG",0.1455,0.14,0.16);
  RooRealVar sigG("sigG","sigG",0.00043,0,0.01);
  RooGaussian sig_pdf_pG("sig_pdf_pG", "signal Gaussian G",delmG,meanG,sigG);

  RooRealVar sig1G("sig1G","sig1G",0.001,0,0.0015);
  RooRealVar sig2G("sig2G","sig2G",0.001,0,0.0015);

//  RooRealVar frac1("frac1","Fraction 1",fr1);
//  RooRealVar frac2("frac2","Fraction 2",fr2);

//  RooFormulaVar sig1G("sig1G","frac1*sigG",RooArgList(frac1,sigG));
//  RooFormulaVar sig2G("sig2G","frac2*sigG",RooArgList(frac2,sigG));

  RooBifurGauss sig_pdf1_pG("sig_pdf1_pG", "signal Gaussian G",delmG,meanG,sig1G,sig2G);
  
  RooRealVar dm0G("dm0G", "threshold valueG",0.13957018);
  RooRealVar a1G("a1G", "coef 1G",0.2,0.0,10.0);
  RooRealVar a2G("a2G", "coef 2G",0.0);
  RooRealVar a3G("a3G", "coef 3G",-0.119,-10.0,10.0);
  RooDstD0BG bkg_pG("bkg_pG","backgroundG", delmG, dm0G,a1G,a2G,a3G);

//  RooRealVar f_SigG("f_SigG","f_SigG",v6);
  RooRealVar f_SigG("f_SigG","f_SigG",0.4,0,1);
//  RooRealVar f_SigG("f_SigG","f_SigG",0.6,0,1);
  RooAddPdf sig_tot_pG("sig_tot_pG","sig_tot_pG",RooArgList(sig_pdf_pG,sig_pdf1_pG),f_SigG);
  
  //____________________________________________________
  RooRealVar ArawG("ArawG","ArawG",0,-1,1);
  RooRealVar N_tG("N_tG","N_tG",400000,0,1000000);
  RooFormulaVar N_nG("N_nG","(0.5)*(1-ArawG)*N_tG",RooArgList(ArawG,N_tG));
  RooFormulaVar N_pG("N_pG","(0.5)*(1+ArawG)*N_tG",RooArgList(ArawG,N_tG));
  //____________________________________________________
  RooRealVar AbkgG("AbkgG","AbkgG",0,-1,1);
  RooRealVar N_tbG("N_tbG","N_tbG",100000,0,800000);
  RooFormulaVar N_nbG("N_nbG","(0.5)*(1-AbkgG)*N_tbG",RooArgList(AbkgG,N_tbG));
  RooFormulaVar N_pbG("N_pbG","(0.5)*(1+AbkgG)*N_tbG",RooArgList(AbkgG,N_tbG));
  //_________________________________fdelm.____________________


  RooAddPdf model_pG("model_pG", "ModelG",RooArgList(sig_tot_pG,bkg_pG), RooArgList(N_pG,N_pbG));
  //_____________________________________________________
  //CHARGE=-1

  RooGaussian sig_pdf_nG("sig_pdf_nG", "signal GaussianG",delmG,meanG,sigG);
  RooBifurGauss sig_pdf1_nG("sig_pdf1_nG", "signal GaussianG",delmG,meanG,sig1G,sig2G);
  RooDstD0BG bkg_nG("bkg_nG","backgroundG", delmG, dm0G,a1G,a2G,a3G);
  RooAddPdf sig_tot_nG("sig_tot_nG","sig_tot_nG",RooArgList(sig_pdf_nG,sig_pdf1_nG),f_SigG);
  RooAddPdf model_nG("model_nG", "ModelG",RooArgList(sig_tot_nG,bkg_nG), RooArgList(N_nG,N_nbG));
   
  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sampleG("sampleG","sampleG");
  sampleG.defineType("posG");
  sampleG.defineType("negG");
  RooDataSet combDataG("combDataG","combDataG",delmG,Index(sampleG),Import("posG",*data_pG),Import("negG",*data_nG));
  
  RooSimultaneous simPdfG("simPdfG","simPdfG",sampleG);
  simPdfG.addPdf(model_pG,"posG");
  simPdfG.addPdf(model_nG,"negG");
  //_____________________________________________________

  simPdfG.fitTo(combDataG);

  //PLOTING
  RooPlot *xframe_1G =delmG.frame(Bins(100),Title("D*^{+}_{s}"));
  combDataG.plotOn(xframe_1G,Cut("sampleG==sampleG::posG"));
  simPdfG.plotOn(xframe_1G,Slice(sampleG,"posG"),ProjWData(sampleG,combDataG));
  simPdfG.plotOn(xframe_1G,Slice(sampleG,"posG"),Components("bkg_pG"),ProjWData(sampleG,combDataG),LineStyle(kDashed));
   model_pG.paramOn(xframe_1G,data_pG);
//   model_pG.paramOn(xframe_1G,Layout(0.85));
  RooPlot *xframe_2G = delmG.frame(Bins(100),Title("D*^{-}_{s}"));
  combDataG.plotOn(xframe_2G,Cut("sampleG==sampleG::negG"));
  simPdfG.plotOn(xframe_2G,Slice(sampleG,"negG"),ProjWData(sampleG,combDataG));
  simPdfG.plotOn(xframe_2G,Slice(sampleG,"negG"),Components("bkg_nG"),ProjWData(sampleG,combDataG),LineStyle(kDashed));
 //    model_n.paramOn(xframe_2,data_nS);  
      TCanvas* c2 = new TCanvas("c2","c2",1200,600) ;

  c2->Divide(2) ;

  c2->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1G->GetYaxis()->SetTitleOffset(1.4) ; xframe_1G->Draw() ;
  c2->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2G->GetYaxis()->SetTitleOffset(1.4) ; xframe_2G->Draw() ;

//char *b=new char[name2.size()+1]
//b[name2.size()]=0;
//memcpy(b,name2.c_str(),name2.size());

  char buffer [50]; int n1;
  n1=sprintf (buffer, "Fit%d.png",num);
  c2->SaveAs(buffer);

Araw[num]=ArawG.getVal();
EAraw[num]=ArawG.getError();
Cost[num]=-1+(m*num)+(m/2);
if(num<n){Cost2[num]=(m*num)+(m/2);}
  delete c2;

//  c[trial][cnt]->SetWindowSize(1200,600); 
//  c[trial][cnt]->Divide(2) ;
//  c[trial][cnt]->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1G->GetYaxis()->SetTitleOffset(1.4) ; xframe_1G->Draw() ;
//  c[trial][cnt]->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2G->GetYaxis()->SetTitleOffset(1.4) ; xframe_2G->Draw() ;


/*
//Getting the parameters
  RooArgSet* params = simPdfG.getVariables() ;
  RooRealVar* paraw = (RooRealVar*) params->find("ArawG") ;

//Float_t asym=p1->getVal();p1->getError();
ofstream fout;
fout.open("Asym1.txt",ios::app);
fout<<paraw->getVal()<<"\t"<<paraw->getError()<<endl;
fout.close();

A[trial][cnt]=paraw->getVal();
B[trial][cnt]=paraw->getError();
*/
}
















void show(Float_t C[11],Float_t D[11])
{
Int_t i,j;

ofstream fout;
fout.open("Asym.txt",ios::app);
for(i=1;i<=10;i++){
fout<<C[i]<<"\t"<<D[i]<<endl;}
fout.close();
}

void bestval(Float_t A[9][11], Float_t B[9][11],Float_t C[11],Float_t D[11])
{
Int_t i,j,tempj;
Float_t tempB;
for(i=1;i<=10;i++){
tempB=10; tempj=1;
for(j=1;j<=8;j++){
if(B[j][i]<=tempB && B[j][i]>=0.001){tempB=B[j][i]; tempj=j;}
}
C[i]=A[tempj][i]; D[i]=B[tempj][i]; BJ[i]=tempj;
}
}//void bestval


void getasy(Float_t C[11],Float_t D[11])
{
Int_t i=0;
Float_t M=0,S=0,D1=0;
for(i=1;i<=5;i++){
C[i]=(C[i]+C[11-i])/2;
D[i]=(1/(D[i]*D[i]))+(1/(D[11-i]*D[11-i]));
D[i]=sqrt(D[i]);
D[i]=1/D[i];
//D[i]=(D[i]*D[i])+(D[11-i]*D[11-i]);
//D[i]=sqrt(D[i]);
D1=D1+(1/(D[i]*D[i]));
M=M+(C[i]/(D[i]*D[i]));
}
S=sqrt(D1);
S=1/S;
M=M/D1;
ofstream fout;
fout.open("Asym.txt",ios::app);
fout<<endl<<"Combined ACP = "<<M<<" Combined Error on ACP = "<<S<<endl;
fout.close();

}

void savecan(TCanvas *c[9][11], Int_t BJ[11])
{//Saving it
//  TCanvas c2;
  char buffer [50]; int n;
Int_t i,j;
for(i=1;i<=10;i++){
  j=BJ[i];
  n=sprintf (buffer, "Fit%dCanvas.png",i);
  c[j][i]->SaveAs(buffer);
}

}

void initcan(TCanvas *c[9][11])
{
  char buffer [50];
int n;

for(int i=1;i<=8;i++){
for(int j=1;j<=10;j++){
  n=sprintf (buffer, "Canvas[%d][%d]", i, j);
  c[i][j]= new TCanvas(buffer);
//  c[i][j]->SetWindowSize(1200,600); 
//  c[i][j]->Divide(2) ;

}}

}

  void delcan(TCanvas *c[9][11])
{
for(int i=1;i<=8;i++){
for(int j=1;j<=10;j++){
delete c[i][j];
}}

}






