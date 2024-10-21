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
  RooRealVar dzero("dzero","dzero",1.68,2.05,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
void Dzero(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("DtoPhiGamma_500Kx2.root");
//  chain->Add("DtoPhiPi0_1Mx2.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, f_helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, Piflag, Mypgenfl,f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz;

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
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
		  dzero.setVal(MD0);

if(Deltam > 0.14 && Deltam < 0.16){
//if(Deltam > 0.142 && Deltam < 0.148){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.062){
if(ProbEta < 0.363){
if(E9E25 > 0.938){
if(Gamener > 0.6){
if(Dstp > 2.55){
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID

data->add(RooArgSet(dzero));


}}}}}}}}}}}}
    }
  chain->Add("DtoPhiGamma_80760x2.root");
//  chain->Add("DtoPhiPi0_161520x2.root");

  Int_t nevt1=(int)chain->GetEntries();
  cout<<"nevt2\t"<<nevt1 <<endl;
  for(int i=nevt;i<nevt1;i++) 
    {
      chain->GetEntry(i);
		  dzero.setVal(MD0);

if(Deltam > 0.14 && Deltam < 0.16){
//if(Deltam > 0.142 && Deltam < 0.148){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.065){
if(ProbEta < 0.37){
if(E9E25 > 0.939){
if(Gamener > 0.490){
if(Dstp > 3.1){
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID

data->add(RooArgSet(dzero));


}}}}}}}}}}}}
    }


  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,1000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,5000);
//Signal
  RooRealVar m("#mu","m",1.86,1.83,1.88);//
  RooRealVar s("#sigma","s",0.01,0,0.04);//
  RooRealVar a("#alpha","a",0.5,0,2);
  RooRealVar n("n","n",2,0,150);


  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n); 
//Background
  RooRealVar p("p","p",0.0);   
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooRealVar c("c","c",-9,-10,0);
  RooExponential bkg2("bkg2", "bkg2", dzero, c);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
  RooAddPdf Bkg("Bkg","Bkg",RooArgList(bkg1,bkg2),RooArgList(f));


//Full Model
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));

//Fit
   RooFitResult* fitRes = Model.fitTo(*data);
//   RooFitResult* fitRes = Sig.fitTo(*data);
  RooPlot* dzero_frame = dzero.frame();
  dzero_frame->SetTitle("Mass of D candidate");


  data->plotOn(dzero_frame);

  Model.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
  Model.paramOn(dzero_frame);
//  Sig.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
//  Sig.paramOn(dzero_frame);


dzero_frame->Draw();


/*
  RooPlot *yframe1 = data->plotOn(dzero.frame(100),MarkerColor(kBlue));

  TGaxis::SetMaxDigits(3);
  yframe1->Draw();
*/

}

