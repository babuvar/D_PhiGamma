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
  RooRealVar delmass("delmass","delmass",0.14,0.16,"GeV");
  RooDataSet* data=new RooDataSet("data_p","data_p",RooArgSet(delmass));
void deltaM(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("DtoPhiGamma_500Kx2.root");
//  chain->Add("DtoPhiPi0_1Mx2.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, f_helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, Piflag, Mypgenfl,f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz ;

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
		  delmass.setVal(Deltam);

if(Deltam > 0.14 && Deltam < 0.16){
if(MD0 < 1.89 && MD0 > 1.84){
//if(MD0 > 1.68 && MD0 < 2.05){
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

data->add(RooArgSet(delmass));


}}}}}}}}}}}}
    }

  chain->Add("DtoPhiGamma_80760x2.root");
//  chain->Add("DtoPhiPi0_161520x2.root");
  Int_t nevt1=(int)chain->GetEntries();
  cout<<"nevt2\t"<<nevt <<endl;
  for(int i=nevt;i<nevt1;i++) 
    {
      chain->GetEntry(i);
		  delmass.setVal(Deltam);

if(Deltam > 0.14 && Deltam < 0.16){
//if(MD0 > 1.68 && MD0 < 2.05){
if(MD0 < 1.89 && MD0 > 1.84){
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

data->add(RooArgSet(delmass));


}}}}}}}}}}}}
    }


  //DEFINE PDF
  //COMMON 
  RooRealVar mean("#mu","mean",0.145428,0.14,0.15);
  RooRealVar sig("#sigma","sig",0.000397327,0,0.0015);//sigma);//0.000397327,0,0.003);
  RooRealVar sig1("#sigma_{L}","sig1",0.00089779,0.0004,0.0050);
  RooRealVar sig2("#sigma_{R}","sig2",0.00109104,0,0.0050);
  RooRealVar f_Sig("frac_{gauss/bifur}","f_Sig",0.6012,0.0,1);
  RooRealVar dm0("dm0", "threshold value",0.13957018);
  RooRealVar a1("#alpha_{threshold}", "coef 1",0.002,0.0,0.015);
  RooRealVar a2("a2", "coef 2",0.0);
  RooRealVar a3("#beta_{threshold}", "coef 3",-0.119,-20.0,20);
  RooRealVar N_sig("N_{sig}","N_sig",100,0,1000000);
  RooRealVar N_bkg("N_{bkg}","N_bkg",100,0,1000000);
//SIGNAL

//Signal 
  //Gauss
  RooGaussian sig_pdf("sig_pdf", "signal Gaussian",delmass,mean,sig);
 
  //Double-Gauss
  RooBifurGauss sig_pdf1("sig_pdf1", "signal Gaussian",delmass,mean,sig1,sig2);//***************8

  //Full Signal
  RooAddPdf sig_tot("sig_tot","sig_tot",RooArgList(sig_pdf,sig_pdf1),f_Sig);

//Background
  RooDstD0BG bkg("bkg","background", delmass, dm0,a1,a2,a3);

//Full Model
  RooAddPdf Model("Model","",RooArgList(sig_tot,bkg),RooArgList(N_sig,N_bkg));

//Fit
   RooFitResult* fitRes = Model.fitTo(*data);
//  RooPlot* delmass_frame = new RooPlot("varg","varg");

  RooPlot* delmass_frame = delmass.frame();
delmass_frame->SetTitle("M_{D*} - M_{D}");
  data.plotOn(delmass_frame);

  Model.plotOn(delmass_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
  Model.paramOn(delmass_frame);



delmass_frame->Draw();


/*
  RooPlot *yframe1 = data->plotOn(delmass.frame(100),MarkerColor(kBlue));

  TGaxis::SetMaxDigits(3);
  yframe1->Draw();
*/

}

