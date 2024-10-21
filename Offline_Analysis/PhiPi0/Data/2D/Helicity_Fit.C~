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
  RooRealVar helicity("helicity","cos #theta_{hel}",-1,1,"");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(helicity));
void fitit(int index);
void Helicity_Fit(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("PhiPi0_4s_Data.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm,f_dstcts,f_dstcharge,f_pi0mom,f_pi0mass,f_phimass,f_dstps,f_d0mass,f_hel1,f_hel2,f_piscost, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;

  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
  chain->SetBranchAddress("Kplushel",&f_hel1);
  chain->SetBranchAddress("Kminheli",&f_hel2);
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
//  for(int i=0;i<100;i++) 
    {
      chain->GetEntry(i);
		  
//Phot1cut
if(f_Phot1the > -60 && f_Phot1the < 67 && f_Photon1e > 0.05){photon1cutflag=1;}
else if((f_Phot1the < -60 || f_Phot1the > 67) && f_Photon1e > 0.1){photon1cutflag=1;}
else{photon1cutflag=0;}

//Phot2cut
if(f_Phot2the > -60 && f_Phot2the < 67 && f_Photon2e > 0.05){photon2cutflag=1;}
else if((f_Phot2the < -60 || f_Phot2the > 67) && f_Photon2e > 0.1){photon2cutflag=1;}
else{photon2cutflag=0;}

if(f_delm > 0.142 && f_delm < 0.149) {
//if(f_delm > 0.15) {
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
		  helicity.setVal(f_hel1);
		  nevt_ac_p++;
data->add(RooArgSet(helicity));

		}
	      else
		{
		  nevt_ac_n++;
		  helicity.setVal(f_hel2);
data->add(RooArgSet(helicity));

		}
}}}}}}}}}}}
    }



  //DEFINE PDF

  RooRealVar p("p","p",0);
  RooChebychev shape1("shape1", "shape1",helicity, RooArgList(p));
  RooRealVar shift1("a_{(shift)}","a_{(shift)}",0,-0.5,0.5);
  RooGenericPdf shape2("shape2","(helicity-a_{(shift)})*(helicity-a_{(shift)})",RooArgSet(helicity,shift1));
  RooRealVar f("frac_{cheby/x^{2}}","f_Sig2",0.6012,0,1);

   RooAddPdf model("model","model",RooArgList(shape1,shape2),f);

 
  //_____________________________________________________

   RooFitResult* fitRes = model.fitTo(*data);

  //PLOTING
  RooPlot* hel_frame = helicity.frame();
  hel_frame->SetTitle("Helicity Distribution");
  data->plotOn(hel_frame);

  model.plotOn(hel_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
  model.paramOn(hel_frame);
  model.plotOn(hel_frame,Components("shape1"),LineStyle(kDashed));


hel_frame->Draw();
/*
  //PLOTING
//Helicity signal projection
//  RooPlot *yframe_1 =helicity.frame(Bins(100),Title("D^{0} #rightarrow #phi #pi^{0}"));
  RooPlot *yframe = data->plotOn(helicity.frame(100),MarkerColor(kBlue));
//  combData.plotOn(yframe_1,Cut("sample==sample::pos"),CutRange("signal"));
//  simPdf.plotOn(yframe_1,Slice(sample,"pos"),ProjWData(sample,combData),ProjectionRange("signal"));
//  simPdf.plotOn(yframe_1,Slice(sample,"pos"),Components("bkgH_p"),ProjWData(sample,combData),LineStyle(kDashed),ProjectionRange("signal"));
//   Model_P.paramOn(yframe_1,data_p);

//  RooPlot *yframe_2 = helicity.frame(Bins(100),Title("D^{0} #rightarrow #phi #pi^{0}"));
//  combData.plotOn(yframe_2,Cut("sample==sample::neg"),CutRange("signal"));
//  simPdf.plotOn(yframe_2,Slice(sample,"neg"),ProjWData(sample,combData),ProjectionRange("signal"));
//  simPdf.plotOn(yframe_2,Slice(sample,"neg"),Components("bkgH_n"),ProjWData(sample,combData),LineStyle(kDashed),ProjectionRange("signal"));

  model.plotOn(yframe, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));




  TCanvas* c3 = new TCanvas("c3","c3",1200,800) ; c3->cd() ; 
  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;

*/

}


