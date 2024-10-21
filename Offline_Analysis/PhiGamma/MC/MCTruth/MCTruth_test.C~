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
void MCTruth_test(void)
{
    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");


    TCanvas *c4 = new TCanvas("myCanvas4","My Canvas4");
    TCanvas *c5 = new TCanvas("myCanvas5","My Canvas5");
    TCanvas *c6 = new TCanvas("myCanvas6","My Canvas6");


    TH1F *h111 = new TH1F("my_hist1","hist1",50,0.14,0.16);  //deltam
    TH1F *h2 = new TH1F("my_hist2","hist2",50,0.14,0.16);  //deltam
    TH1F *h3 = new TH1F("my_hist3","hist3",50,0.14,0.16);  //deltam
    TH1F *h4 = new TH1F("my_hist4","hist4",50,0.14,0.16);  //deltam
    TH1F *ha = new TH1F("my_hista","hista",50,0.14,0.16);  //deltam

    TH1F *h5 = new TH1F("my_hist5","hist5",50,1.65,2.05);  //d0 mass
    TH1F *h6 = new TH1F("my_hist6","hist6",50,1.65,2.05);  //d0 mass
    TH1F *h7 = new TH1F("my_hist7","hist7",50,1.65,2.05);  //d0 mass
    TH1F *h8 = new TH1F("my_hist8","hist8",50,1.65,2.05);  //d0 mass
    TH1F *hb = new TH1F("my_histb","histb",50,1.65,2.05);  //deltam

    TH1F *h9 = new TH1F("my_hist9","hist9",50,-1,1);  //helicity
    TH1F *h10 = new TH1F("my_hist10","hist10",50,-1,1);  //helicity
    TH1F *h11 = new TH1F("my_hist11","hist11",50,-1,1);  //helicity
    TH1F *h12 = new TH1F("my_hist12","hist12",50,-1,1);  //helicity
    TH1F *hc = new TH1F("my_histc","histc",50,-1,1);  //helicity


  h111->SetFillColor(4);
  h2->SetFillColor(8);
  h3->SetFillColor(2);
  h4->SetFillColor(6);

  h5->SetFillColor(4);
  h6->SetFillColor(8);
  h7->SetFillColor(2);
  h8->SetFillColor(6);

  h9->SetFillColor(4);
  h10->SetFillColor(8);
  h11->SetFillColor(2);
  h12->SetFillColor(6);

  ha->SetFillColor(5);
  hb->SetFillColor(5);
  hc->SetFillColor(5);

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

  chain->Add("FullGMC4S.root");
  
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t f_Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, f_helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, Piflag, Mypgenfl, Etaflag ;

  h1->SetBranchAddress("Deltam",&f_Deltam);
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
  h1->SetBranchAddress("Etaflag",&Etaflag);
  h1->SetBranchAddress("Dstarcha",&Dstarcha);

  THStack *hs1 = new THStack("hs1","M(D*) - M(D^{0})");
  THStack *hs2 = new THStack("hs2","D^{0} Mass");
  THStack *hs3 = new THStack("hs3","Helicity Distribution");
 

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);

if(f_Deltam > 0.14 && f_Deltam < 0.16){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.055){
if(ProbEta < 0.289){
if(E9E25 > 0.944){
if(Gamener > 0.6){
if(Dstp > 2.57){


if(Mygenfla==1){h111->Fill(f_Deltam); h5->Fill(MD0); h9->Fill(f_helicity);}
else if(Etaflag==1){h2->Fill(f_Deltam); h6->Fill(MD0); h10->Fill(f_helicity);}
else if(Mypgenfl==1){h3->Fill(f_Deltam); h7->Fill(MD0); h11->Fill(f_helicity);}
else if(Piflag==1){ha->Fill(f_Deltam); hb->Fill(MD0); hc->Fill(f_helicity);}
else{h4->Fill(f_Deltam); h8->Fill(MD0); h12->Fill(f_helicity);}

}}}}}}}}
    }

  chain->Add("FullGMC5S.root");
  Int_t nevt1=(int)chain->GetEntries();
  for(int i=nevt;i<nevt1;i++) 
    {
      chain->GetEntry(i);

if(f_Deltam > 0.14 && f_Deltam < 0.16){
if(MD0 > 1.68 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.056){
//if(ProbEta < 0.285){
if(E9E25 > 0.947){
if(Gamener > 0.59){
if(Dstp > 3.1){


if(Mygenfla==1){h111->Fill(f_Deltam); h5->Fill(MD0); h9->Fill(f_helicity);}
else if(Etaflag==1){h2->Fill(f_Deltam); h6->Fill(MD0); h10->Fill(f_helicity);}
else if(Piflag==1){h3->Fill(f_Deltam); h7->Fill(MD0); h11->Fill(f_helicity);}
else{h4->Fill(f_Deltam); h8->Fill(MD0); h12->Fill(f_helicity);}


}}}}}}}//}
    }




  hs1->Add(h4);
  hs1->Add(h2);
  hs1->Add(h3);
  hs1->Add(h111);
  hs1->Add(ha);


  hs2->Add(h8);
  hs2->Add(h6);
  hs2->Add(h7);
  hs2->Add(h5);
  hs2->Add(hb);

  hs3->Add(h12);
  hs3->Add(h10);
  hs3->Add(h11);
  hs3->Add(h9);
  hs3->Add(hc);

  leg = new TLegend(0.6,0.7,0.89,0.89);
//  leg = new TLegend();
  leg->SetHeader("Reconstructed Object");
  leg->AddEntry(h111,"Signal","f");
  leg->AddEntry(h2,"#eta contamination","f");
  leg->AddEntry(h3,"#pi^{0}  contamination","f");
  leg->AddEntry(h4,"combinatorial bkg","f");
  leg->AddEntry(ha,"new one","f");

c1->cd(); hs1->Draw();  leg->Draw();
c2->cd(); hs2->Draw();  leg->Draw();
c3->cd(); hs3->Draw();  leg->Draw();

c4->cd(); ha->Draw();  
c5->cd(); hb->Draw();  
c6->cd(); hc->Draw();  




}

