void MD(){

  int nsignal[121][201]={0},nback[121][201]={0},bestk=100,bestk2=100;
  float cut,cut2,fom[121][201],bestfom=100,x[121][201],y[121][201];

  TFile *f = new TFile("FullGMC.root");

// gStyle->SetOptStat(0);
//   gStyle->SetOptFit();
//gStyle->SetPalette(1);

//TGraph2D *g = new TGraph2D();

  TTree *t1 = (TTree*)f->Get("h1");
 float Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25;

  t1->SetBranchAddress("Deltam",&Deltam);
  t1->SetBranchAddress("Mygenfla",&Mygenfla);
  t1->SetBranchAddress("Dstps",&Dstp);
  t1->SetBranchAddress("Kplushel",&helicity);
  t1->SetBranchAddress("Gamenerg",&Gamener);
  t1->SetBranchAddress("P_piz",&P_piz);
  t1->SetBranchAddress("P_eta",&P_eta);
  t1->SetBranchAddress("Vetom1",&Vetom1);
  t1->SetBranchAddress("Vetom2",&Vetom2); 
  t1->SetBranchAddress("P_piz_2",&ProbPi);
  t1->SetBranchAddress("P_eta_2",&ProbEta);
  t1->SetBranchAddress("Phimass",&MPhi);
  t1->SetBranchAddress("Dzeromas",&MD0);
  t1->SetBranchAddress("Ga_e9e25",&E9E25);

/*for(int k=1;k<=10;k++){    // for E(Gam)
for(int k2=1;k2<=10;k2++){ //
cut = (float)k*0.01+1.86;  //upper cut
cut2 = 1.86-((float)(k2*0.01));  //lower cut
cout<<cut<<"\t"<<cut2<<endl;
}}*/
  Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    t1->GetEntry(i);

//if(Deltam > 0.14 && Deltam < 0.16){
if(fabs(Deltam-0.1454) < 0.0017) {
if(MD0 > 1.65 && MD0 < 2.05){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.055){
if(ProbEta < 0.289){
if(E9E25 > 0.945){
if(Gamener > 0.65 && Dstp > 2.8){//Satisfies Pi0 Veto


for(int k=1;k<=10;k++){    // for E(Gam)
for(int k2=1;k2<=10;k2++){ //for P*(D*)
cut = (float)k*0.01;  //upper cut
cut  = cut+1.86;
cut2 = (float)(k2*0.01);  //lower cut
cut2  = 1.86-cut2;
if(MD0 < cut && MD0 > cut2){//Satisfies Pi0 Veto
if(Mygenfla==1){nsignal[k][k2]++;}
if(Mygenfla!=1){nback[k][k2]++;}
}
}}//k-k2 loops

}}}}}}}

}
int count=0;

for(int k=1;k<=10;k++){    
for(int k2=1;k2<=10;k2++){ 
if(nsignal[k][k2]!=0){
//fom[k][k2]=(float)sqrt(nsignal[k][k2]+nback[k][k2])/nsignal[k][k2];
fom[k][k2]=(float)sqrt((nsignal[k][k2]/20)+nback[k][k2])/(nsignal[k][k2]/20);
//cout<<"fom is "<<fom[k][k2]<<endl;
count++;
if(fom[k][k2] <= bestfom){bestfom=fom[k][k2];bestk=k; bestk2=k2;}
}//if
}}

cut = (float)bestk*0.01;  //upper cut
cut  = cut+1.86;
cut2 = (float)(bestk2*0.01);  //lower cut
cut2  = 1.86-cut2;


//cout<<bestk<<"\t"<<bestk2<<endl;
cout<<fom[bestk][bestk2]<<"\t"<<nsignal[bestk][bestk2]/20<<"\t"<<nback[bestk][bestk2]<<"\t upper cut = "<<cut<<"MeV\t lower cut = "<<cut2<<"GeV"<<endl;



//  g->Draw("TRI1");
//  g->Draw("SAME P0");
//  g->Draw("surf1");
//  g->GetZaxis()->SetRangeUser(tempminz,tempmaxz);





}
