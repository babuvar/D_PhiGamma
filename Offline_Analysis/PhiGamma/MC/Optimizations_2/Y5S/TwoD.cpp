//Here we include D-> Phi Pi0 as Background
void TwoD(){

//  int nsignal[121][201]={0},nback[121][201]={0},bestk,bestk2;
//  float cut,cut2,fom[121][201],bestfom=100,x[121][201],y[121][201];
  int nsignal[300][300]={0},nback[300][300]={0},bestk,bestk2;
  float cut,cut2,fom[300][300],bestfom=100,x[300][300],y[300][300];


  TFile *f = new TFile("PhiGamma_5s_FullMC.root");

 gStyle->SetOptStat(0);
   gStyle->SetOptFit();
gStyle->SetPalette(1);

TGraph2D *g = new TGraph2D();
	g->SetTitle("2-D Optimization for #Upsilon(5S)");

  TTree *t1 = (TTree*)f->Get("h1");
 float Deltam, Dzeromas, Dstarmas, Mygenfla, Dstarcha, Dstp, helicity, Gamener, P_piz, P_eta, Vetom1, Vetom2, ProbPi, ProbEta, MPhi, MD0 ,  E9E25, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz;

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
  t1->SetBranchAddress("Kppid",&f_Kppid);
  t1->SetBranchAddress("Kmpid",&f_Kmpid);
  t1->SetBranchAddress("Pispid",&f_Pispid);
  t1->SetBranchAddress("Kpdr",&f_Kpdr);
  t1->SetBranchAddress("Kmdr",&f_Kmdr);
  t1->SetBranchAddress("Pisdr",&f_Pisdr);
  t1->SetBranchAddress("Kpdz",&f_Kpdz);
  t1->SetBranchAddress("Kmdz",&f_Kmdz);
  t1->SetBranchAddress("Pisdz",&f_Pisdz);


  Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    t1->GetEntry(i);

if(Deltam > 0.142 && Deltam < 0.149) {
if(MD0 > 1.8 && MD0 < 1.9){
if(MPhi > 1.01 && MPhi < 1.03){
if(ProbPi < 0.065){
if(ProbEta < 0.37){
if(E9E25 > 0.939){
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID

for(int k=0;k<=200;k++){    // for E(Gam)
for(int k2=0;k2<=100;k2++){ //for P*(D*)
cut = (float)k*0.01+0.2;  //200 MeV
cut2 = (float)k2*0.01+3.1;  //2 GeV
if(Gamener > cut && Dstp > cut2){//Satisfies Pi0 Veto
if(Mygenfla==1){nsignal[k][k2]++;}
if(Mygenfla!=1){nback[k][k2]++;}
}
}}//k-k2 loops

}}}}}}}}}}

}
int count=0;

for(int k=0;k<=200;k++){    
for(int k2=0;k2<=100;k2++){ 
//fom[k][k2]=(float)sqrt(nsignal[k][k2]+nback[k][k2])/nsignal[k][k2];
fom[k][k2]=(float)sqrt((nsignal[k][k2]/100.0)+(nback[k][k2]/6.0))/(nsignal[k][k2]/100.0);
x[k][k2]=(float)k*10+200;
y[k][k2]=(float)k2*0.01+3.1;
g->SetPoint(count,x[k][k2],y[k][k2],fom[k][k2]);
count++;
if(fom[k][k2] <= bestfom){bestfom=fom[k][k2];bestk=k; bestk2=k2;}
}}

cout<<fom[bestk][bestk2]<<"\t"<<nsignal[bestk][bestk2]/100.0<<"\t"<<nback[bestk][bestk2]/20<<"\t best E(Gam) = "<<(float)bestk*10+200<<"MeV\t best P*(D*) = "<<(float)bestk2*0.01+3.1<<"GeV"<<endl;


TCanvas *c1 = new TCanvas("c","2-D Optimization");
c1->cd();

gStyle->SetPalette(1);

//  g->Draw("TRI1");
//  g->Draw("SAME P0");
//  g->Draw("surf1");
  g->Draw("colz");
//  g->GetZaxis()->SetRangeUser(tempminz,tempmaxz);





}

//0.26972	66.36	76	 best E(Gam) = 590MeV	 best P*(D*) = 3.1GeV


