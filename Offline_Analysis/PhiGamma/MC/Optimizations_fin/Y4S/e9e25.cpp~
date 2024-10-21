//Here we include D-> Phi Pi0 as Background
void e9e25(){

  int nsignal[101]={0},nback[101]={0},bestk;
  float cut,fom[101],bestfom=100;

  TFile *f = new TFile("PhiGamma_4s_FullMC.root");

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
if(ProbPi < 0.056){
if(ProbEta < 0.288){
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID

for(int k=1;k<=99;k++){
cut = (float)k*0.001+0.9;
if(E9E25 > cut){//Satisfies Pi0 Veto
if(Mygenfla==1){nsignal[k]++;}
if(Mygenfla!=1){nback[k]++;}
}
}//k loop

}}}}}}}}}

}

for(int k=1;k<=99;k++){
//fom[k]=(float)sqrt(nsignal[k]+nback[k])/nsignal[k];
cout<<"nsig = "<<nsignal[k]<<"\t nbkg = "<<nback[k]<<endl;
if(nsignal[k] != 0){fom[k]=(float)sqrt((nsignal[k]/20.0)+nback[k])/(nsignal[k]/20.0);}else{fom[k]==1000;}
cout<<"fom = "<<fom[k]<<endl;
if(fom[k] <= bestfom){bestfom=fom[k];bestk=k;}
}

cout<<fom[bestk]<<"\t"<<nsignal[bestk]/20.0<<"\t best veto = "<<(float)bestk*0.001+0.9<<endl;


}
