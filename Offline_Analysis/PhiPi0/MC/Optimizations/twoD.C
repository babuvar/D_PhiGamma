void twoD()
{

TCanvas *c1 = new TCanvas("c","2-D Optimization",0,0,700,600);
TGraph2D *dt = new TGraph2D();
	dt->SetTitle("2-D Optimization");

  TChain* chain=new TChain("h1");

  chain->Add("PhiPi0_4s_GMC_0.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;
  Float_t f_delm, f_kid, f_pi0mom, f_pi0mass, f_phimass, f_d0mass, f_dstps, f_dr, f_dz, f_phihel, f_dstcts, f_dstcharge, f_mygenflag, NS[500][50]={0}, NB[500][50]={0} ,sig[500][50],i,j,bestsig=1000000,besti=3,bestj=20, f_Kppid,f_Kmpid,f_Pispid,f_Kpdr,f_Kmdr,f_Pisdr,f_Kpdz,f_Kmdz,f_Pisdz,f_Phot1the,f_Phot2the,f_Photon1e,f_Photon2e;

  chain->SetBranchAddress("Deltam",&f_delm);
  chain->SetBranchAddress("Pizmomen",&f_pi0mom);
  chain->SetBranchAddress("Pizmassn",&f_pi0mass);
  chain->SetBranchAddress("Phimass",&f_phimass);
  chain->SetBranchAddress("Dzeromas",&f_d0mass);
  chain->SetBranchAddress("Dstps",&f_dstps);
  chain->SetBranchAddress("Kplushel",&f_phihel);
  chain->SetBranchAddress("Dstcts",&f_dstcts);
  chain->SetBranchAddress("Dstarcha",&f_dstcharge);
  chain->SetBranchAddress("Mygenfla",&f_mygenflag);
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



  Int_t flag;
  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);

int photon1cutflag=0, photon2cutflag=0;

  for(int ii=0;ii<nevt;ii++) 
    {
      flag=0;
      chain->GetEntry(ii);

//Phot1cut
if(f_Phot1the > -60 && f_Phot1the < 67 && f_Photon1e > 0.05){photon1cutflag=1;}
else if((f_Phot1the < -60 || f_Phot1the > 67) && f_Photon1e > 0.1){photon1cutflag=1;}
else{photon1cutflag=0;}

//Phot2cut
if(f_Phot2the > -60 && f_Phot2the < 67 && f_Photon2e > 0.05){photon2cutflag=1;}
else if((f_Phot2the < -60 || f_Phot2the > 67) && f_Photon2e > 0.1){photon2cutflag=1;}
else{photon2cutflag=0;}

if(fabs(f_delm-0.1454) < 0.0017) {
if(f_phimass > 1.01 && f_phimass < 1.03) {
if(f_pi0mass > 0.119 && f_pi0mass < 0.151) {
if(f_d0mass > 1.83 && f_d0mass < 1.89) {
if(fabs(f_Kpdr) < 1.0 && fabs(f_Kmdr) < 1.0 && fabs(f_Pisdr) < 1.0){//dr
if(fabs(f_Kpdz) < 3.0 && fabs(f_Kmdz) < 3.0 && fabs(f_Pisdz) < 3.0){//dz
if(fabs(f_Kppid) > 0.1 && fabs(f_Kmpid) > 0.1){//K-PID
if(fabs(f_Pispid) < 0.9){//Pi-PID
if(photon1cutflag == 1 && photon1cutflag == 1){

Float_t a,b,c,d;
a=100*f_pi0mom;
b=10*f_dstps;
c=abs(a);
d=abs(b);

for(i=30;i<=c;i++){
for(j=20;j<=d;j++){

if(f_mygenflag == 1){NS[i][j]++;}else{NB[i][j]++;}

}}

}}}}}}}}}

//Float_t bestsig1,besti1,bestj1;
    }//  for(int i=0;i<nevt;i++)
Int_t count=-1;
Float_t bestsig1,besti1,bestj1;
for(j=25;j<=35;j++){
bestsig1=1000000; besti1=3; bestj1=20;
for(i=30;i<=100;i++){
sig[i][j]=NS[i][j]+NB[i][j];
sig[i][j]=sqrt(sig[i][j]);
sig[i][j]=sig[i][j]/NS[i][j];
count++;
Float_t ab,cd,ef;
ab=i/100;
cd=j/10;
ef=sig[i][j];
dt->SetPoint(count,ab,cd,ef);
//dt->SetPoint(count,ab,cd,ef);
if(sig[i][j]<=bestsig){bestsig=sig[i][j];besti=i;bestj=j;}
if(sig[i][j]<=bestsig1){bestsig1=sig[i][j];besti1=i;bestj1=j;}
}
//cout<<"best i = "<<besti1<<" best j = "<<bestj1<<endl;
}

//cout<<"for mass width = "<<f<<" :"<<endl;
cout<<"bestest i = "<<besti<<" bestest j = "<<bestj<<endl;
cout<<"count = "<<count<<endl;
//}

////////////////////////////////////////////////////////////////////////

gStyle->SetPalette(1);
//dt->Draw("surf1");
dt->Draw("colz0");

}


