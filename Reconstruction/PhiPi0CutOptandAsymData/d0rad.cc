//
#include "particle/Particle.h"
#include "particle/PID.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "kid/atc_pid.h"
#include "panther/panther.h"
#include "mdst/mdst.h"
#include "mdst/findKs.h"
#include "ip/IpProfile.h"
#include "belle.h"

#include "geninfo.h"
#include "d0rad.h"

//#include "d0rad_skim/d0rad.h"
#include "benergy/BeamEnergy.h"

// system include files

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include<sstream>

#include "basf/basfshm.h"
#include "basf/basfout.h"
#include "tuple/BelleTupleManager.h"
#include "panther/panther.h"

#include HEPEVT_H
#include BELLETDF_H
#include MDST_H
#include EVTCLS_H
#include "belleutil/debugout.h"

int A[100]={0},B[100]={0},n=0,C[100]={0},D[100]={0},nK=0,count1=0,multiplicity,pmultiplicity,nmultiplicity,eventnum=0,blueevent;

using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

 
  extern "C" Module_descr *mdcl_d0rad()
  {
    d0rad *module = new d0rad;
    Module_descr *dscr = new Module_descr ( "d0rad", module );

    BeamEnergy::define_global(dscr);

    return dscr;
  }

  
  
  d0rad::d0rad( void ) {
    
    //  std::strcpy( m_SkimFileName,   "d0rad.index" );
    
  }
  
  void d0rad::disp_stat( const char* ) {
}



void d0rad::end_run( BelleEvent*, int* ) {

}

void d0rad::other( int*, BelleEvent*, int* ) {
}

//Particle Types
  const Ptype d0rad::m_ptypeD0("D0");
  const Ptype d0rad::m_ptypeD0B("D0B");
  const Ptype d0rad::m_ptypeDstarP("D*+");
  const Ptype d0rad::m_ptypeDstarM("D*-");


//Cuts
const double DeltaMLowCut=0.13;//*
const double DeltaMHighCut=0.165;//*
const double DrCut=1.0;//*
const double DzCut=3.0;//*
const double PhiMassLowCut=0.9;//*
const double PhiMassHighCut=1.06;//*
const double PiZeroMassWidthCut=0.025;//*
const double DStarCMMomentumCut=2.0; //2.9*;
const double PiZeroMomentumCut=0.3; //0.75*;
//int n;


void d0rad::init ( int * ) {
  //  extern BasfOutputManager* BASF_Output;
  //  m_SkimFile = BASF_Output->open ( m_SkimFileName );

}

void d0rad::begin_run( BelleEvent*, int* ) {
  
  IpProfile::begin_run();
  BeamEnergy::begin_run();
  //  m_SkimFile->write();

}

void d0rad::term( void ) {

  // delete m_SkimFile;
char c1[15],c2[15];
string s1,s2,s3,s4;
s1="PID/";
s3=".txt";

int filenum;
ifstream fin;
fin.open("PID/num.txt");
fin>>s2;
fin.close();

strcpy(c2, s2.c_str());
filenum=atoi(c2);
filenum++;
stringstream ss;
ss << filenum;
s2 = ss.str();
s4=s1+s2+s3;
strcpy(c1, s4.c_str());

ofstream fout;
fout.open("PID/num.txt");
fout<<filenum;
fout.close();

fout.open(c1);
for(int i=1;i<=n;i++)
{
fout<<A[i]<<"\t"<<B[i]<<endl;
}
fout.close();

}

void d0rad::hist_def(void)
{
  // no histograms

  extern BelleTupleManager *BASF_Histogram;

  std::string title, title2;
  title = "information on D candidates";
  title2 = "information on D candidates2";

nt_check = BASF_Histogram->ntuple( title.c_str(),"Event BlueEvent Multiplicity PMultiplicity NMultiplicity PiSnSVD  PiSCT KpPID KmPID PisPID PiZMassNmc KMinHelicity KPlusHelicity ThetaHel PiZMomentum PisMom PisCosT Photon1Energy Photon2Energy Phot1Theta Phot2Theta PhiMass DZeroMass DStarCharge DStarMass DeltaM DSTID  D0ID  PiSID PhiID Pi0ID KpID KmID DSTF D0F  PiSF PhiF Pi0F KpF KmF MyGenFlag GenFlag2  DSTPS DSTCTS KpDR KpDZ KmDR KmDZ PisDR PisDZ", 1);

nt_mult =  BASF_Histogram->ntuple( title2.c_str(),"Multiplicity PMultiplicity NMultiplicity",2);


//  nt_d0kpkm = BASF_Histogram->ntuple( title2.c_str(),"KPlusHelicity KMinusMass KID PiZMassNmc PiZMomentum PhotonEnergy PhiMass DZeroMass DStarCharge DStarMass DeltaM DSTID DSTID1 DSTID2 DSTF DSTPS DR DZ", 2);

}

void d0rad::event( BelleEvent *, int *status ) {

  *status = 0;

        //To store beam info
        double E_LER = BeamEnergy::E_LER();          // e+ beam Energy 
        double E_HER = BeamEnergy::E_HER();          // e- beam energy
        double theta = BeamEnergy::Cross_angle();        // crossing angle

        double e,mass,px,py,pz;
        e=sqrt(E_LER*E_HER);

//For FillIDHEP
//int A[100],B[100]; 

  // ---------------- particle reconstruction and selection 

  vector<Particle>   k_p, k_m, pi_p, pi_m;
  makeKPi(k_p, k_m, pi_p, pi_m, 0);

//I think vertex factor constraints are added here. Have to double check
  with_imp_cut(k_p);
  with_imp_cut(k_m);
  with_imp_cut(pi_p);
  with_imp_cut(pi_m);

  withKaonId(k_p,0.6,3,1,5);
  withKaonId(k_m,0.6,3,1,5);

//--------------------------

  vector<Particle>  pi0;
  makePi0(pi0);

  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
    if(i->mdstPi0().gamma(0).ecl().energy()<0.05||
       i->mdstPi0().gamma(1).ecl().energy()<0.05||
       fabs(i->mdstPi0().mass()-.135)>PiZeroMassWidthCut || i->ptot()<PiZeroMomentumCut) {
//       fabs(i->mdstPi0().mass()-.135)>0.025 || i->ptot()<0.3) {
      pi0.erase(i); 
      --i;
    }

//  withPCut( pi0, 0.75);

//--------------------------  



  vector<Particle> k_s;
  Mdst_vee2_Manager &m = Mdst_vee2_Manager::get_manager();
  for (Mdst_vee2_Manager::iterator it = m.begin(); it != m.end(); it++) {
    Mdst_vee2* k0s = &(*it);
    if(k0s->kind() == 1)  // kind == 1 are K0s                                                                                                             
      {
	FindKs ffangks;
	
	ffangks.candidates(*k0s, IpProfile::position());
	
        Particle Kshort(*it);
	
	if (ffangks.goodKs()==1)
          k_s.push_back(Kshort);
      }
  }

  vector<Particle> gamma;
  vector<Particle> gamma_low;

  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin(); i != gamma_mag.end(); ++i){
    Particle gammaP(*i);
    
    if(gammaP.e()<0.050)
      continue;
    
    gamma_low.push_back(gammaP);

    if(gammaP.e()>0.200)
      gamma.push_back(gammaP);
  }

//-------------------------- 


  vector<Particle> phi;
  combination( phi, Ptype(333),k_p,k_m);
  withMassCut( phi, PhiMassLowCut, PhiMassHighCut);
//  withMassCut( phi, 0.9, 1.06);
  
//-------------------------- 

  vector<Particle> eta;
  combination( eta, Ptype(221),gamma_low,gamma_low);
  withMassCut( eta, 0.45, 0.65);
  withPCut( eta, 0.3);

  vector<Particle> D0;
  vector<Particle> D0K;

  // D0 -> gamma gamma
  //  combination(D0, m_ptypeD0, gamma, gamma , 0.2);
  // D0 -> phi pi0
    combination(D0, m_ptypeD0, phi, pi0 , 0.2);
  // D0 -> phi eta
  // combination(D0, m_ptypeD0, phi, eta, 0.2);
  // D0 -> phi gamma
  //  combination(D0, m_ptypeD0, phi, gamma, 0.2);
  // D0 -> K0s gamma
  //  combination(D0, m_ptypeD0, k_s, gamma, 0.2);
  // D0 -> K0s pi0
  //  combination(D0, m_ptypeD0, k_s, pi0, 0.2);
  // D0 -> pi0 pi0
  //combination(D0, m_ptypeD0, pi0, pi0, 0.2);
  // D0(K) -> K+ K-
    combination(D0K, m_ptypeD0, k_p, k_m, 0.2);

//----------------------------------------------

cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  NO  PROBLEMO  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

  // D*+ -> D0 pi-
  vector<Particle> DstP;
  vector<Particle> DstM;

//  vector<Particle> DstPK;
//  vector<Particle> DstMK;

  combination( DstP, m_ptypeDstarP, D0, pi_p);
  combination( DstM, m_ptypeDstarM, D0, pi_m);
  withMassDifCut( DstP, DeltaMLowCut, DeltaMHighCut, 0);
  withMassDifCut( DstM, DeltaMLowCut, DeltaMHighCut, 0);
  withPSCut( DstP, DStarCMMomentumCut);
  withPSCut( DstM, DStarCMMomentumCut);


    Belle_event_Manager& EvtMgr=Belle_event_Manager::get_manager();
    Belle_event& Evt = *EvtMgr.begin();
    int ExpNo = Evt.ExpNo();
    int RunNo = Evt.RunNo();
    int EvtNo = Evt.EvtNo();

eventnum++;
multiplicity= 0;
pmultiplicity= 0;
nmultiplicity= 0;
blueevent=0;//for correct signal events(blue colour in my plots)

for (int i=0;i<DstP.size();i++) {
AddMult(DstP[i]);
}
for (int i=0;i<DstM.size();i++) {
AddMult(DstM[i]);
}

    nt_mult->column("Multiplicity",multiplicity);
    nt_mult->column("PMultiplicity",pmultiplicity);
    nt_mult->column("NMultiplicity",nmultiplicity);
    nt_mult->dumpData();

  for (int i=0;i<DstP.size();i++) {
FillTuple(DstP[i],nt_check,A,B,n);
}
  for (int i=0;i<DstM.size();i++) {
FillTuple(DstM[i],nt_check,A,B,n);
}

    *status = 1;
  
}

void d0rad::AddMult(Particle& Dst_cand){
HepLorentzVector DSTPPS;
double delm,phimass,pi0mass,pi0mom,dstps,d0mass;
//f(f_delm > 0.14 && f_delm < 0.16) {
//if(f_phimass > 1.01 && f_phimass < 1.03) {
//if(f_pi0mass > 0.119 && f_pi0mass < 0.151) {
//if(f_pi0mom > 0.55) {//0.75
//if(f_dstps > 2.5) {//2.9
//if(f_d0mass > 1.65 && f_d0mass < 2.05) {
delm=Dst_cand.mass()-Dst_cand.child(0).mass();
phimass=Dst_cand.child(0).child(0).mass();
pi0mass=Dst_cand.child(0).child(1).mass();
d0mass=Dst_cand.child(0).mass();
pi0mom=Dst_cand.child(0).child(1).p().vect().mag();
DSTPPS=pStar(Dst_cand);
dstps=DSTPPS.vect().mag();
if(delm > 0.14 && delm < 0.16) {
if(phimass > 1.01 && phimass < 1.03) {
if(pi0mass > 0.119 && pi0mass < 0.151) {
if(pi0mom > 0.55) {//0.75
if(dstps > 2.5) {//2.9
//if(d0mass > 1.65 && d0mass < 2.05) {
if(d0mass > 1.83 && d0mass < 1.89) {

multiplicity++;
if(Dst_cand.charge() == 1){pmultiplicity++;}
if(Dst_cand.charge() == -1){nmultiplicity++;}

if(getMCtruthFlag(Dst_cand)==1 && fabs(IDhep(Dst_cand))==413){blueevent=1;}

}}}}}}

}

void d0rad::FillTuple(Particle& Dst_cand, BelleTuple *nt, int A[100], int B[100], int & n){
HepLorentzVector DSTPPS,pizero4V,gamma4V,dzero4V,phi4V,kp4V,km4V,pis4v;
Hep3Vector pizero_boost,phi_boost;

double massdiff,kid,dstpps,dstcts,cos_helicity,cos_helicity2,cos_helicity3,pis3vmag,pis3vcost;

//D0ID PiSID PhiID Pi0ID KpID KmID

      setMCtruth(Dst_cand);
    nt->column("Event",eventnum);
    nt->column("BlueEvent",blueevent);
    nt->column("Multiplicity",multiplicity);
    nt->column("PMultiplicity",pmultiplicity);
    nt->column("NMultiplicity",nmultiplicity);


    nt->column("DStarMass",Dst_cand.mass());
    nt->column("DStarCharge",Dst_cand.charge());
      nt->column("DSTID",IDhep(Dst_cand));
      nt->column("D0ID",IDhep(Dst_cand.child(0)));
      nt->column("PiSID",IDhep(Dst_cand.child(1)));
      nt->column("PhiID",IDhep(Dst_cand.child(0).child(0)));
      nt->column("Pi0ID",IDhep(Dst_cand.child(0).child(1)));
      nt->column("KpID",IDhep(Dst_cand.child(0).child(0).child(0)));
      nt->column("KmID",IDhep(Dst_cand.child(0).child(0).child(1)));
      nt->column("DSTF",getMCtruthFlag(Dst_cand));
      nt->column("D0F",getMCtruthFlag(Dst_cand.child(0)));
      nt->column("PiSF",getMCtruthFlag(Dst_cand.child(1)));
      nt->column("PhiF",getMCtruthFlag(Dst_cand.child(0).child(0)));
      nt->column("Pi0F",getMCtruthFlag(Dst_cand.child(0).child(1)));
      nt->column("KpF",getMCtruthFlag(Dst_cand.child(0).child(0).child(0)));
      nt->column("KmF",getMCtruthFlag(Dst_cand.child(0).child(0).child(1)));

      nt->column("PiSnSVD",(Dst_cand.child(1).mdstCharged().trk().mhyp(2).nhits(3)+Dst_cand.child(1).mdstCharged().trk().mhyp(2).nhits(4)));//SVD hit by soft Pi



int ID=IDhep(Dst_cand); 
FillIDHEP(ID,A,B,n);

//CMF Momentum of DSTP
DSTPPS=pStar(Dst_cand);
dstpps=DSTPPS.vect().mag();
      nt->column("DSTPS",dstpps);
dstcts=DSTPPS.vect().cosTheta();
      nt->column("DSTCTS",dstcts);

//Pi slow Cos Theta
      nt->column("PiSCT",Dst_cand.child(1).p().vect().cosTheta());


pis4v=Dst_cand.child(1).p();
pis3vmag=pis4v.vect().mag();
      nt->column("PisMom",pis3vmag);
pis3vcost=pis4v.vect().cosTheta();
      nt->column("PisCosT",pis3vcost);

//helicity dist. of K+
dzero4V=Dst_cand.child(0).p();
phi4V=Dst_cand.child(0).child(0).p();
kp4V=Dst_cand.child(0).child(0).child(0).p();
phi_boost=-1*phi4V.boostVector();
dzero4V.boost(phi_boost);
kp4V.boost(phi_boost);
cos_helicity2=cos(dzero4V.angle(kp4V));

//helicity dist. of K-
km4V=Dst_cand.child(0).child(0).child(1).p();
km4V.boost(phi_boost);
cos_helicity3=cos(dzero4V.angle(km4V));
double thetahel=dzero4V.angle(kp4V);

      nt->column("KPlusHelicity",cos_helicity2);
      nt->column("KMinHelicity",cos_helicity3);
      nt->column("ThetaHel",thetahel);

//Particle ID
massdiff=Dst_cand.mass()-Dst_cand.child(0).mass();
    nt->column("DeltaM",massdiff);

kid=kaonId(Dst_cand.child(0).child(0).child(0));
    nt->column("KpPID",kid);

kid=kaonId(Dst_cand.child(0).child(0).child(1));
    nt->column("KmPID",kid);

kid=kaonId(Dst_cand.child(1));
    nt->column("PisPID",kid);

//Photon thetas
double theta1, theta2, px1, py1, pz1, px2, py2, pz2;
px1=Dst_cand.child(0).child(1).child(0).px();
py1=Dst_cand.child(0).child(1).child(0).py();
pz1=Dst_cand.child(0).child(1).child(0).pz();

px2=Dst_cand.child(0).child(1).child(1).px();
py2=Dst_cand.child(0).child(1).child(1).py();
pz2=Dst_cand.child(0).child(1).child(1).pz();

theta1=pz1/sqrt(px1*px1+py1*py1);
theta1=(atan(theta1)*180)/3.14159265359;
theta2=pz2/sqrt(px2*px2+py2*py2);
theta2=(atan(theta2)*180)/3.14159265359;



      nt->column("Photon1Energy",Dst_cand.child(0).child(1).child(0).e());
      nt->column("Photon2Energy",Dst_cand.child(0).child(1).child(1).e());
      nt->column("Phot1Theta",theta1);
      nt->column("Phot2Theta",theta2);

    nt->column("PiZMomentum",Dst_cand.child(0).child(1).ptot());
    nt->column("PiZMassNmc",Dst_cand.child(0).child(1).mdstPi0().mass());
    nt->column("DZeroMass",Dst_cand.child(0).mass());
    nt->column("PhiMass",Dst_cand.child(0).child(0).mass());
//DR and DZ
        HepVector b(param_at_ip(Dst_cand.child(0).child(0).child(0)));
        HepVector b1(param_at_ip(Dst_cand.child(0).child(0).child(1)));
        HepVector b2(param_at_ip(Dst_cand.child(1)));
        double dr,dz;

    dr=b[0];  dz=b[3];
    nt->column("KpDR",dr);
    nt->column("KpDZ",dz);

    dr=b1[0];  dz=b1[3];
    nt->column("KmDR",dr);
    nt->column("KmDZ",dz);

    dr=b2[0];  dz=b2[3];
    nt->column("PisDR",dr);
    nt->column("PisDZ",dz);







int Temp[10]={0};

//My Gen Logic
int mygenflag=0,genflag2=0,allflags=0;
if(getMCtruthFlag(Dst_cand)==1 && getMCtruthFlag(Dst_cand.child(0))==1 && getMCtruthFlag(Dst_cand.child(1))==1 && getMCtruthFlag(Dst_cand.child(0).child(1))==1 && getMCtruthFlag(Dst_cand.child(0).child(0))==1 && getMCtruthFlag(Dst_cand.child(0).child(0).child(1))==1 && getMCtruthFlag(Dst_cand.child(0).child(0).child(0))==1){allflags=1;}

Gen_hepevt KpGen,KmGen,PhiGen,PizGen,DzGen,PisGen,DstGen,Gam1Gen,Gam2Gen;
//Setting possible links to all the different reconstructed candidates
KpGen= Dst_cand.child(0).child(0).child(0).genHepevt();
KmGen= Dst_cand.child(0).child(0).child(1).genHepevt();
PhiGen= Dst_cand.child(0).child(0).genHepevt();
PizGen= Dst_cand.child(0).child(1).genHepevt();
DzGen= Dst_cand.child(0).genHepevt();
PisGen= Dst_cand.child(1).genHepevt();
DstGen= Dst_cand.genHepevt();
Gam1Gen= Dst_cand.child(0).child(1).child(0).genHepevt();
Gam2Gen= Dst_cand.child(0).child(1).child(1).genHepevt();


if(KpGen && KmGen && PhiGen && PizGen && DzGen && PisGen && DstGen && Gam1Gen && Gam2Gen){//if all links exist

Temp[0]=DstGen.idhep();
Temp[1]=PisGen.idhep();
Temp[2]=DzGen.idhep();
Temp[3]=PizGen.idhep();
Temp[4]=PhiGen.idhep();
Temp[5]=KmGen.idhep();
Temp[6]=KpGen.idhep();

if(KpGen.idhep()==321 && KmGen.idhep()==-321 && PhiGen.idhep()==333 && PizGen.idhep()==111){//Are all the particles correctly identified?
if( (PisGen.idhep()==211 && DstGen.idhep()==413 && DzGen.idhep()==421 ) || (PisGen.idhep()==-211 && DstGen.idhep()==-413 && DzGen.idhep()==-421 )  ){
genflag2=1;

//No extra FSP
if((DstGen.daLast()-DstGen.daFirst())==1 && (DzGen.daLast()-DzGen.daFirst())==1 && (PhiGen.daLast()-PhiGen.daFirst())==1 && (PizGen.daLast()-PizGen.daFirst())==1){

//Check all mother daughter relationships
if(KpGen.mother().get_ID()== PhiGen.get_ID() && KmGen.mother().get_ID()== PhiGen.get_ID()){
if(PhiGen.mother().get_ID()== DzGen.get_ID() && PizGen.mother().get_ID()== DzGen.get_ID()){
if(DzGen.mother().get_ID()==DstGen.get_ID() && PisGen.mother().get_ID()==DstGen.get_ID() ){
if(Gam1Gen.mother().get_ID()==PizGen.get_ID() && Gam2Gen.mother().get_ID()==PizGen.get_ID() ){
//cout<<"+-+-+-+-+-+-+-+-+-+-+-+  SUCCESS!!!  +-+-+-+-+-+-+--+-+-+-+-+-+-+-"<<endl;
mygenflag=1;
}}}}

}//if((DstGen.daLast()-DstGen.daFirst())==1...
}//Correct idheps
}//Are all the particles correctly identified?


}//if all links exist

if(mygenflag==1 && allflags!=1){PrintGenTable(Dst_cand,8,count1);}
if(mygenflag!=1 && allflags==1){PrintGenTable(Dst_cand,9,count1);}

    nt->column("MyGenFlag",mygenflag);
    nt->column("GenFlag2",genflag2);

    nt->dumpData();
}//  void d0rad::FillTuple(Particle& Dst_cand, BelleTuple *nt)

void d0rad::PrintGenTable(Particle& Dst_cand,int GenTableFlag, int & counts)
{
int ref=GenTableFlag;
//Genstuff
int maxdau=1;
ofstream fout1;
int a1,a2,a3,a4,a5,a6,a7;
a1=Dst_cand.genHepevt().get_ID();
a2=Dst_cand.child(0).genHepevt().get_ID();
a3=Dst_cand.child(1).genHepevt().get_ID();
a4=Dst_cand.child(0).child(0).genHepevt().get_ID();
a5=Dst_cand.child(0).child(1).genHepevt().get_ID();
a6=Dst_cand.child(0).child(0).child(0).genHepevt().get_ID();
a7=Dst_cand.child(0).child(0).child(1).genHepevt().get_ID();


switch (ref)
{
case 1:  fout1.open("GenTables1.txt", ios_base::app);
    break;
case 2:  fout1.open("GenTables2.txt", ios_base::app);
    break;
case 3:  fout1.open("GenTables3.txt", ios_base::app);
    break;
case 4:  fout1.open("GenTables4.txt", ios_base::app);
    break;
case 5:  fout1.open("GenTables5.txt", ios_base::app);
    break;
case 6:  fout1.open("GenTables6.txt", ios_base::app);
    break;
case 7:  fout1.open("GenTablesTest1.txt", ios_base::app);
    break;
case 8:  fout1.open("GenTablesTest2.txt", ios_base::app);
    break;
case 9:  fout1.open("GenTablesTest3.txt", ios_base::app);
    break;
default: break;
}

counts++;
fout1<<"start\tX\tX\tX\tX\tX\tX\tX"<<endl;
//fout1<<"----------------------------------------------------------------------------"<<endl;
        Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
fout1<<"$\tID\tIDHEP\tDa1\tDa2\tMo1\tMo2\tPartFlag"<<endl;

int count,ParticipantFlag;//, flag;

for(count=1;count<=500;count++){
ParticipantFlag=0;

if(count>maxdau){break;}
                const Gen_hepevt& Guy = GenMgr(Panther_ID(count));
if(count==a1 || count==a2 || count==a3 || count==a4 || count==a5 || count==a6 || count==a7){ParticipantFlag=1;}

fout1<<"@\t"<<Guy.get_ID()<<"\t"<<Guy.idhep()<<"\t"<<Guy.daFirst()<<"\t"<<Guy.daLast()<<"\t"<<Guy.moFirst()<<"\t"<<Guy.moLast()<<"\t"<<ParticipantFlag<<endl;
if(Guy.daLast()>maxdau){maxdau=Guy.daLast();}

}
fout1<<"end\tX\tX\tX\tX\tX\tX\tX"<<endl;
fout1.close();
}



void d0rad::FillIDHEP(int a, int A[100], int B[100], int & n)
{
int i,st=0;
for(i=1;i<=n;i++){
if(a==A[i]){
st=1;
B[i]++;}
}

if(st==0){
n++;
A[n]=a;
B[n]++;
}

}//void d0rad::FillIDHEP(int a, int A[100], int B[100], int & n)


  void d0rad::with_imp_cut(std::vector<Particle> &list) {
    for(int i=0;i<(int)list.size();++i){
      if(list[i].mdstCharged()){
	HepVector a(param_at_ip(list[i]));
	if (!(abs(a[0])<DrCut && abs(a[3])<DzCut)){
	  list.erase(list.begin()+i);
	  --i;
	}
      }
    }
  }


  HepVector d0rad::param_at_ip(Particle &p){

    const Mdst_charged charged(p.mdstCharged());

    double thisMass = p.mass();

    int hyp = 4;
    if(thisMass < 0.005){ // e = 0.000511                                                                                         
      hyp = 0;
    }else if(thisMass < 0.110){ // mu = 0.1056                                                                                               
      hyp = 1;
    }else if(thisMass < 0.200){ // pi = 0.13956                                                                                                 
      hyp = 2;
    }else if(thisMass < 0.5){ // K = 0.4936                                                                                                          
      hyp = 3;
    }
    const HepPoint3D pivot(charged.trk().mhyp(hyp).pivot_x(),
			   charged.trk().mhyp(hyp).pivot_y(),
			   charged.trk().mhyp(hyp).pivot_z());

    HepVector  a(5);
    a[0] = charged.trk().mhyp(hyp).helix(0);
    a[1] = charged.trk().mhyp(hyp).helix(1);
    a[2] = charged.trk().mhyp(hyp).helix(2);
    a[3] = charged.trk().mhyp(hyp).helix(3);
    a[4] = charged.trk().mhyp(hyp).helix(4);
    HepSymMatrix Ea(5,0);
    Ea[0][0] = charged.trk().mhyp(hyp).error(0);
    Ea[1][0] = charged.trk().mhyp(hyp).error(1);
    Ea[1][1] = charged.trk().mhyp(hyp).error(2);
    Ea[2][0] = charged.trk().mhyp(hyp).error(3);
    Ea[2][1] = charged.trk().mhyp(hyp).error(4);
    Ea[2][2] = charged.trk().mhyp(hyp).error(5);
    Ea[3][0] = charged.trk().mhyp(hyp).error(6);
    Ea[3][1] = charged.trk().mhyp(hyp).error(7);
    Ea[3][2] = charged.trk().mhyp(hyp).error(8);
    Ea[3][3] = charged.trk().mhyp(hyp).error(9);
    Ea[4][0] = charged.trk().mhyp(hyp).error(10);
    Ea[4][1] = charged.trk().mhyp(hyp).error(11);
    Ea[4][2] = charged.trk().mhyp(hyp).error(12);
    Ea[4][3] = charged.trk().mhyp(hyp).error(13);
    Ea[4][4] = charged.trk().mhyp(hyp).error(14);
    Helix helix(pivot, a, Ea);

    const Hep3Vector&   IP     = IpProfile::position();
    if (IP.mag())
      helix.pivot(IP);
    return helix.a();
  }


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
