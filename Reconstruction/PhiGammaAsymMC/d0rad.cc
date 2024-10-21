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
#include "eid/eid.h"

#include "pi0eta_prob.h"
#include "geninfo.h"
#include "d0rad.h"

//#include "d0rad_skim/d0rad.h"
#include "benergy/BeamEnergy.h"

// system include files

#include <cmath>
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

int A[100]={0},B[100]={0},n=0,C[100]={0},D[100]={0},nK=0,IdMatrix[10000][10]={0},nIM=0,count1=0,goodcandflag;

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
const double DeltaMLowCut=0.13;
const double DeltaMHighCut=0.165;
const double DrCut=1.0;
const double DzCut=3.0;
const double PhiMassLowCut=0.9;
const double PhiMassHighCut=1.06;
const double DStarCMMomentumCut=2.0;

//int n;


void d0rad::init ( int * ) {
  //  extern BasfOutputManager* BASF_Output;
  //  m_SkimFile = BASF_Output->open ( m_SkimFileName );

}

void d0rad::begin_run( BelleEvent*, int* ) {
  
  IpProfile::begin_run();
  BeamEnergy::begin_run();
        eid::init_data();
  //  m_SkimFile->write();
// Data or MC?
    Belle_runhead_Manager &runhead_m = Belle_runhead_Manager::get_manager();
    Belle_runhead &runhead = runhead_m( ( Panther_ID) 1 );

    if(runhead)
      {
        if(runhead.ExpMC() == 1)
          {
            dataType = 0;
            std::cout << ">>>User_ana[begin_run]: running on real data" << std::endl;
          }
        else
          {
            dataType = 1;
            std::cout << ">>>User_ana[begin_run]: running on Monte Carlo" << std::endl;
          }
      }

}

void d0rad::term( void ) {

  // delete m_SkimFile;

}

void d0rad::hist_def(void)
{
  // no histograms

  extern BelleTupleManager *BASF_Histogram;

  std::string title;
  title = "information on D candidates";

nt_check = BASF_Histogram->ntuple( title.c_str()," XVertex YVertex Zvertex ConvPhoM EnergyFlag PiSMom PiSTheta MyPGenFlag PiFlag EtaFlag KID VetoM1 VetoM1_2 VetoM2 VetoM2_2 Heli2 Heli3 KPlusHelicity PhiMass DZeroMass DStarCharge DStarMass DeltaM MGamEnergy RMGamEnergy GamEnergy MuConVX MuConVY MuConVZ EConVX EConVY EConVZ Ga_e9e25 DSTID D0ID PiSID PhiID KpID KmID GamID DSTF D0F  PiSF PhiF KpF KmF GamF MyGenFlag DSTPS DSTCTS DR DZ P_Piz P_Piz_2 P_Eta P_Eta_2 PiAngle EtaAngle CPiMom CEtaMom KmPID KpPID PisPID PhotonEnergy PhotTheta KpDR KpDZ KmDR KmDZ PisDR PisDZ", 1);


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
//       fabs(i->mdstPi0().mass()-.135)>PiZeroMassWidthCut || i->ptot()<PiZeroMomentumCut) {
       fabs(i->mdstPi0().mass()-.135)>0.025 || i->ptot()<0.3) {
      pi0.erase(i); 
      --i;
    }


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
//  vector<Particle> D0K;

  // D0 -> gamma gamma
  //  combination(D0, m_ptypeD0, gamma, gamma , 0.2);
  // D0 -> phi pi0
  // combination(D0, m_ptypeD0, phi, pi0 , 0.2);
  // D0 -> phi eta
  // combination(D0, m_ptypeD0, phi, eta, 0.2);
  // D0 -> phi gamma
    combination(D0, m_ptypeD0, phi, gamma, 0.2);
  // D0 -> K0s gamma
  //  combination(D0, m_ptypeD0, k_s, gamma, 0.2);
  // D0 -> K0s pi0
  //  combination(D0, m_ptypeD0, k_s, pi0, 0.2);
  // D0 -> pi0 pi0
  //combination(D0, m_ptypeD0, pi0, pi0, 0.2);
  // D0(K) -> K+ K-
  // combination(D0K, m_ptypeD0, k_p, k_m, 0.2);

//----------------------------------------------


// D*+ -> D0 pi-
  vector<Particle> DstP;
  vector<Particle> DstM;

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


  for (int i=0;i<DstP.size();i++) {
FillTuple(DstP[i],nt_check,A,B,n,IdMatrix,nIM);
}

  for (int i=0;i<DstM.size();i++) {
FillTuple(DstM[i],nt_check,A,B,n,IdMatrix,nIM);
}

    *status = 1;

  
}

void d0rad::FillTuple(Particle& Dst_cand, BelleTuple *nt, int A[100], int B[100], int & n, int IdMatrix[10000][10], int & nIM){
HepLorentzVector DSTPPS,pizero4V,gamma4V,dzero4V,phi4V,kp4V;
Hep3Vector pizero_boost,phi_boost, d0_boost;

goodcandflag=0;

double massdiff,kid,dstpps,dstcts,cos_helicity,cos_helicity2;

      nt->column("DStarMass",Dst_cand.mass());
      nt->column("DStarCharge",Dst_cand.charge());
      nt->column("GamEnergy",Dst_cand.child(0).child(1).e());



//if(getMCtruthFlag(Dst_cand)==1 && fabs(IDhep(Dst_cand))!=413){PrintGenTable();}

//int ID=IDhep(Dst_cand); 
//FillIDHEP(ID,A,B,n);

//CMF Momentum of DSTP
DSTPPS=pStar(Dst_cand);
dstpps=DSTPPS.vect().mag();
      nt->column("DSTPS",dstpps);
dstcts=DSTPPS.vect().cosTheta();
      nt->column("DSTCTS",dstcts);

      //e9/e25
      Mdst_ecl_aux_Manager &eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
      Mdst_ecl_aux &aux = eclaux_mgr( Panther_ID(Dst_cand.child(0).child(1).mdstGamma().ecl().get_ID()));
      nt->column("Ga_e9e25",aux.e9oe25());

//Koppenburg Veto
double ppiz,peta,mm1,mm2;
Mdst_gamma OUT1,OUT2;
ppiz=Closest_Pi0_Probability(Dst_cand.child(0).child(1).mdstGamma(),OUT1,mm1);
peta=Closest_Eta_Probability(Dst_cand.child(0).child(1).mdstGamma(),OUT2,mm2);

      nt->column("P_Piz",ppiz);
      nt->column("P_Eta",peta);
      nt->column("VetoM1",mm1);
      nt->column("VetoM2",mm2);

//Extended Koppenburg Veto
/*//Step 1: Getting the electrons
        std::vector<Particle> elec_list;
        elec_list.clear();
        Mdst_charged_Manager &chg_mgr = Mdst_charged_Manager::get_manager();
        for( std::vector<Mdst_charged>::const_iterator i = chg_mgr.begin();
                        i!= chg_mgr.end(); i++){
                const Mdst_charged &chg = *i;
                eid sel_e(chg);
                double prob_e = sel_e.prob( 3, -1, 5 );
//                if(prob_e >= 0.9){    Anze's skim value
                if(prob_e >= 0.01){    //Making it Nice and Loose
                        Particle elec(chg, Ptype(chg.charge()>0 ? "E+" : "E-") );
                        elec_list.push_back( elec );
                }
			}

  vector<Particle>   e_p, e_m;
  e_p.clear(); e_m.clear();
  vector<Particle>   conv_gamma;
  conv_gamma.clear();

for (int i=0;i<elec_list.size();i++) {//e+
if(elec_list[i].charge()==1){e_p.push_back(elec_list[i]);}
else if(elec_list[i].charge()==-1){e_m.push_back(elec_list[i]);}
}

  combination(conv_gamma, Ptype(22),e_p,e_m);
*/
//Step 1: Converted gamma list
  vector<Particle> conv_gamma;
  conv_gamma.clear();
  Mdst_vee2_Manager &m = Mdst_vee2_Manager::get_manager();
  for (Mdst_vee2_Manager::iterator it = m.begin(); it != m.end(); it++) {
    Mdst_vee2* congam = &(*it);
    if(congam->kind() == 4)  // kind == 4 are converted gammas                                                                                  
      {
        Particle ConvGamm(*it);
if(ConvGamm.mass() < 0.1) {conv_gamma.push_back(ConvGamm);}
             
      }
  }
                                                                    
//Step 2: Getting the full gamma list
  vector<Particle>   Full_Gamma_List;
  Full_Gamma_List.clear();
  Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();
  for (std::vector<Mdst_gamma>::iterator it1= GammaMgr.begin(); it1 !=GammaMgr.end(); ++it1){
Full_Gamma_List.push_back(*it1);
}

Full_Gamma_List.insert(Full_Gamma_List.end(), conv_gamma.begin(), conv_gamma.end()); //append converted photons

//Step 3: Applying the Veto
double convertedphotonmass=0,xvertex=0,yvertex=0,zvertex=0;
  HepLorentzVector Reco1( Dst_cand.child(0).child(1).mdstGamma().px(), Dst_cand.child(0).child(1).mdstGamma().py(), Dst_cand.child(0).child(1).mdstGamma().pz() );
  Reco1= HepLorentzVector( Reco1.px(), Reco1.py(), Reco1.pz(), Reco1.rho()) ;
  double bp = -100., mm = -1;//Pi0
  double bp_e = -100., mm_e = -1;//Eta

for(int i=0;i<Full_Gamma_List.size();i++){

    HepLorentzVector PG( Full_Gamma_List[i].px(), Full_Gamma_List[i].py(), Full_Gamma_List[i].pz()  ) ;
    PG = HepLorentzVector( PG.px(), PG.py(), PG.pz(), PG.rho()) ;
    double m = (Reco1 + PG).m();
    double p = Pi0_Eta_Prob(1, m, PG.e(), PG.theta() );//Pi0 Probability
    if (p>bp) {
      bp = p ;
      mm = m ;
if(Full_Gamma_List[i].mdstVee2()){convertedphotonmass=Full_Gamma_List[i].mass();
xvertex=Full_Gamma_List[i].mdstVee2().vx();
yvertex=Full_Gamma_List[i].mdstVee2().vy();
zvertex=Full_Gamma_List[i].mdstVee2().vz();
}
    }
    double p_e = Pi0_Eta_Prob(2, m, PG.e(), PG.theta() );//Eta Probability
    if (p_e>bp_e) {
      bp_e = p_e ;
      mm_e = m ;
    }
  }
nt->column("ConvPhoM",convertedphotonmass);
nt->column("XVertex",xvertex);
nt->column("YVertex",yvertex);
nt->column("ZVertex",zvertex);




      nt->column("P_Piz_2",bp);
      nt->column("VetoM1_2",mm);
      nt->column("P_Eta_2",bp_e);
      nt->column("VetoM2_2",mm_e);

//helicity dist. of K+
dzero4V=Dst_cand.child(0).p();
phi4V=Dst_cand.child(0).child(0).p();
kp4V=Dst_cand.child(0).child(0).child(0).p();
phi_boost=-1*phi4V.boostVector();
d0_boost=-1*dzero4V.boostVector();

dzero4V.boost(phi_boost);
kp4V.boost(phi_boost);
cos_helicity2=cos(dzero4V.angle(kp4V));

    nt->column("KPlusHelicity",cos_helicity2);
/*
cout<<"unboosted phi 4-mom is"<<phi4V.e()<<"\t"<<phi4V.px()<<"\t"<<phi4V.py()<<"\t"<<phi4V.pz()<<endl;
phi4V = phi4V+phi4V;
phi4V.boost(phi_boost);
cout<<"boosted phi 4-mom is"<<phi4V.e()<<"\t"<<phi4V.px()<<"\t"<<phi4V.py()<<"\t"<<phi4V.pz()<<endl<<endl;
    nt->column("KPlusHelicity",cos_helicity2);
*/

//Alternate helicity defn
cos_helicity2=cos(phi4V.angle(kp4V));
    nt->column("Heli2",cos_helicity2);

//Alternate helicity defn 3
phi4V.boost(d0_boost);
cos_helicity2=cos(phi4V.angle(kp4V));
    nt->column("Heli3",cos_helicity2);

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
px1=Dst_cand.child(0).child(1).px();
py1=Dst_cand.child(0).child(1).py();
pz1=Dst_cand.child(0).child(1).pz();


theta1=sqrt(px1*px1+py1*py1)/pz1;
theta1=atan(theta1);

      nt->column("PhotonEnergy",Dst_cand.child(0).child(1).e());
      nt->column("PhotTheta",theta1);

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





massdiff=Dst_cand.mass()-Dst_cand.child(0).mass();
    nt->column("DeltaM",massdiff);

   kid=kaonId(Dst_cand.child(0).child(0).child(0));
    nt->column("KID",kid);

    nt->column("DZeroMass",Dst_cand.child(0).mass());
     nt->column("PhiMass",Dst_cand.child(0).child(0).mass());

//DR and DZ(They are taken only for K+ candidates)
//        HepVector b(param_at_ip(Dst_cand.child(0).child(0).child(0)));
//        double dr,dz;
        dr=b[0];
        dz=b[3];

    nt->column("DR",dr);
    nt->column("DZ",dz);

    nt->column("PiSMom",Dst_cand.child(1).p().vect().mag());
    nt->column("PiSTheta",Dst_cand.child(1).p().vect().theta());


int Temp[10]={0};

if(dataType){
      setMCtruth(Dst_cand);
      nt->column("DSTID",IDhep(Dst_cand));
      nt->column("D0ID",IDhep(Dst_cand.child(0)));
      nt->column("PiSID",IDhep(Dst_cand.child(1)));
      nt->column("PhiID",IDhep(Dst_cand.child(0).child(0)));
      nt->column("KpID",IDhep(Dst_cand.child(0).child(0).child(0)));
      nt->column("KmID",IDhep(Dst_cand.child(0).child(0).child(1)));
      nt->column("GamID",IDhep(Dst_cand.child(0).child(1)));
      nt->column("DSTF",getMCtruthFlag(Dst_cand));
      nt->column("D0F",getMCtruthFlag(Dst_cand.child(0)));
      nt->column("PiSF",getMCtruthFlag(Dst_cand.child(1)));
      nt->column("PhiF",getMCtruthFlag(Dst_cand.child(0).child(0)));
      nt->column("KpF",getMCtruthFlag(Dst_cand.child(0).child(0).child(0)));
      nt->column("KmF",getMCtruthFlag(Dst_cand.child(0).child(0).child(1)));
      nt->column("GamF",getMCtruthFlag(Dst_cand.child(0).child(1) ));

//My Gen Logic
int mygenflag=0,mypgenflag=0,genflag2=0,allflags=0,piflag=0,etaflag=0;
//allflags

Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

Gen_hepevt KpGen,KmGen,PhiGen,GamGen,DzGen,PisGen,DstGen;
//Setting possible links to all the different reconstructed candidates
KpGen= Dst_cand.child(0).child(0).child(0).genHepevt();
KmGen= Dst_cand.child(0).child(0).child(1).genHepevt(); 
PhiGen= Dst_cand.child(0).child(0).genHepevt();
GamGen= Dst_cand.child(0).child(1).genHepevt();
DzGen= Dst_cand.child(0).genHepevt();
PisGen= Dst_cand.child(1).genHepevt();
DstGen= Dst_cand.genHepevt();


//logic to get info from pi0 and eta contamination
int energyflag=0;
double cpimom=0,piangle=0,mgamenergy=0,rmgamenergy=0,econvx=0,econvy=0,econvz=0,muconvx=0,muconvy=0,muconvz=0,cetamom=0,etaangle=0;

if(GamGen){
if(GamGen.idhep()==22){
if(GamGen.mother().idhep()==111 && ( GamGen.mother().daLast()-GamGen.mother().daFirst() )==1){piflag=1;}
if(GamGen.mother().idhep()==221 && ( GamGen.mother().daLast()-GamGen.mother().daFirst() )==1){etaflag=1;}
if(piflag==1){
double px1,py1,pz1,p1,e1;
px1=GamGen.mother().PX();
py1=GamGen.mother().PY();
pz1=GamGen.mother().PZ();
p1=sqrt((px1*px1)+(py1*py1)+(pz1*pz1));
//nt->column("CPiMom",p1);//Contaminating Pi***********************
cpimom=p1;//*******************************************************
//opening angle
Gen_hepevt& child1 = GenMgr(Panther_ID(GamGen.mother().daFirst()));
Gen_hepevt& child2 = GenMgr(Panther_ID(GamGen.mother().daLast()));
Hep3Vector A1(child1.PX(),child1.PY(),child1.PZ()); Hep3Vector B1(child2.PX(),child2.PY(),child2.PZ());
piangle=A1.angle(B1);
//nt->column("PiAngle",piangle);//***********************************

//nt->column("MGamEnergy",child1.E()+child2.E()-Dst_cand.child(0).child(1).e());//Missing gamma energy//*******************
mgamenergy=child1.E()+child2.E()-Dst_cand.child(0).child(1).e();

if(GamGen.get_ID()==child1.get_ID()){//child2 is missing
//nt->column("RMGamEnergy",child2.E());//Real Missing gamma energy//********************************
rmgamenergy=child2.E();
if(child2.E() > child1.E()){energyflag=1;}//More energetic photon missing
else{energyflag=-1;}//Less energetic photon missing
//nt->column("EnergyFlag",energyflag);********************** **************************
//Converted Photons
if(child2.daLast()-child2.daFirst()==1)
{
Gen_hepevt& grandchild1 = GenMgr(Panther_ID(child2.daFirst()));
Gen_hepevt& grandchild2 = GenMgr(Panther_ID(child2.daLast()));
if(fabs(grandchild1.idhep())==11 && fabs(grandchild2.idhep())==11){
econvx=grandchild1.VX();
econvy=grandchild1.VY();
econvz=grandchild1.VZ();
//nt->column("EConVX",grandchild1.VX());**************************************
//nt->column("EConVY",grandchild1.VY());****************************************8
//nt->column("EConVZ",grandchild1.VZ());*******************************************
}
if(fabs(grandchild1.idhep())==13 && fabs(grandchild2.idhep())==13){
muconvx=grandchild1.VX();
muconvy=grandchild1.VY();
muconvz=grandchild1.VZ();
}
}

}
else{//child1 is missing
//nt->column("RMGamEnergy",child1.E());//Real Missing gamma energy****************************
rmgamenergy=child1.E();
if(child1.E() > child2.E()){energyflag=1;}//More energetic photon missing
else{energyflag=-1;}//Less energetic photon missing
//nt->column("EnergyFlag",energyflag);**********************************************************
//Converted Photons
if(child1.daLast()-child1.daFirst()==1)
{
Gen_hepevt& grandchild3 = GenMgr(Panther_ID(child1.daFirst()));
Gen_hepevt& grandchild4 = GenMgr(Panther_ID(child1.daLast()));
if(fabs(grandchild3.idhep())==11 && fabs(grandchild4.idhep())==11){
//nt->column("EConVX",grandchild3.VX());***********************************
//nt->column("EConVY",grandchild3.VY());************************************88
//nt->column("EConVZ",grandchild3.VZ());************************************
econvx=grandchild3.VX();
econvy=grandchild3.VY();
econvz=grandchild3.VZ();
}
if(fabs(grandchild3.idhep())==13 && fabs(grandchild4.idhep())==13){
muconvx=grandchild3.VX();
muconvy=grandchild3.VY();
muconvz=grandchild3.VZ();
}


}

}
}
if(etaflag==1){
double px2,py2,pz2,p2,e2;
px2=GamGen.mother().PX();
py2=GamGen.mother().PY();
pz2=GamGen.mother().PZ();
p2=sqrt((px2*px2)+(py2*py2)+(pz2*pz2));
//nt->column("CEtaMom",p2);//Contaminating Eta//***********************
cetamom=p2;
//opening angle
Gen_hepevt& child3 = GenMgr(Panther_ID(GamGen.mother().daFirst()));
Gen_hepevt& child4 = GenMgr(Panther_ID(GamGen.mother().daLast()));
Hep3Vector A2(child3.PX(),child3.PY(),child3.PZ()); Hep3Vector B2(child4.PX(),child4.PY(),child4.PZ());
etaangle=A2.angle(B2);
//nt->column("EtaAngle",etaangle);//**************************************

}

}}

nt->column("CPiMom",cpimom);
nt->column("PiAngle",piangle);
nt->column("MGamEnergy",mgamenergy);
nt->column("RMGamEnergy",rmgamenergy);
nt->column("EnergyFlag",energyflag);
nt->column("EConVX",econvx);
nt->column("EConVY",econvy);
nt->column("EConVZ",econvz);
nt->column("MuConVX",muconvx);
nt->column("MuConVY",muconvy);
nt->column("MuConVZ",muconvz);
nt->column("CEtaMom",cetamom);
nt->column("EtaAngle",etaangle);




//good candidates
if(massdiff<0.16 && massdiff >0.14){
if(Dst_cand.child(0).child(0).mass() >1.01 && Dst_cand.child(0).child(0).mass() <1.03){
if(dstpps > 2.43){
if(Dst_cand.child(0).mass() > 1.65 && Dst_cand.child(0).mass() < 2.05){
if(Dst_cand.child(0).child(1).e() > 0.47){
if(aux.e9oe25() > 0.93){
if(ppiz < 0.2 && peta < 0.2){
if(piflag==1){

goodcandflag=1;

}}}}}}}}


      nt->column("PiFlag",piflag);
      nt->column("EtaFlag",etaflag);

//
if(KpGen && KmGen && PhiGen && GamGen && DzGen && PisGen && DstGen){//if all links exist

if(KpGen.idhep()==321 && KmGen.idhep()==-321 && PhiGen.idhep()==333 && GamGen.idhep()==22){//Are all the particles correctly identified?
if( (PisGen.idhep()==211 && DstGen.idhep()==413 && DzGen.idhep()==421 ) || (PisGen.idhep()==-211 && DstGen.idhep()==-413 && DzGen.idhep()==-421 )  ){

//Check all mother daughter relationships for phi gamma
if(KpGen.mother().get_ID()== PhiGen.get_ID() && KmGen.mother().get_ID()== PhiGen.get_ID()){
if(PhiGen.mother().get_ID()== DzGen.get_ID() && GamGen.mother().mother().get_ID()== DzGen.get_ID()){
if(DzGen.mother().get_ID()==DstGen.get_ID() && PisGen.mother().get_ID()==DstGen.get_ID() ){

//cout<<"+-+-+-+-+-+-+-+-+-+-+-+  SUCCESS!!!  +-+-+-+-+-+-+--+-+-+-+-+-+-+-"<<endl;
if( GamGen.mother().idhep()== 111){mypgenflag=1;}

}}}


//Check all mother daughter relationships for phi gamma
if(KpGen.mother().get_ID()== PhiGen.get_ID() && KmGen.mother().get_ID()== PhiGen.get_ID()){
if(PhiGen.mother().get_ID()== DzGen.get_ID() && GamGen.mother().get_ID()== DzGen.get_ID()){
if(DzGen.mother().get_ID()==DstGen.get_ID() && PisGen.mother().get_ID()==DstGen.get_ID() ){

//cout<<"+-+-+-+-+-+-+-+-+-+-+-+  SUCCESS!!!  +-+-+-+-+-+-+--+-+-+-+-+-+-+-"<<endl;
mygenflag=1;
}}}

}
}//Are all the particles correctly identified?


}//if all links exist
    if(mypgenflag!=1 && goodcandflag==1){PrintGenTable(Dst_cand,1,count1);}
    nt->column("MyPGenFlag",mypgenflag);
    nt->column("MyGenFlag",mygenflag);
}//if(dataType)

    nt->dumpData();
}//  void d0rad::FillTuple(Particle& Dst_cand, BelleTuple *nt)





void d0rad::PrintGenTable(Particle& Dst_cand,int GenTableFlag, int & counts)
{
int ref=GenTableFlag;
//Genstuff
int maxdau=1;
ofstream fout1;
int a1,a2,a3,a4,a5,a6,a7;
a1=Dst_cand.genHepevt().get_ID();                              //D*
a2=Dst_cand.child(0).genHepevt().get_ID();                     //D0   
a3=Dst_cand.child(1).genHepevt().get_ID();                     //sPi
a4=Dst_cand.child(0).child(0).genHepevt().get_ID();            //Phi
a5=Dst_cand.child(0).child(1).genHepevt().get_ID();            //Gamma
a6=Dst_cand.child(0).child(0).child(0).genHepevt().get_ID();   //K+
a7=Dst_cand.child(0).child(0).child(1).genHepevt().get_ID();   //K-


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

void d0rad::Transfer(int Temp[10],int IdMatrix[10000][10],int & nIM)
{
int i,j,k;     //i=0 no match , i=1 match

if(nIM==0){nIM++;

for(j=0;j<=6;j++){
IdMatrix[nIM][j]=Temp[j];
}
IdMatrix[nIM][7]=1;}
//---------------------------------------------------------------
else{// else1
int flag=0;
for(k=1;k<=nIM;k++){
i=0;
for(j=0;j<=6;j++){
if(Temp[j]==IdMatrix[k][j]){i++;} // no of matches increased
}// j loop

if(i==7){IdMatrix[k][7]=IdMatrix[k][7]+1; flag=1;} // all entries match

}// k loop

if(flag==0){nIM++;// else2

for(j=0;j<=6;j++){
IdMatrix[nIM][j]=Temp[j];
}
IdMatrix[nIM][7]=1;

 }// else2

}// else1

}

void d0rad::Append(int tmpmat[10000][10], int & tmpn, int IdMatrix[10000][10],int & nIM)
{
int i,j,k,k1;     //i=0 no match , i=1 match

for(k1=1;k1<=tmpn;k1++){ //looping over tmpmat sequences

int flag=0;
for(k=1;k<=nIM;k++){ //looping over IdMatrix sequences
i=0;
for(j=0;j<=6;j++){
if(tmpmat[k1][j]==IdMatrix[k][j]){i++;} // no of matches increased
}// j loop

if(i==7){IdMatrix[k][7]=IdMatrix[k][7]+tmpmat[k1][7]; flag=1;} // all entries match

}// k loop

if(flag==0){nIM++;// else2

for(j=0;j<=6;j++){
IdMatrix[nIM][j]=tmpmat[k1][j];
}
IdMatrix[nIM][7]=tmpmat[k1][7];

 }// else2

}//k1 loop

}


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

void d0rad::Print(int IdMatrix[10000][10],int & nIM)
{
int jobsdone;
ifstream fin; ofstream fout;
fin.open("matrix.txt");
fin>>jobsdone;
fin.close();
if(jobsdone==0){
jobsdone++;
fout.open("matrix.txt");
fout<<jobsdone<<endl;
fout<<"first\t"<<"second\t"<<"third\t"<<"fourth\t"<<"fifth\t"<<"sixth\t"<<"seventh\t"<<"frequency"<<endl;

for(int i=1;i<=nIM;i++){
for(int j=0;j<=7;j++){
fout<<IdMatrix[i][j]<<"\t";
}
fout<<endl;
}
fout.close();
}
else{//if jobsdone != 0
int tmpmat[10000][10],tmpn;
read(tmpmat,tmpn);
Append(tmpmat,tmpn,IdMatrix,nIM);
jobsdone++;
fout.open("matrix.txt");
fout<<jobsdone<<endl;
fout<<"first\t"<<"second\t"<<"third\t"<<"fourth\t"<<"fifth\t"<<"sixth\t"<<"seventh\t"<<"frequency"<<endl;

for(int i=1;i<=nIM;i++){
for(int j=0;j<=7;j++){
fout<<IdMatrix[i][j]<<"\t";
}
fout<<endl;
}
fout.close();
}//else
}

void d0rad::read(int tmpmat[10000][10],int & tmpn)
{
int count=0;
ifstream fin; string s;
fin.open("matrix.txt");
fin>>s>>s>>s>>s>>s>>s>>s>>s>>s;
while(!fin.eof()){
count++;
fin>>tmpmat[count][0]>>tmpmat[count][1]>>tmpmat[count][2]>>tmpmat[count][3]>>tmpmat[count][4]>>tmpmat[count][5]>>tmpmat[count][6]>>tmpmat[count][7];
}
tmpn=count-1;
fin.close();
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
