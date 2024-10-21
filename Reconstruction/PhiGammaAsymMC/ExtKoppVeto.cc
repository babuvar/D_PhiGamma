#include "belle.h"
//
// P.Koppenburg - 2003 08 29
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include "pi0eta_prob.h"
//#include "ana_utils.h"
//#include "properties.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//Extended Pi0 Veto
double Closest_Pi0_Probability_Ext
(int what, const Mdst_gamma& IN, Mdst_gamma& OUT, double& mm, Mdst_gamma* VETO)
{

  HepLorentzVector Reco1( IN.px(), IN.py(), IN.pz()  ) ;
  Reco1= HepLorentzVector( Reco1.px(), Reco1.py(), Reco1.pz(), Reco1.rho()) ;
  double bp = -100.;
  mm = -1;

  Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

  for (std::vector<Mdst_gamma>::iterator it1= GammaMgr.begin(); it1 !=GammaMgr.end(); ++it1){
    Mdst_gamma* D = &(*it1);
    if (VETO!=NULL){ if (D->get_ID()==VETO->get_ID()) continue ;}
    if (D->get_ID()==IN.get_ID()) continue ;

    HepLorentzVector PG( D->px(), D->py(), D->pz()  ) ;
    PG = HepLorentzVector( PG.px(), PG.py(), PG.pz(), PG.rho()) ;
    double m = (Reco1 + PG).m();
    double p = Pi0_Eta_Prob(what, m, PG.e(), PG.theta() );

    if (p>bp) {
      bp = p ;
      mm = m ;
      OUT = *D ;
    }
  }
  return bp;
}




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



                                    
#if defined(BELLE_NAMESPACE)
}// namecpace Belle
#endif

   
