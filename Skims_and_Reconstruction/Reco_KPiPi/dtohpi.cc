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
#include "dtohpi.h"
#include "benergy/BeamEnergy.h"
#include "eid/eid.h"

#include "pi0eta_prob.h"
#include "geninfo.h"


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

using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


  extern "C" Module_descr *mdcl_dtohpi()
  {
    dtohpi *module = new dtohpi;
    Module_descr *dscr = new Module_descr ( "dtohpi", module );
    return dscr;
  }



dtohpi::dtohpi( void ) {

//  std::strcpy( m_SkimFileName,   "dtohpi.index" );

}

void dtohpi::disp_stat( const char* ) {
}

void dtohpi::end_run( BelleEvent*, int* ) {
}

void dtohpi::other( int*, BelleEvent*, int* ) {
}

//Particle Types
  const Ptype dtohpi::m_ptypeDP("D+");
  const Ptype dtohpi::m_ptypeDM("D-");
  const Ptype dtohpi::m_ptypeDstarP("D*+");
  const Ptype dtohpi::m_ptypeDstarM("D*-");

//  const Ptype dtohpi::m_ptypeDsP(431);
//  const Ptype dtohpi::m_ptypeDsM(-431);
//  const Ptype dtohpi::m_ptypeDs_starP("433");
//  const Ptype dtohpi::m_ptypeDs_starM("-433");

void dtohpi::init ( int * ) {
//  extern BasfOutputManager* BASF_Output;
//  m_SkimFile = BASF_Output->open ( m_SkimFileName );

}

void dtohpi::begin_run( BelleEvent*, int* ) {
  
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


void dtohpi::term( void ) {

//  delete m_SkimFile;

}

void dtohpi::hist_def(void)
{
  // no histograms
  extern BelleTupleManager *BASF_Histogram;
  std::string title;


  nt_dpipi_small = BASF_Histogram->ntuple("information for mode d to pipi with deltam","DMass DCharge DeltaM Gam1The Gam2The DStarMass PstDst PstD PPiZs PstPiZs GamHeli PizHeli EgamAsy Categ EGam1s EGam2s  PiZsMass Pi1Mom Pi2Mom KMom DstID DID DstF DF", 1);


 

}

void dtohpi::event( BelleEvent *, int *status ) {

  *status = 0;
  // ---------------- particle reconstruction and selection 

  vector<Particle>   k_p, k_m, pi_p, pi_m;
  makeKPi(k_p, k_m, pi_p, pi_m, 0);

  with_imp_cut(k_p);
  with_imp_cut(k_m);
  with_imp_cut(pi_p);
  with_imp_cut(pi_m);

  withKaonId(k_p,0.6,3,1,5);
  withKaonId(k_m,0.6,3,1,5);
//  withPionId(pi_p,0.6,3,1,5);
//  withPionId(pi_m,0.6,3,1,5);

//  withPCut(pi_p, 0.75);
//  withPCut(pi_m, 0.75);

  vector<Particle>  pi0;
  makePi0(pi0);

  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
    if(i->mdstPi0().gamma(0).ecl().energy()<0.03||
       i->mdstPi0().gamma(1).ecl().energy()<0.03||
       fabs(i->mdstPi0().mass()-.135)>0.025 ) {
      pi0.erase(i); 
      --i;
    }
  
                                                                  

  vector<Particle> D_Plus, D_Minus;


  // D+  ->   pi0  h+
  combination(D_Plus, Ptype(411), k_m, pi_p, pi_p, 0.08);
  // D-  ->   pi0  h-
  combination(D_Minus, Ptype(-411), k_p, pi_m, pi_m, 0.08);


  // D* -> D pi0
  vector<Particle> DstP, DstM;

  combination( DstP, Ptype(413), D_Plus, pi0);
  combination( DstM, Ptype(-413), D_Minus, pi0);
  withMassDifCut( DstP, 0.13, 0.16, 0);
  withMassDifCut( DstM, 0.13, 0.16, 0);
  withPSCut( DstP, 2.5);
  withPSCut( DstM, 2.5);






  for (int i=0;i<DstP.size();i++) {
FillTuple_small_D(DstP[i],nt_dpipi_small);

}
  for (int i=0;i<DstM.size();i++) {
FillTuple_small_D(DstM[i],nt_dpipi_small);

}




//    *status = 1;
  
}

void dtohpi::FillTuple_small_D(Particle& Dst_cand, BelleTuple *nt){


//Fill stuff
double massdiff, E_Asy;

      nt->column("DStarMass",Dst_cand.mass());
      nt->column("DCharge",Dst_cand.child(0).charge());
      massdiff=Dst_cand.mass()-Dst_cand.child(0).mass();
      nt->column("DeltaM",massdiff);
      nt->column("DMass",Dst_cand.child(0).mass());
      nt->column("PstDst",pStar(Dst_cand).vect().mag());
      nt->column("PstD",pStar(Dst_cand.child(0)).vect().mag());

      nt->column("KPID",kaonId(Dst_cand.child(0).child(0)));
      nt->column("Pi1PID",kaonId(Dst_cand.child(0).child(1)));
      nt->column("Pi2PID",kaonId(Dst_cand.child(0).child(2)));
      nt->column("EGam1s",Dst_cand.child(1).child(0).e());
      nt->column("EGam2s",Dst_cand.child(1).child(1).e());


      nt->column("PiZsMass",Dst_cand.child(1).mdstPi0().mass());
      nt->column("Pi1Mom",Dst_cand.child(0).child(1).ptot());
      nt->column("Pi2Mom",Dst_cand.child(0).child(1).ptot());
      nt->column("KMom",Dst_cand.child(0).child(1).ptot());

      nt->column("PstPiZs",pStar(Dst_cand.child(1)).vect().mag());
      nt->column("PPiZs",Dst_cand.child(1).p().vect().mag());

//Gamma Helicity
HepLorentzVector dst4V, piz4V, gam4V, dst4V_st;
Hep3Vector piz_boost, dst_boost; double cos_helicity;

dst4V=Dst_cand.p();		//D*  mom
piz4V=Dst_cand.child(1).p();	//pi0 mom
gam4V=Dst_cand.child(1).child(0).p();   //1st Gamma mom
piz_boost=-1*piz4V.boostVector();
dst4V.boost(piz_boost);
gam4V.boost(piz_boost);
cos_helicity=cos(dst4V.angle(gam4V));
      nt->column("GamHeli",cos_helicity);

//Pi0 Helicity
dst4V=Dst_cand.p();
dst_boost=-1*dst4V.boostVector();
dst4V_st=pStar(Dst_cand);
piz4V=Dst_cand.child(1).p();
piz4V.boost(dst_boost);
cos_helicity=cos(dst4V_st.angle(piz4V));
      nt->column("PizHeli",cos_helicity);


//Photon energy asymmetry
double f_EgamAsy;
double factor = Dst_cand.child(1).e() / Dst_cand.child(1).ptot();

f_EgamAsy=( Dst_cand.child(1).child(0).e() - Dst_cand.child(1).child(1).e() ) / ( Dst_cand.child(1).child(0).e() + Dst_cand.child(1).child(1).e() );
f_EgamAsy=f_EgamAsy*factor;

      nt->column("EgamAsy",f_EgamAsy);

//Categories
//Photon thetas
double theta1, theta2; int theta1_flag, theta2_flag, category;


Hep3Vector The1(Dst_cand.child(1).child(0).mdstGamma().px(),Dst_cand.child(1).child(0).mdstGamma().py(),Dst_cand.child(1).child(0).mdstGamma().pz());
theta1=The1.theta()*180/3.14159265359;
      nt->column("Gam1The",theta1);

Hep3Vector The2(Dst_cand.child(1).child(1).mdstGamma().px(),Dst_cand.child(1).child(1).mdstGamma().py(),Dst_cand.child(1).child(1).mdstGamma().pz());
theta2=The2.theta()*180/3.14159265359;
      nt->column("Gam2The",theta2);


if(theta1 < 32.2){theta1_flag=-1;} //FW
else if(theta1 > 128.7){theta1_flag=1;} //BW
else{theta1_flag=0;} //BARREL


if(theta2 < 32.2){theta2_flag=-1;} //FW
else if(theta2 > 128.7){theta2_flag=1;} //BW
else{theta2_flag=0;} //BARREL


category =-1; //Default

if(theta1_flag == -1 && theta2_flag == -1){category =1;}
if(theta1_flag == -1 && theta2_flag == 0){category =2;}
if(theta1_flag == -1 && theta2_flag == 1){category =3;}
if(theta1_flag == 0 && theta2_flag == -1){category =4;}
if(theta1_flag == 0 && theta2_flag == 0){category =5;}
if(theta1_flag == 0 && theta2_flag == 1){category =6;}
if(theta1_flag == 1 && theta2_flag == -1){category =7;}
if(theta1_flag == 1 && theta2_flag == 0){category =8;}
if(theta1_flag == 1 && theta2_flag == 1){category =9;}

      nt->column("Categ",category);



if(dataType){
      setMCtruth(Dst_cand);
      nt->column("DstID",IDhep(Dst_cand));
      nt->column("DID",IDhep(Dst_cand.child(0)));
      nt->column("DstF",getMCtruthFlag(Dst_cand));
      nt->column("DF",getMCtruthFlag(Dst_cand.child(0)));
}


      nt->dumpData();
}



  void dtohpi::with_imp_cut(std::vector<Particle> &list) {
    for(int i=0;i<(int)list.size();++i){
      if(list[i].mdstCharged()){
	HepVector a(param_at_ip(list[i]));
	if (!(abs(a[0])<1.0 && abs(a[3])<3.0)){
	  list.erase(list.begin()+i);
	  --i;
	}
      }
    }
  }


  HepVector dtohpi::param_at_ip(Particle &p){

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
