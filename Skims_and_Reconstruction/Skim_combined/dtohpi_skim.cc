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
#include "dtohpi_skim.h"
#include "benergy/BeamEnergy.h"
#include "eid/eid.h"

//#include "pi0eta_prob.h"
//#include "geninfo.h"



// system include files

#include <cmath>
#include <string>
#include <cstring>
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


  extern "C" Module_descr *mdcl_dtohpi_skim()
  {
    dtohpi_skim *module = new dtohpi_skim;
    Module_descr *dscr = new Module_descr ( "dtohpi_skim", module );
    dscr->define_param ( "parm", "parameter to select skim and optimized cuts ",  &module->parm_ );
    dscr->define_param ( "parm2", "parameter to select skim and optimized cuts 2" , &module->parm2_ );
    dscr->define_param ( "parm3", "parameter to select skim and optimized cuts 3" , &module->parm3_ );
//    dscr->define_param ( "parm4", "parameter to select skim and optimized cuts 4" ,10, module->parm4_ );	//MC
//    dscr->define_param ( "parm5", "parameter to select skim and optimized cuts 5" , &module->parm5_ ); 	//MC

    return dscr;
  }




dtohpi_skim::dtohpi_skim( void ) {

  std::strcpy( m_SkimFileName,   "dtohpi_skim.index" );

}

void dtohpi_skim::disp_stat( const char* ) {
}

void dtohpi_skim::end_run( BelleEvent*, int* ) {
}

void dtohpi_skim::other( int*, BelleEvent*, int* ) {
}

//Particle Types
  const Ptype dtohpi_skim::m_ptypeDP("D+");
  const Ptype dtohpi_skim::m_ptypeDM("D-");
  const Ptype dtohpi_skim::m_ptypeDstarP("D*+");
  const Ptype dtohpi_skim::m_ptypeDstarM("D*-");


void dtohpi_skim::init ( int * ) {
//Data
//char * str1 = "index_Data_Y5S/exp";
char * str1 = "index_Data/exp";
char * str2 = "_";
char * str3 = "_Data.index";
  stringstream strs;
  strs <<str1<<parm3_<<str2<<parm_<<str2<<parm2_<<str3;


//MC
/*
char * str1 = "index_MC";
char * str1p = "/exp";
char * str2 = "_";
char * str3 = "_MC.index";
  stringstream strs;
  strs <<str1<<parm5_<<str1p<<parm3_<<str2<<parm_<<str2<<parm2_<<str2<<parm4_<<str3;
*/

//Common
  string temp_str = strs.str();
  char* str = (char*) temp_str.c_str();
  std::strcpy( m_SkimFileName, str );

cout<<" name of output file is "<<m_SkimFileName<<endl; 
  extern BasfOutputManager* BASF_Output;
  m_SkimFile = BASF_Output->open ( m_SkimFileName );
}

void dtohpi_skim::begin_run( BelleEvent*, int* ) {
  
  IpProfile::begin_run();
    m_SkimFile->write();

}


void dtohpi_skim::term( void ) {

  delete m_SkimFile;

}

void dtohpi_skim::hist_def(void)
{
  // no histograms
}

void dtohpi_skim::event( BelleEvent *, int *status ) {

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




  vector<Particle>    pi_p_H, pi_m_H; //High momemtum charged pion lists

  for(vector<Particle>::iterator i=pi_p.begin(); i!=pi_p.end();++i){
    Particle pi_chP(*i);
//if(pi_chP.ptot() > 0.75){pi_p_H.push_back(pi_chP);}
if(pi_chP.ptot() > 0.5){pi_p_H.push_back(pi_chP);}
}
  for(vector<Particle>::iterator i=pi_m.begin(); i!=pi_m.end();++i){
    Particle pi_chM(*i);
//if(pi_chM.ptot() > 0.75){pi_m_H.push_back(pi_chM);}
if(pi_chM.ptot() > 0.5){pi_m_H.push_back(pi_chM);}
}


  withPionId(pi_p_H,0.6,3,1,5);
  withPionId(pi_m_H,0.6,3,1,5);



  vector<Particle>  pi0;
  makePi0(pi0);

  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i)
    if(i->mdstPi0().gamma(0).ecl().energy()<0.03||
       i->mdstPi0().gamma(1).ecl().energy()<0.03||
       fabs(i->mdstPi0().mass()-.135)>0.025 ) {
      pi0.erase(i); 
      --i;
    }
                                                                      

  vector<Particle>     pi0_H; //High momemtum neutral pion list

  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i){
    Particle pi0H(*i);
if(pi0H.ptot() > 0.75 && pi0H.mdstPi0().gamma(0).ecl().energy()>0.05 && pi0H.mdstPi0().gamma(1).ecl().energy()>0.05){pi0_H.push_back(pi0H);}
}
   




  vector<Particle> D_Plus_pp, D_Plus_kpp, D_Minus_pp, D_Minus_kpp;


  // D+  ->   pi0  pi+
  combination(D_Plus_pp, m_ptypeDP, pi0_H, pi_p_H, 0.2);
  // D-  ->   pi0  pi-
  combination(D_Minus_pp, m_ptypeDM, pi0_H, pi_m_H, 0.2);


  withPSCut( D_Plus_pp, 2.0);
  withPSCut( D_Minus_pp, 2.0);



// D+  ->   K- pi+ pi+  
  combination(D_Plus_kpp, m_ptypeDP, k_m, pi_p, pi_p, 0.08);
  // D-  ->   K+ pi- pi-
  combination(D_Minus_kpp, m_ptypeDM, k_p, pi_m, pi_m, 0.08);


  vector<Particle> DstP_kpp, DstM_kpp;

  combination( DstP_kpp, m_ptypeDstarP, D_Plus_kpp, pi0);
  combination( DstM_kpp, m_ptypeDstarM, D_Minus_kpp, pi0);
  withMassDifCut( DstP_kpp, 0.13, 0.16, 0);
  withMassDifCut( DstM_kpp, 0.13, 0.16, 0);
  withPSCut( DstP_kpp, 2.5);
  withPSCut( DstM_kpp, 2.5);



//cout<<"General Event"<<endl;


if(D_Plus_pp.size()+D_Minus_pp.size()+DstP_kpp.size()+DstM_kpp.size() > 0){
    m_SkimFile->write();
    *status = 1;
//cout<<"Skimmed Event"<<endl;
}

  


}





  void dtohpi_skim::with_imp_cut(std::vector<Particle> &list) {
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


  HepVector dtohpi_skim::param_at_ip(Particle &p){

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
