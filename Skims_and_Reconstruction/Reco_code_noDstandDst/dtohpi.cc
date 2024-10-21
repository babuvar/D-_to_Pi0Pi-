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

    BeamEnergy::define_global(dscr);

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


  nt_dpipi_small = BASF_Histogram->ntuple("information for mode d to pipi with deltam","DMass DCharge DeltaM DStarMass PstDst CThetDst PstD PstPiZs Heli Gam1hThe Gam2hThe Gam1sThe Gam2sThe Categ GamHeli PizHeli EgamAsy EGamma1 EGamma2 EGam1s EGam2s AsyEGam PiZMass PiZsMass PiZMom  PiZsMom PiMom GenPiMom GenPiZMom PiCThet PiDR PiDZ PiPID DstID DID PiZ0ID PiZ1ID PiID DstF DF PiZ0F PiZ1F PiF MultFlag TruEvt Exp Run Evt GenPizID GenPiID ", 1);


  nt_dpipi_big = BASF_Histogram->ntuple("information for mode d to pipi","DMass DCharge  PstD  Gam1hThe Gam2hThe EGamma1 EGamma2 PiZMass  PiZMom  PiMom  PiCThet PiDR PiDZ PiPID  DID PiZID PiID  DF PiZF  PiF MultFlag TruEvt  Exp Run Evt GenPizID GenPiID ", 2);


  nt_mult_before =  BASF_Histogram->ntuple("Multiciplity info before selection criteria","PMult NMult TMult", 3);
  nt_mult_after =  BASF_Histogram->ntuple("Multiciplity info after selection criteria","PMult NMult TMult", 4);

  d_mult_after =  BASF_Histogram->ntuple("D multiciplity info after selection criteria","PMult NMult TMult", 5);



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
  withPionId(pi_p,0.6,3,1,5);
  withPionId(pi_m,0.6,3,1,5);

  withPCut(pi_p, 0.75);
  withPCut(pi_m, 0.75);
//  withPCut(pi_p, 0.5);
//  withPCut(pi_m, 0.5);



  vector<Particle>  pi0, pi0_high;
  makePi0(pi0);

  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i){
    if(//i->mdstPi0().gamma(0).ecl().energy()<0.05||
       //i->mdstPi0().gamma(1).ecl().energy()<0.05||
       //fabs(i->mdstPi0().mass()-0.135)>0.025 || i->ptot()<0.75) {
//       fabs(i->mdstPi0().mass()-0.135)>0.025 ) {
       fabs(i->mdstPi0().mass()-0.135)>0.016 ) {
      pi0.erase(i); 
      --i;
    }}
  
  for(vector<Particle>::iterator i=pi0.begin(); i!=pi0.end();++i){
    Particle pi0H(*i);
if(pi0H.ptot() > 0.75){pi0_high.push_back(pi0H);}
//pi0_high.push_back(pi0H);
}


/*  vector<Particle> k_s;
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
  i}*/

  vector<Particle> gamma;
  vector<Particle> gamma_low;

  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin(); i != gamma_mag.end(); ++i){
    Particle gammaP(*i);
    
    if(gammaP.e()<0.050)
      continue;
    
    gamma_low.push_back(gammaP);

    if(gammaP.e()>0.050)
      gamma.push_back(gammaP);
  }

//  vector<Particle> phi;
//  combination( phi, Ptype(333),k_p,k_m);
//  withMassCut( phi, 0.9, 1.06);
  

//  vector<Particle> eta;
//  combination( eta, Ptype(221),gamma_low,gamma_low);
//  withMassCut( eta, 0.45, 0.65);
//  withPCut( eta, 0.3);

  vector<Particle> D_Plus_Pi, D_Plus_K, D_Minus_Pi, D_Minus_K;
  vector<Particle> Ds_Plus_Pi, Ds_Minus_Pi, Ds_Plus_K, Ds_Minus_K;

  // D+  ->   pi0  h+
  combination(D_Plus_Pi, Ptype(411), pi0_high, pi_p, 0.2);
  // D-  ->   pi0  h-
  combination(D_Minus_Pi, Ptype(-411), pi0_high, pi_m, 0.2);





  // D* -> D pi0
  vector<Particle> DstP_Pi, DstM_Pi;

  combination( DstP_Pi, Ptype(413), D_Plus_Pi, pi0);
  combination( DstM_Pi, Ptype(-413), D_Minus_Pi, pi0);
//  withMassDifCut( DstP_Pi, 0.135, 0.155, 0);
//  withMassDifCut( DstM_Pi, 0.135, 0.155, 0);
  withMassDifCut( DstP_Pi, 0.139, 0.143, 0);
  withMassDifCut( DstM_Pi, 0.139, 0.143, 0);
//  withPSCut( DstP_Pi, 2.0);
//  withPSCut( DstM_Pi, 2.0);
  withPSCut( DstP_Pi, 2.5);
  withPSCut( DstM_Pi, 2.5);

/*
//FillTuple without Best Candidate Selection
  for (int i=0;i<DstP_Pi.size();i++) {
FillTuple_small_D(DstP_Pi[i],nt_dpipi_small);
}
  for (int i=0;i<DstM_Pi.size();i++) {
FillTuple_small_D(DstM_Pi[i],nt_dpipi_small);
}

//    *status = 1;
*/


if(DstM_Pi.size()+DstP_Pi.size() > 0){
nt_mult_before->column("PMult",DstP_Pi.size());
nt_mult_before->column("NMult",DstM_Pi.size());
nt_mult_before->column("TMult",DstP_Pi.size()+DstM_Pi.size());
nt_mult_before->dumpData();
}



//Cuts before BCS
int GoodCand_P[50]={0}, GoodCand_M[50]={0};
int photon1cutflag, photon2cutflag, sphoton1cutflag, sphoton2cutflag;
double f_Deltam,f_dzero,f_Pizsmass,f_Pizmass,f_Pizmom,f_Pimom,f_PD,f_Categ,f_Gam1thet,f_Egamma1,f_Gam2thet,f_Egamma2,f_Egam1s,f_Egam2s;
//For P
  for (int i=0;i<DstP_Pi.size();i++) {
 photon1cutflag=0; photon2cutflag=0; sphoton1cutflag=0; sphoton2cutflag=0;

//Defining Variables
f_Deltam=DstP_Pi[i].mass()-DstP_Pi[i].child(0).mass();
f_dzero=DstP_Pi[i].child(0).mass();
f_Pizsmass=DstP_Pi[i].child(1).mdstPi0().mass();
f_Pizmass=DstP_Pi[i].child(0).child(0).mdstPi0().mass();
f_Pizmom=DstP_Pi[i].child(0).child(0).ptot();
f_Pimom=DstP_Pi[i].child(0).child(1).ptot();
f_PD=pStar(DstP_Pi[i]).vect().mag(); 
f_Egamma1=DstP_Pi[i].child(0).child(0).child(0).e();
f_Egamma2=DstP_Pi[i].child(0).child(0).child(1).e();
f_Egam1s=DstP_Pi[i].child(1).child(0).e();
f_Egam2s=DstP_Pi[i].child(1).child(1).e();

int theta1_flag, theta2_flag, category;
double theta1,theta2;

Hep3Vector The1h(DstP_Pi[i].child(0).child(0).child(0).mdstGamma().px(),DstP_Pi[i].child(0).child(0).child(0).mdstGamma().py(),DstP_Pi[i].child(0).child(0).child(0).mdstGamma().pz());
theta1=The1h.theta()*180/3.14159265359;
f_Gam1thet=theta1;
Hep3Vector The2h(DstP_Pi[i].child(0).child(0).child(1).mdstGamma().px(),DstP_Pi[i].child(0).child(0).child(1).mdstGamma().py(),DstP_Pi[i].child(0).child(0).child(1).mdstGamma().pz());
theta2=The2h.theta()*180/3.14159265359;
f_Gam2thet=theta2;


Hep3Vector The1s(DstP_Pi[i].child(1).child(0).mdstGamma().px(),DstP_Pi[i].child(1).child(0).mdstGamma().py(),DstP_Pi[i].child(1).child(0).mdstGamma().pz());
theta1=The1s.theta()*180/3.14159265359;
Hep3Vector The2s(DstP_Pi[i].child(1).child(1).mdstGamma().px(),DstP_Pi[i].child(1).child(1).mdstGamma().py(),DstP_Pi[i].child(1).child(1).mdstGamma().pz());
theta2=The2s.theta()*180/3.14159265359;



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

f_Categ=category;


if(f_Deltam > 0.139 && f_Deltam < 0.142){       //optimized for M_D fitting
//if(f_Deltam > 0.135 && f_Deltam < 0.155){
//if(f_dzero >1.68  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
 if(f_dzero >1.68  && f_dzero < 2.06){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){

//Hard Photon 1
if(f_Gam1thet < -60 && f_Egamma1 > 0.150){photon1cutflag=1;}  //BWD
else if(f_Gam1thet > 73 && f_Egamma1 > 0.100){photon1cutflag=1;}  //FWD
else if(f_Gam1thet > -60 && f_Gam1thet < 73 && f_Egamma1 > 0.050){photon1cutflag=1;}  //BARREL
else{photon1cutflag=0;}

//Hard Photon 2
if(f_Gam2thet < -60 && f_Egamma2 > 0.150){photon2cutflag=1;}  //BWD
else if(f_Gam2thet > 73 && f_Egamma2 > 0.100){photon2cutflag=1;}  //FWD
else if(f_Gam2thet > -60 && f_Gam2thet < 73 && f_Egamma2 > 0.050){photon2cutflag=1;}  //BARREL
else{photon2cutflag=0;}

//Soft Photon 1 
if(f_Categ == 5 && f_Egam1s > 0.046){sphoton1cutflag=1;}
else if(f_Categ == 2 && f_Egam1s > 0.068){sphoton1cutflag=1;}
else if(f_Categ == 6 && f_Egam1s > 0.030){sphoton1cutflag=1;}
else{sphoton1cutflag=0;}

//Soft Photon 2 
if(f_Categ == 5 && f_Egam2s > 0.046){sphoton2cutflag=1;}
else if(f_Categ == 2 && f_Egam2s > 0.036){sphoton2cutflag=1;}
else if(f_Categ == 6 && f_Egam2s > 0.044){sphoton2cutflag=1;}
else{sphoton2cutflag=0;}


//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){
if(sphoton1cutflag == 1 && sphoton2cutflag == 1){


GoodCand_P[i]=1;


}}//Photon cuts

}}}}}}}

}//DstP Loop







//For N
for(int i=0;i<DstM_Pi.size();i++) {
   photon1cutflag=0; photon2cutflag=0; sphoton1cutflag=0; sphoton2cutflag=0;

   //Defining Variables
   f_Deltam=DstM_Pi[i].mass()-DstM_Pi[i].child(0).mass();
   f_dzero=DstM_Pi[i].child(0).mass();
   f_Pizsmass=DstM_Pi[i].child(1).mdstPi0().mass();
   f_Pizmass=DstM_Pi[i].child(0).child(0).mdstPi0().mass();
   f_Pizmom=DstM_Pi[i].child(0).child(0).ptot();
   f_Pimom=DstM_Pi[i].child(0).child(1).ptot();
   f_PD=pStar(DstM_Pi[i]).vect().mag();
   f_Egamma1=DstM_Pi[i].child(0).child(0).child(0).e();
   f_Egamma2=DstM_Pi[i].child(0).child(0).child(1).e();
   f_Egam1s=DstM_Pi[i].child(1).child(0).e();
   f_Egam2s=DstM_Pi[i].child(1).child(1).e();

   int theta1_flag, theta2_flag, category;
   double theta1,theta2;

   Hep3Vector The1h(DstM_Pi[i].child(0).child(0).child(0).mdstGamma().px(),DstM_Pi[i].child(0).child(0).child(0).mdstGamma().py(),DstM_Pi[i].child(0).child(0).child(0).mdstGamma().pz());
   theta1=The1h.theta()*180/3.14159265359;
   f_Gam1thet=theta1;
   Hep3Vector The2h(DstM_Pi[i].child(0).child(0).child(1).mdstGamma().px(),DstM_Pi[i].child(0).child(0).child(1).mdstGamma().py(),DstM_Pi[i].child(0).child(0).child(1).mdstGamma().pz());
   theta2=The2h.theta()*180/3.14159265359;
   f_Gam2thet=theta2;


   Hep3Vector The1s(DstM_Pi[i].child(1).child(0).mdstGamma().px(),DstM_Pi[i].child(1).child(0).mdstGamma().py(),DstM_Pi[i].child(1).child(0).mdstGamma().pz());
   theta1=The1s.theta()*180/3.14159265359;
   Hep3Vector The2s(DstM_Pi[i].child(1).child(1).mdstGamma().px(),DstM_Pi[i].child(1).child(1).mdstGamma().py(),DstM_Pi[i].child(1).child(1).mdstGamma().pz());
   theta2=The2s.theta()*180/3.14159265359;


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

   f_Categ=category;


   if(f_Deltam > 0.139 && f_Deltam < 0.142){       //optimized for M_D fitting
//   if(f_Deltam > 0.135 && f_Deltam < 0.155){
//   if(f_dzero >1.68  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
   if(f_dzero >1.68  && f_dzero < 2.06){ 
  if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
   if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
   if(f_Pizmom > 1.06 ){
   if(f_Pimom > 0.84 ){
   if(f_PD > 2.5){


   //Hard Photon 1
   if(f_Gam1thet < -60 && f_Egamma1 > 0.150){photon1cutflag=1;}  //BWD
   else if(f_Gam1thet > 73 && f_Egamma1 > 0.100){photon1cutflag=1;}  //FWD
   else if(f_Gam1thet > -60 && f_Gam1thet < 73 && f_Egamma1 > 0.050){photon1cutflag=1;}  //BARREL
   else{photon1cutflag=0;}

   //Hard Photon 2
   if(f_Gam2thet < -60 && f_Egamma2 > 0.150){photon2cutflag=1;}  //BWD
   else if(f_Gam2thet > 73 && f_Egamma2 > 0.100){photon2cutflag=1;}  //FWD
   else if(f_Gam2thet > -60 && f_Gam2thet < 73 && f_Egamma2 > 0.050){photon2cutflag=1;}  //BARREL
   else{photon2cutflag=0;}

   //Soft Photon 1 
   if(f_Categ == 5 && f_Egam1s > 0.046){sphoton1cutflag=1;}
   else if(f_Categ == 2 && f_Egam1s > 0.068){sphoton1cutflag=1;}
   else if(f_Categ == 6 && f_Egam1s > 0.030){sphoton1cutflag=1;}
   else{sphoton1cutflag=0;}

   //Soft Photon 2 
   if(f_Categ == 5 && f_Egam2s > 0.046){sphoton2cutflag=1;}
   else if(f_Categ == 2 && f_Egam2s > 0.036){sphoton2cutflag=1;}
   else if(f_Categ == 6 && f_Egam2s > 0.044){sphoton2cutflag=1;}
   else{sphoton2cutflag=0;}


   //Photon cuts
   if(photon1cutflag == 1 && photon2cutflag == 1){
   if(sphoton1cutflag == 1 && sphoton2cutflag == 1){


   GoodCand_M[i]=1;


   }}//Photon cuts

   }}}}}}}

   }//DstM Loop



//FillTuple with temporary Best Candidate Selection
//Step 1: Select the Best Candidate
double DeltaM, DeltaM_FOM, DeltaM_FOM_best=1.0;
int PorM_Flag=0, i_best; // 1 indicates DstP, -1 indicates DstM
int GoodPGuys=0, GoodMGuys=0;

 TrueCand_Evt=0;

//DstP loop
  for (int i=0;i<DstP_Pi.size();i++) {
//See if any true D* candidate is there in the event
      setMCtruth(DstP_Pi[i]);
      if(getMCtruthFlag(DstP_Pi[i])==1 && GoodCand_P[i]==1){TrueCand_Evt=1;}


if(GoodCand_P[i]==1){GoodPGuys++;
DeltaM= DstP_Pi[i].mass() - DstP_Pi[i].child(0).mass();
DeltaM_FOM = fabs(DeltaM-0.14069);
//Check best
if(DeltaM_FOM < DeltaM_FOM_best){
DeltaM_FOM_best= DeltaM_FOM;
i_best=i;
PorM_Flag=1;
}
}//Goodcand
}//DstP loop

//DstM loop
  for (int i=0;i<DstM_Pi.size();i++) {
//See if any true D* candidate is there in the event
      setMCtruth(DstM_Pi[i]);
      if(getMCtruthFlag(DstM_Pi[i])==1 &&  GoodCand_M[i]==1){TrueCand_Evt=1;}


if(GoodCand_M[i]==1){GoodMGuys++;
DeltaM= DstM_Pi[i].mass() - DstM_Pi[i].child(0).mass();
DeltaM_FOM = fabs(DeltaM-0.14069);
//Check best
if(DeltaM_FOM < DeltaM_FOM_best){
DeltaM_FOM_best= DeltaM_FOM;
i_best=i;
PorM_Flag=-1;
}
}//Goodcand
}//DstM loop

if(GoodPGuys > 0 || GoodMGuys > 0){
nt_mult_after->column("PMult",GoodPGuys);
nt_mult_after->column("NMult",GoodMGuys);
nt_mult_after->column("TMult",GoodPGuys+GoodMGuys);
nt_mult_after->dumpData();
}


mult_flag=GoodPGuys+GoodMGuys;

//Step 2: Fill nTuple 
if(PorM_Flag == 1){FillTuple_small_D(DstP_Pi[i_best],nt_dpipi_small);}
else if(PorM_Flag == -1){FillTuple_small_D(DstM_Pi[i_best],nt_dpipi_small);}



//-------------------------------------------------------------------------------
if(PorM_Flag == 0){
//D meson decay
  withPSCut( D_Plus_Pi, 2.5);
  withPSCut( D_Minus_Pi, 2.5);


//Cuts before BCS
//For P
  for (int i=0;i<D_Plus_Pi.size();i++) {
photon1cutflag=0; photon2cutflag=0;

//Defining Variables
f_dzero=D_Plus_Pi[i].mass();
f_Pizmass=D_Plus_Pi[i].child(0).mdstPi0().mass();
f_Pizmom=D_Plus_Pi[i].child(0).ptot();
f_Pimom=D_Plus_Pi[i].child(1).ptot();
f_PD=pStar(D_Plus_Pi[i]).vect().mag();
f_Egamma1=D_Plus_Pi[i].child(0).child(0).e();
f_Egamma2=D_Plus_Pi[i].child(0).child(1).e();

double theta1,theta2;


Hep3Vector The1h(D_Plus_Pi[i].child(0).child(0).mdstGamma().px(),D_Plus_Pi[i].child(0).child(0).mdstGamma().py(),D_Plus_Pi[i].child(0).child(0).mdstGamma().pz());
theta1=The1h.theta()*180/3.14159265359;
f_Gam1thet=theta1;

Hep3Vector The2h(D_Plus_Pi[i].child(0).child(1).mdstGamma().px(),D_Plus_Pi[i].child(0).child(1).mdstGamma().py(),D_Plus_Pi[i].child(0).child(1).mdstGamma().pz());
theta2=The2h.theta()*180/3.14159265359;
f_Gam2thet=theta2;

 if(f_dzero >1.68  && f_dzero < 2.06){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){

//Hard Photon 1
if(f_Gam1thet < -60 && f_Egamma1 > 0.150){photon1cutflag=1;}  //BWD
else if(f_Gam1thet > 73 && f_Egamma1 > 0.100){photon1cutflag=1;}  //FWD
else if(f_Gam1thet > -60 && f_Gam1thet < 73 && f_Egamma1 > 0.050){photon1cutflag=1;}  //BARREL
else{photon1cutflag=0;}

//Hard Photon 2
if(f_Gam2thet < -60 && f_Egamma2 > 0.150){photon2cutflag=1;}  //BWD
else if(f_Gam2thet > 73 && f_Egamma2 > 0.100){photon2cutflag=1;}  //FWD
else if(f_Gam2thet > -60 && f_Gam2thet < 73 && f_Egamma2 > 0.050){photon2cutflag=1;}  //BARREL
else{photon2cutflag=0;}

//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){

GoodCand_P[i]=1;

}//Photon cuts


}}}}}

}//D_Plus_Pi Loop



//For N
  for (int i=0;i<D_Minus_Pi.size();i++) {
photon1cutflag=0; photon2cutflag=0;

//Defining Variables
f_dzero=D_Minus_Pi[i].mass();
f_Pizmass=D_Minus_Pi[i].child(0).mdstPi0().mass();
f_Pizmom=D_Minus_Pi[i].child(0).ptot();
f_Pimom=D_Minus_Pi[i].child(1).ptot();
f_PD=pStar(D_Minus_Pi[i]).vect().mag();
f_Egamma1=D_Minus_Pi[i].child(0).child(0).e();
f_Egamma2=D_Minus_Pi[i].child(0).child(1).e();

double theta1,theta2;


Hep3Vector The1h(D_Minus_Pi[i].child(0).child(0).mdstGamma().px(),D_Minus_Pi[i].child(0).child(0).mdstGamma().py(),D_Minus_Pi[i].child(0).child(0).mdstGamma().pz());
theta1=The1h.theta()*180/3.14159265359;
f_Gam1thet=theta1;

Hep3Vector The2h(D_Minus_Pi[i].child(0).child(1).mdstGamma().px(),D_Minus_Pi[i].child(0).child(1).mdstGamma().py(),D_Minus_Pi[i].child(0).child(1).mdstGamma().pz());
theta2=The2h.theta()*180/3.14159265359;
f_Gam2thet=theta2;

 if(f_dzero >1.68  && f_dzero < 2.06){
 if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
 if(f_Pizmom > 1.06 ){
 if(f_Pimom > 0.84 ){
 if(f_PD > 2.5){

//Hard Photon 1
 if(f_Gam1thet < -60 && f_Egamma1 > 0.150){photon1cutflag=1;}  //BWD
 else if(f_Gam1thet > 73 && f_Egamma1 > 0.100){photon1cutflag=1;}  //FWD
 else if(f_Gam1thet > -60 && f_Gam1thet < 73 && f_Egamma1 > 0.050){photon1cutflag=1;}  //BARREL
 else{photon1cutflag=0;}

//Hard Photon 2
 if(f_Gam2thet < -60 && f_Egamma2 > 0.150){photon2cutflag=1;}  //BWD
 else if(f_Gam2thet > 73 && f_Egamma2 > 0.100){photon2cutflag=1;}  //FWD
 else if(f_Gam2thet > -60 && f_Gam2thet < 73 && f_Egamma2 > 0.050){photon2cutflag=1;}  //BARREL
else{photon2cutflag=0;}

//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){

GoodCand_M[i]=1;

}//Photon cuts


}}}}}

}//D_Minus_Pi Loop


 
//BCS for D candidate
double PiZMass, PiZMass_FOM, PiZMass_FOM_best=100000;
int PorM_Flag_2=0;
GoodPGuys=0; GoodMGuys=0;

 TrueCand_Evt_D=0;


//D+ loop
      for (int i=0;i<D_Plus_Pi.size();i++) {
//See if any true D* candidate is there in the event
      setMCtruth(D_Plus_Pi[i]);
      if(getMCtruthFlag(D_Plus_Pi[i])==1 && GoodCand_P[i]==1){TrueCand_Evt_D=1;}


if(GoodCand_P[i]==1){GoodPGuys++;
PiZMass= D_Plus_Pi[i].child(0).mdstPi0().mass();
PiZMass_FOM = fabs(PiZMass-0.1349766);
//Check best
if(PiZMass_FOM < PiZMass_FOM_best){
PiZMass_FOM_best= PiZMass_FOM;
i_best=i;
PorM_Flag_2=1;
}
}//Goodcand
}//D+ loop


//D- loop
      for (int i=0;i<D_Minus_Pi.size();i++) {
//See if any true D* candidate is there in the event
      setMCtruth(D_Minus_Pi[i]);
      if(getMCtruthFlag(D_Minus_Pi[i])==1 && GoodCand_M[i]==1){TrueCand_Evt_D=1;}


      if(GoodCand_M[i]==1){GoodMGuys++;
      PiZMass= D_Minus_Pi[i].child(0).mdstPi0().mass();;
      PiZMass_FOM = fabs(PiZMass-0.1349766);
      //Check best
      if(PiZMass_FOM < PiZMass_FOM_best){
      PiZMass_FOM_best= PiZMass_FOM;
      i_best=i;
      PorM_Flag_2=-1;
      }
      }//Goodcand
      }//D- loop

if(GoodPGuys > 0 || GoodMGuys > 0){
d_mult_after->column("PMult",GoodPGuys);
d_mult_after->column("NMult",GoodMGuys);
d_mult_after->column("TMult",GoodPGuys+GoodMGuys);
d_mult_after->dumpData();
}

mult_flag_D=GoodPGuys+GoodMGuys;


//If there are multiple best candidates
// M side double check
if(PorM_Flag_2 == -1){
      for (int i=0;i<D_Minus_Pi.size();i++) {
      if(GoodCand_M[i]==1){
      if(D_Minus_Pi[i].child(0).mdstPi0().mass() == D_Minus_Pi[i_best].child(0).mdstPi0().mass()){//If the pi0 is common
	double dr1,dr2;
        HepVector ip_par1(param_at_ip(D_Minus_Pi[i_best].child(1)));
        HepVector ip_par2(param_at_ip(D_Minus_Pi[i].child(1)));

	dr1 = ip_par1[0];//current best candidate
	dr2 = ip_par2[0];

	if(dr2 < dr1){i_best = i;}
cout<<"BCS switch happened"<<endl;
	}//If the pi0 is common
}}
}

// P side double check
 if(PorM_Flag_2 == 1){
 for (int i=0;i<D_Plus_Pi.size();i++) {
 if(GoodCand_P[i]==1){
 if(D_Plus_Pi[i].child(0).mdstPi0().mass() == D_Plus_Pi[i_best].child(0).mdstPi0().mass()){//If the pi0 is common
 double dr1,dr2;
 HepVector ip_par1(param_at_ip(D_Plus_Pi[i_best].child(1)));
 HepVector ip_par2(param_at_ip(D_Plus_Pi[i].child(1)));

 dr1 = ip_par1[0];//current best candidate
 dr2 = ip_par2[0];
 if(dr2 < dr1){i_best = i;}
cout<<"BCS switch happened"<<endl;
  }//If the pi0 is common
  }}
  }
                                                                           


//if(PorM_Flag_2 == 1){FillTuple_big(D_Plus_Pi[i_best],nt_dpipi_big);}
//else if(PorM_Flag_2 == -1){FillTuple_big(D_Minus_Pi[i_best],nt_dpipi_big);}

/*
      for (int i=0;i<D_Plus_Pi.size();i++) {
      FillTuple_big(D_Plus_Pi[i],nt_dpipi_big);
          *status = 1;
          }
            for (int i=0;i<D_Minus_Pi.size();i++) {
            FillTuple_big(D_Minus_Pi[i],nt_dpipi_big);
                *status = 1;
                }
*/

}

}





void dtohpi::FillTuple_small_D(Particle& Dst_cand, BelleTuple *nt){




//Fill stuff
    Belle_event_Manager& EvtMgr=Belle_event_Manager::get_manager();
    Belle_event& Evt = *EvtMgr.begin();
    int ExpNo = Evt.ExpNo();
    int RunNo = Evt.RunNo();
    int EvtNo = Evt.EvtNo();

    nt->column("Exp",(float)ExpNo);
    nt->column("Run",(float)RunNo);
    nt->column("Evt",(float)EvtNo);


double massdiff, E_Asy;

      nt->column("DStarMass",Dst_cand.mass());
      nt->column("DCharge",Dst_cand.child(0).charge());
      massdiff=Dst_cand.mass()-Dst_cand.child(0).mass();
      nt->column("DeltaM",massdiff);
      nt->column("DMass",Dst_cand.child(0).mass());

      nt->column("PstDst",pStar(Dst_cand).vect().mag());
      nt->column("CThetDst",pStar(Dst_cand).cosTheta());
      nt->column("PstD",pStar(Dst_cand.child(0)).vect().mag());
      nt->column("PiPID",kaonId(Dst_cand.child(0).child(1)));
      nt->column("EGamma1",Dst_cand.child(0).child(0).child(0).e());
      nt->column("EGamma2",Dst_cand.child(0).child(0).child(1).e());
      nt->column("EGam1s",Dst_cand.child(1).child(0).e());
      nt->column("EGam2s",Dst_cand.child(1).child(1).e());

	double dr, dz;
        HepVector ip_par(param_at_ip(Dst_cand.child(0).child(1)));
	dr=ip_par[0];
        dz=ip_par[3];
      nt->column("PiDR",dr);
      nt->column("PiDZ",dz);





E_Asy = ( Dst_cand.child(0).child(0).child(0).e() - Dst_cand.child(0).child(0).child(1).e() ) / ( Dst_cand.child(0).child(0).child(0).e() + Dst_cand.child(0).child(0).child(1).e() );

      nt->column("AsyEGam",E_Asy);
      nt->column("PiZMass",Dst_cand.child(0).child(0).mdstPi0().mass());
      nt->column("PiZsMass",Dst_cand.child(1).mdstPi0().mass());
      nt->column("PiZMom",Dst_cand.child(0).child(0).ptot());
      nt->column("PiMom",Dst_cand.child(0).child(1).ptot());
      nt->column("PiCThet",Dst_cand.child(0).child(1).p().cosTheta());
      nt->column("PiZsMom",Dst_cand.child(1).ptot());

      nt->column("MultFlag",mult_flag);
      nt->column("TruEvt",TrueCand_Evt);



//Helicity
HepLorentzVector dstar4V, dzero4V, piz4V;
Hep3Vector dzero_boost; double cos_helicity;

dstar4V=Dst_cand.p();
dzero4V=Dst_cand.child(0).p();
piz4V=Dst_cand.child(0).child(0).p();
dzero_boost=-1*dzero4V.boostVector();
dstar4V.boost(dzero_boost);
piz4V.boost(dzero_boost);
cos_helicity=cos(dstar4V.angle(piz4V));

      nt->column("Heli",cos_helicity);

//hard-pi0 Photon thetas
double theta1, theta2;

Hep3Vector The1h(Dst_cand.child(0).child(0).child(0).mdstGamma().px(), Dst_cand.child(0).child(0).child(0).mdstGamma().py(), Dst_cand.child(0).child(0).child(0).mdstGamma().pz());
theta1=The1h.theta()*180/3.14159265359;
      nt->column("Gam1hThe",theta1);

Hep3Vector The2h(Dst_cand.child(0).child(0).child(1).mdstGamma().px(), Dst_cand.child(0).child(0).child(1).mdstGamma().py(), Dst_cand.child(0).child(0).child(1).mdstGamma().pz());
theta2=The2h.theta()*180/3.14159265359;
      nt->column("Gam2hThe",theta2);

//Gamma Helicity^Z
HepLorentzVector dst4V, gam4V, dst4V_st;
Hep3Vector piz_boost, dst_boost; 

dst4V=Dst_cand.p();             //D*  mom
piz4V=Dst_cand.child(1).p();    //pi0 mom
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
//soft-pi0 Photon thetas
int theta1_flag, theta2_flag, category;


Hep3Vector The1s(Dst_cand.child(1).child(0).mdstGamma().px(),Dst_cand.child(1).child(0).mdstGamma().py(),Dst_cand.child(1).child(0).mdstGamma().pz());
theta1=The1s.theta()*180/3.14159265359;
      nt->column("Gam1sThe",theta1);

Hep3Vector The2s(Dst_cand.child(1).child(1).mdstGamma().px(),Dst_cand.child(1).child(1).mdstGamma().py(),Dst_cand.child(1).child(1).mdstGamma().pz());
theta2=The2s.theta()*180/3.14159265359;
      nt->column("Gam2sThe",theta2);


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
      nt->column("PiZ0ID",IDhep(Dst_cand.child(1)));
      nt->column("PiZ1ID",IDhep(Dst_cand.child(0).child(0)));
      nt->column("PiID",IDhep(Dst_cand.child(0).child(1)));

      nt->column("DstF",getMCtruthFlag(Dst_cand));
      nt->column("DF",getMCtruthFlag(Dst_cand.child(0)));
      nt->column("PiZ0F",getMCtruthFlag(Dst_cand.child(1)));
      nt->column("PiZ1F",getMCtruthFlag(Dst_cand.child(0).child(0)));
      nt->column("PiF",getMCtruthFlag(Dst_cand.child(0).child(1)));

double gen_momt1=-1, gen_momt2=-1;


if(Dst_cand.child(0).child(1).genHepevt()){
if(Dst_cand.child(0).child(0).genHepevt()){

Hep3Vector Gen_Pi_Mom(Dst_cand.child(0).child(1).genHepevt().PX(),Dst_cand.child(0).child(1).genHepevt().PY(),Dst_cand.child(0).child(1).genHepevt().PZ());

Hep3Vector Gen_Piz_Mom(Dst_cand.child(0).child(0).genHepevt().PX(),Dst_cand.child(0).child(0).genHepevt().PY(),Dst_cand.child(0).child(0).genHepevt().PZ());

gen_momt1=Gen_Pi_Mom.mag();
gen_momt2=Gen_Piz_Mom.mag();
}}

      nt->column("GenPiMom",gen_momt1);
      nt->column("GenPiZMom",gen_momt2);


	//Looking for SCF
        Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
if(Dst_cand.child(0).genHepevt()){

	Gen_hepevt& D_child1 = GenMgr(Panther_ID(Dst_cand.child(0).genHepevt().da(0)));
	Gen_hepevt& D_child2 = GenMgr(Panther_ID(Dst_cand.child(0).genHepevt().da(1)));

if(D_child1 && D_child2){
      	nt->column("GenPizID",D_child1.idhep());
      	nt->column("GenPiID",D_child2.idhep());
}
else{
           nt->column("GenPizID",-1);
           nt->column("GenPiID",-1);
}

}
else{
           nt->column("GenPizID",-1);
           nt->column("GenPiID",-1);
}


}

      nt->dumpData();

}



void dtohpi::FillTuple_big(Particle& D_cand, BelleTuple *nt){

//Fill stuff

    Belle_event_Manager& EvtMgr=Belle_event_Manager::get_manager();
    Belle_event& Evt = *EvtMgr.begin();
    int ExpNo = Evt.ExpNo();
    int RunNo = Evt.RunNo();
    int EvtNo = Evt.EvtNo();

    nt->column("Exp",(float)ExpNo);
    nt->column("Run",(float)RunNo);
    nt->column("Evt",(float)EvtNo);

    nt->column("DCharge",D_cand.charge());
    nt->column("DMass",D_cand.mass());
    nt->column("PstD",pStar(D_cand).vect().mag());
    nt->column("PiPID",kaonId(D_cand.child(1)));
    nt->column("EGamma1",D_cand.child(0).child(0).e());
    nt->column("EGamma2",D_cand.child(0).child(1).e());
    double dr, dz;
    HepVector ip_par(param_at_ip(D_cand.child(1)));
    dr=ip_par[0];
    dz=ip_par[3];
    nt->column("PiDR",dr);
    nt->column("PiDZ",dz);
    nt->column("PiZMass",D_cand.child(0).mdstPi0().mass());
    nt->column("PiZMom",D_cand.child(0).ptot());
    nt->column("PiMom",D_cand.child(1).ptot());

      nt->column("MultFlag",mult_flag_D);
      nt->column("TruEvt",TrueCand_Evt_D);

//hard-pi0 Photon thetas
double theta1, theta2;

Hep3Vector The1h(D_cand.child(0).child(0).mdstGamma().px(), D_cand.child(0).child(0).mdstGamma().py(), D_cand.child(0).child(0).mdstGamma().pz());
theta1=The1h.theta()*180/3.14159265359;
      nt->column("Gam1hThe",theta1);

      Hep3Vector The2h(D_cand.child(0).child(1).mdstGamma().px(), D_cand.child(0).child(1).mdstGamma().py(), D_cand.child(0).child(1).mdstGamma().pz());
      theta2=The2h.theta()*180/3.14159265359;
            nt->column("Gam2hThe",theta2);




if(dataType){


      setMCtruth(D_cand);
      nt->column("DID",IDhep(D_cand));
      nt->column("PiZID",IDhep(D_cand.child(0)));
      nt->column("PiID",IDhep(D_cand.child(1)));
      nt->column("DF",getMCtruthFlag(D_cand));
      nt->column("PiZF",getMCtruthFlag(D_cand.child(0)));
      nt->column("PiF",getMCtruthFlag(D_cand.child(1)));

        //Looking for SCF
        Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
if (D_cand.genHepevt())
{
           Gen_hepevt& D_child1 = GenMgr(Panther_ID(D_cand.genHepevt().da(0)));
           Gen_hepevt& D_child2 = GenMgr(Panther_ID(D_cand.genHepevt().da(1)));
if(D_child1 && D_child2){
           nt->column("GenPizID",D_child1.idhep());
           nt->column("GenPiID",D_child2.idhep());
}
else{
           nt->column("GenPizID",-1);
           nt->column("GenPiID",-1);
}

}
else{
           nt->column("GenPizID",-1);
           nt->column("GenPiID",-1);
}


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
