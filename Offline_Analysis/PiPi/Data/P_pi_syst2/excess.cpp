#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TChain.h"
using namespace RooFit;


int bin_theta(float theta);
int bin_mom(float mom);


void excess(void)
{
float nsig_kspi[15][15]={0}, nsig_pipi[15][15]={0}, excess[15][15]={0};
int nsig_kspi_tot=0, nsig_pipi_tot=0;

int binx, biny;

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");



  chain->Add("pipi_MC0_Y4s.root");
  chain->Add("pipi_MC1_Y4s.root");
  chain->Add("pipi_SMC_2M_4s.root");


  Int_t nevt=(int)chain->GetEntries();
//  cout<<"nevt\t"<<nevt <<endl;

  chain->Add("KsPi_GMC0_4s.root");
  chain->Add("KsPi_GMC1_4s.root");
Int_t nevt2=(int)chain->GetEntries();



 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df, f_Ksmom, f_Phimom, f_Picthet;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);
  h1->SetBranchAddress("Picthet",&f_Picthet);

  h1->SetBranchAddress("Ksmom",&f_Ksmom);
  h1->SetBranchAddress("Phimom",&f_Phimom);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);


if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
//if(f_Pimom > 0.84){
if(f_Pimom > 0.85 && f_Pimom < 1.82){
if(f_PD > 2.95){
if(f_Df == 1){
//if(f_Dstf == 1){



//if(f_Categ == 5){
//if(f_Categ != 5){


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

binx=bin_theta(f_Picthet);
biny=bin_mom(f_Pimom);
//cout<<"binx = "<<binx<<"biny = "<<biny<<endl;
nsig_pipi[binx][biny]++;
nsig_pipi_tot++;


}}//Photon cuts

}}}}}}}}//}
    }




  for(int i=nevt;i<nevt2;i++) 
    {


      chain->GetEntry(i);


if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Ksmom > 1.06 ){				//***********
//if(f_Pimom > 0.84 ){
if(f_Pimom > 0.85 && f_Pimom < 1.82){
if(f_PD > 2.95){
if(f_Df == 1){

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
if(sphoton1cutflag == 1 && sphoton2cutflag == 1){


binx=bin_theta(f_Picthet);
biny=bin_mom(f_Pimom);
//cout<<"binx = "<<binx<<"biny = "<<biny<<endl;
nsig_kspi[binx][biny]++;
nsig_kspi_tot++;


//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}//Photon cuts

}}}}}}}//}

    }
//cout<<"nsig_pipi_tot = "<<nsig_pipi_tot<<endl;
//cout<<"nsig_kspi_tot = "<<nsig_kspi_tot<<endl;


for(int i=1; i<=10; i++){
for(int j=1; j<=6; j++){
//cout<<(nsig_pipi[i][j]*100)/69865<<"\t";
//cout<<(nsig_kspi[i][j]*100)/95187<<"\t";
cout<<((nsig_kspi[i][j]*100)/95187)-((nsig_pipi[i][j]*100)/69865)<<"\t";

}
cout<<endl;}



}




int bin_theta(float theta)
{int bin;
if(theta  < -0.45){bin =1;}
else if(theta >=-0.45  && theta < -0.18 ){bin =2;}
else if(theta >=-0.18  && theta < 0.05 ){bin =3;}
else if(theta >=0.05  && theta < 0.24 ){bin =4;}
else if(theta >= 0.24 && theta < 0.40 ){bin =5;}
else if(theta >= 0.40 && theta < 0.53 ){bin =6;}
else if(theta >= 0.53 && theta < 0.64 ){bin =7;}
else if(theta >= 0.64 && theta < 0.74 ){bin =8;}
else if(theta >= 0.74 && theta < 0.84 ){bin =9;}
else {bin =10;}

return bin;

}




int bin_mom(float mom)
{int bin;
if(mom >= 0.85 && mom < 1.00 ){bin = 1;}
else if(mom >= 1.00 && mom < 1.15 ){bin = 2;}
else if(mom >= 1.15 && mom < 1.30 ){bin = 3;}
else if(mom >= 1.30 && mom < 1.45 ){bin = 4;}
else if(mom >= 1.45 && mom < 1.61 ){bin = 5;}
else if(mom >= 1.61 && mom < 1.82 ){bin = 6;}


return bin;

}













