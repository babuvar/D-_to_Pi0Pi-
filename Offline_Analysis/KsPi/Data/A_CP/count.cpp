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
  RooRealVar dzero("dzero","dzero",1.80,1.94,"GeV");
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(dzero));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(dzero));


void count()
{


//  TFile *file = new TFile("KsPi_data_4s.root");
//  TFile *file = new TFile("KsPi_data_5s.root");
//  TFile *file = new TFile("KsPi_data_4s_continuum.root");
//  TTree *t1 = (TTree*)file->Get("h1");
//  Int_t nevt = (Int_t)t1->GetEntries();

//LOAD DATA FILE
  TChain* chain=new TChain("h1");
  chain->Add("KsPi_data_4s.root");
  chain->Add("KsPi_data_5s.root");
  chain->Add("KsPi_data_4s_continuum.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;


 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Ksmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;

cout<<"done 1"<<endl;

  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Ksmom",&f_Ksmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);





int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {


      h1->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.84  && f_dzero < 1.90){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
//if(f_Pizmass > 0.119  && f_Pizmass < 0.151){	//***********
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
//if(f_Pimom > 0.84 &&  f_Pimom < 1.74){



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



if(f_PD > 2.5 ){
num_pass++;

}



}//Photon cuts

}}}}}//}//}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;





}












