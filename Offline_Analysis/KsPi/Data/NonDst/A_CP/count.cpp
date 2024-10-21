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

TH1D* data_p=new TH1D("pdat", "pdat", 100, 1.8, 1.94);
TH1D* data_n=new TH1D("ndat", "ndat", 100, 1.8, 1.94);

void count(void)
{

//LOAD DATA FILE
  TChain* chain=new TChain("h2");
  chain->Add("kspi_4sd.root");
  chain->Add("kspi_5sd.root");
  chain->Add("kspi_contd.root");

//  chain->Add("kspi_4smc.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Ksmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;

cout<<"done 1"<<endl;


  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Ksmom",&f_Ksmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);



int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
    {
      h2->GetEntry(i);
		  dzero.setVal(f_dzero);

if(f_dzero >1.84  && f_dzero < 1.90){       //~3sigma range to estimate F.O.M.
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
if(f_PD > 2.65){


num_pass++;


}}}}

    }

cout<<"num_pass"<<num_pass<<endl;


}




