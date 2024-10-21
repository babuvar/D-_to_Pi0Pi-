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
#include "TColor.h"
#include "TAttFill.h"
using namespace RooFit;
void Purity_Bin(void)
{


float  cut[300], sig[300]={0}, bkg[300]={0}, tot[300]={0}, e_sig[300], e_bkg[300], e_cut[300]={0}, purity[300], e_purity[300];  //cut in terms of sigma
float numbins=25;
float width=(5.0 - 2.5)/numbins;
cout<<"width = "<<width<<endl;

int stopbin;
float fullarea=0;


for(int i=0;i<numbins;i++){cut[i]=2.5+((i+0.5)*width);}// cout<<cut[i]<<endl;}


  //LOAD DATA FILE
  TChain* chain=new TChain("h1");



//chain->Add("pipi_MC.root");
//chain->Add("pipi_MC2.root");
//chain->Add("pipi_MC4.root");//pipi_MC4.root has had best candidate selection and Mom_pi+ > 0.5 (instead of 0.75)

chain->Add("pipi_MC7.root");




  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);

int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
    {


      chain->GetEntry(i);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.90){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizmom > 0.95 ){
if(f_Pimom > 0.74 ){
if(f_PD > 2.5 ){

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
if(f_Categ == 5 && f_Egam1s > 0.044){sphoton1cutflag=1;}
else if(f_Categ == 2 && f_Egam1s > 0.066){sphoton1cutflag=1;}
else if(f_Categ == 6 && f_Egam1s > 0.044){sphoton1cutflag=1;}
else if(f_Categ == 1 && f_Egam1s > 0.050){sphoton1cutflag=1;}else{sphoton1cutflag=0;}

//Soft Photon 2 
if(f_Categ == 5 && f_Egam2s > 0.044){sphoton2cutflag=1;}
else if(f_Categ == 2 && f_Egam2s > 0.036){sphoton2cutflag=1;}
else if(f_Categ == 6 && f_Egam2s > 0.054){sphoton2cutflag=1;}
else if(f_Categ == 1 && f_Egam2s > 0.050){sphoton2cutflag=1;}else{sphoton2cutflag=0;}


//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){
if(sphoton1cutflag == 1 && sphoton2cutflag == 1){


stopbin=floor((f_PD-2.5)/width);


if(f_Dstf == 1){
fullarea++;
sig[stopbin]++; 

}
else{
bkg[stopbin]++; 
}

}}}}}}}}
    }


//for(int i=0;i<numbins;i++){cout<<"sig["<<i<<"] = "<<sig[i]<<"\t bkg["<<i<<"] = "<<bkg[i]<<"\t cut["<<i<<"] = "<<cut[i]<<endl;}




//cout<<"Done : 1"<<endl;



float bestfom=0, bestcut_sigma, bestcut, besti; 


//Calculating FOM and Efficiency
for(int i=0;i<numbins;i++){
sig[i]=sig[i]*0.477;
e_sig[i]=sqrt(sig[i]);
e_bkg[i]=sqrt(bkg[i]);

tot[i]=sig[i]+bkg[i];
purity[i] = sig[i]/tot[i];
e_purity[i] = purity[i]*sqrt((1/sig[i])+(1/tot[i]));
}



//Plotting
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;


   TGraphErrors *gr1 = new TGraphErrors(numbins,cut,sig,e_cut,e_sig);
   TGraphErrors *gr2 = new TGraphErrors(numbins,cut,bkg,e_cut,e_bkg);


   	gr1->SetFillColor(kBlue);
   	gr2->SetFillColor(kRed);

     TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Signal-Background overlay");
     mg->Add(gr2);
     mg->Add(gr1);

   cnv1->cd(); 
	mg->Draw("AB");
//		     	gr2->Draw("AB");
//			gr1->Draw("AB, same");


     	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;

   TGraphErrors *gr3 = new TGraphErrors(numbins,cut,purity,e_cut,e_purity);
   	gr3->SetFillColor(kGreen);
   	gr3->SetTitle("Purity");

   cnv2->cd(); gr3->Draw("AB");




}



