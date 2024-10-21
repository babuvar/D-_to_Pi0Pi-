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
using namespace RooFit;
void EGams_opt(void)
{


//float fom[250], eff[250], cut[250], sig[250]={0}, bkg[250]={0}, e_cut[250]={0}, e_eff[250], e_fom[250], siba_MC[250]={0}, siba_data[250]={0}, ratio[250], e_ratio[250];  //cut in terms of sigma

float fom[250], eff[250], cut[250], sig[250]={0}, bkg[250]={0}, e_cut[250]={0}, e_eff[250]={0}, e_fom[250]={0}, siba_MC[250]={0}, siba_data[250]={0}, ratio[250], e_ratio[250]={0}, e_sig[250]={0}, e_bkg[250]={0};  //cut in terms of sigma

float numbins=200;
float width=(0.430 -0.030)/numbins;
//float width=(0.230 -0.030)/numbins;

int stopbin, stopbin2, stopbin_chosen;
float fullarea=0;

//for(int i=0;i<numbins;i++){cut[i]=0.030+(i*width);}
for(int i=0;i<numbins;i++){cut[i]=0.030+(i*width);}// cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}


  //LOAD DATA FILE
  TChain* chain=new TChain("h1");



chain->Add("pipi_MC6.root");
chain->Add("pipi_MC7.root");

  Int_t nevt=(int)chain->GetEntries();

chain->Add("pipi_data.root");

  Int_t nevt2=(int)chain->GetEntries();
//  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Pizmass, f_Dstf, f_Df, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Gam1thet, f_Gam2thet;
int photon1cutflag=0, photon2cutflag=0;

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
  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);

int photon1cutflag=0, photon2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
    {

      chain->GetEntry(i);
stopbin=floor((f_Egam1s-0.03)/width);
stopbin2=floor((f_Egam2s-0.03)/width);
if(stopbin < stopbin2){stopbin_chosen=stopbin;}
else{stopbin_chosen=stopbin2;}

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

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Egam1s > 0.03  && f_Egam1s <  0.43){
if(f_Egam2s > 0.03  && f_Egam2s <  0.43){
//if(f_Categ == 1){
if(f_Categ == 5){
//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){

if(f_dzero >1.80  && f_dzero < 1.90){       //~3sigma M_D range to estimate F.O.M.

if(f_Df == 1){
fullarea++;
for(int j=0;j<stopbin_chosen;j++){sig[j]++; } 
//for(int j=0;j<stopbin;j++){sig[j]++; }
}
else{
for(int j=0;j<stopbin_chosen;j++){bkg[j]++; }
//for(int j=0;j<stopbin;j++){bkg[j]++;} 
}

}//~3sigma M_D range to estimate F.O.M.

//Fill sideband events for MC
if( (f_dzero >1.70  && f_dzero < 1.76) || (f_dzero >1.92  && f_dzero < 2.0) ){
for(int j=0;j<stopbin_chosen;j++){siba_MC[j]++;}
}



}}}}}}}
    }

//Fill sideband events for Data
  for(int i=nevt;i<nevt2;i++) 
//for(int i=0;i<50000;i++)
    {

      chain->GetEntry(i);
stopbin=floor((f_Egam1s-0.03)/width);
stopbin2=floor((f_Egam2s-0.03)/width);
if(stopbin < stopbin2){stopbin_chosen=stopbin;}
else{stopbin_chosen=stopbin2;}

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

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Egam1s > 0.03  && f_Egam1s <  0.43){
if(f_Egam2s > 0.03  && f_Egam2s <  0.43){
//if(f_Categ == 1){
if(f_Categ == 5){
//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){

//Fill sideband events
if( (f_dzero >1.70  && f_dzero < 1.76) || (f_dzero >1.92  && f_dzero < 2.0) ){
for(int j=0;j<stopbin_chosen;j++){siba_data[j]++;}
}



}}}}}}}
    }

//for(int i=0;i<numbins;i++){cout<<"for i = "<<i<<", siba_data ="<<siba_data[i]<<"\t siba_MC ="<<siba_MC[i]<<endl;}



//cout<<"Done : 1"<<endl;



float bestfom=0, bestcut_sigma, bestcut, besti; 


//Calculating FOM, Efficiency and data-MC ratio
for(int i=0;i<numbins;i++){
eff[i] =sig[i]/fullarea;
if(siba_MC[i] != 0){ratio[i] = (siba_data[i]*2)/siba_MC[i];} else{ratio[i] =1;}

sig[i]=sig[i]*0.477*0.5;
bkg[i]=bkg[i]*ratio[i]*0.5;
fom[i] = sig[i]/sqrt(sig[i]+bkg[i]); //cout<<sig[i]<<endl;

if(fom[i] > bestfom){bestfom =fom[i]; besti=i;}

}




//Plotting
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;
	cnv1->Divide(2,2);

//   TGraphErrors *gr1 = new TGraphErrors(numbins,cut,eff,e_cut,e_eff);
   TGraphErrors *gr1 = new TGraphErrors(60,cut,eff,e_cut,e_eff);
   gr1->SetTitle("Signal Efficiency");
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(0.7);
   gr1->SetMarkerColor(kBlue);
   cnv1->cd(3); gr1->Draw("AP");

//   TGraphErrors *gr2 = new TGraphErrors(numbins,cut,ratio,e_cut,e_ratio);
   TGraphErrors *gr2 = new TGraphErrors(60,cut,ratio,e_cut,e_ratio);
   gr2->SetTitle("Data/MC Background Ratio");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(0.7);
   gr2->SetMarkerColor(kBlue);
   cnv1->cd(2); gr2->Draw("AP");


//   TGraphErrors *gr3 = new TGraphErrors(numbins,cut,fom,e_cut,e_fom);
   TGraphErrors *gr3 = new TGraphErrors(60,cut,fom,e_cut,e_fom);
   gr3->SetTitle("Figure of Merit");
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(0.7);
   gr3->SetMarkerColor(kBlue);
   cnv1->cd(4); gr3->Draw("AP");


//   TGraphErrors *gr4 = new TGraphErrors(numbins,cut,sig,e_cut,e_sig);
   TGraphErrors *gr4 = new TGraphErrors(60,cut,sig,e_cut,e_sig);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(0.7);
   gr4->SetMarkerColor(kBlue);


//   TGraphErrors *gr5 = new TGraphErrors(numbins,cut,bkg,e_cut,e_bkg);
   TGraphErrors *gr5 = new TGraphErrors(60,cut,bkg,e_cut,e_bkg);
   gr5->SetMarkerStyle(20);
   gr5->SetMarkerSize(0.7);
   gr5->SetMarkerColor(kRed);


   TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Signal and Background levels");
     mg->Add(gr4);
     mg->Add(gr5);

   cnv1->cd(1); mg->Draw("AP");


      	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;
        cnv2->cd(4); gr3->Draw("AP");



cout<<"best f.o.m. ="<<bestfom<<endl;
cout<<"best cut = +/- "<<cut[besti]<<endl;



}



