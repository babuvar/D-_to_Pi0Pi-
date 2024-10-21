//To find the next significant P*(D*) bin
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
void PstD_opt2(void)
{


float fom[400], eff[400], cut[400], sig[400]={0}, bkg[400]={0}, e_cut[400]={0}, e_eff[400]={0}, e_fom[400]={0}, siba_MC[400]={0}, siba_data[400]={0}, ratio[400], e_ratio[250]={0}, e_sig[250]={0}, e_bkg[250]={0};  //cut in terms of sigma
float cutoff=5.0;
//float cutoff=2.95;
//float cutoff=2.5;

float width=0.01;
float numbins=(cutoff -2.5)/width;
cout<<"width = "<<width<<endl;

int stopbin;
float fullarea=0;


for(int i=0;i<numbins;i++){cut[i]=2.5+(i*width);}// cout<<cut[i]<<endl;}


  //LOAD DATA FILE
  TChain* chain=new TChain("h2");

chain->Add("pipi_MC_4s.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

chain->Add("pipi_data_4s.root");
  Int_t nevt2=(int)chain->GetEntries();



 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Df, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass;



  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Df",&f_Df);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);

int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
    {


      chain->GetEntry(i);
stopbin=floor((f_PD-2.5)/width);





if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5 && f_PD < cutoff){

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


if(f_dzero >1.80  && f_dzero < 1.90){       //~3sigma range to estimate F.O.M.

if(f_Df == 1 || f_Df == 10){
fullarea++;
for(int j=0;j<stopbin;j++){sig[j]++; } 
}
else{
for(int j=0;j<stopbin;j++){bkg[j]++; }
}

}       //~3sigma range to estimate F.O.M.

//Fill sideband events
if( (f_dzero >1.75  && f_dzero < 1.80) || (f_dzero >1.90  && f_dzero < 1.95) ){//MD-sideband
//if( f_dzero >1.92  && f_dzero < 2.06 ){//MD-sideband
for(int j=0;j<stopbin;j++){
siba_MC[j]++; }
}//MD-sideband


}


}}}}
    }


 //Fill sideband events for Data
  for(int i=nevt;i<nevt2;i++)
//   for(int i=0;i<50000;i++)
   {

      chain->GetEntry(i);
stopbin=floor((f_PD-2.5)/width);


if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5 && f_PD < cutoff){

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



//Fill sideband events
if( (f_dzero >1.75  && f_dzero < 1.80) || (f_dzero >1.90  && f_dzero < 1.95) ){//MD-sideband
//if( f_dzero >1.92  && f_dzero < 2.06 ){//MD-sideband
for(int j=0;j<stopbin;j++){
siba_data[j]++; }
}//MD-sideband


}


}}}}

    }

float bestfom=0, bestcut_sigma, bestcut, besti; 


//Calculating FOM, Efficiency and data-MC ratio
for(int i=0;i<numbins;i++){
eff[i] =sig[i]/fullarea;
if(siba_MC[i] != 0){ratio[i] = (siba_data[i])/siba_MC[i];} else{ratio[i] =1;}

sig[i]=sig[i]*0.477;
//sig[i]=sig[i]*0.69;

bkg[i]=bkg[i]*ratio[i];
fom[i] = sig[i]/sqrt(sig[i]+bkg[i]); //cout<<sig[i]<<endl;

if(fom[i] > bestfom){bestfom =fom[i]; besti=i;}

}




//Plotting
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;
      	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;
      	TCanvas* cnv3 = new TCanvas("cnv3","cnv3") ;
      	TCanvas* cnv4 = new TCanvas("cnv4","cnv4") ;
	cnv1->Divide(2,2);

   TGraphErrors *gr1 = new TGraphErrors(numbins,cut,eff,e_cut,e_eff);
//   TGraphErrors *gr1 = new TGraphErrors(200,cut,eff,e_cut,e_eff);
   gr1->SetTitle("Signal Efficiency");
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(0.7);
   gr1->SetMarkerColor(kBlue);
   cnv1->cd(3); gr1->Draw("AP");



   TGraphErrors *gr2 = new TGraphErrors(numbins,cut,ratio,e_cut,e_ratio);
//   TGraphErrors *gr2 = new TGraphErrors(200,cut,ratio,e_cut,e_ratio);
   gr2->SetTitle("Data/MC Background Ratio");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(0.7);
   gr2->SetMarkerColor(kBlue);
   cnv1->cd(2); gr2->Draw("AP");


   TGraphErrors *gr2p = new TGraphErrors(50,cut,ratio,e_cut,e_ratio);
   gr2p->SetTitle("Data/MC Background Ratio");
   gr2p->SetMarkerStyle(20);
   gr2p->SetMarkerSize(0.7);
   gr2p->SetMarkerColor(kBlue);
   cnv3->cd(); gr2p->Draw("AP");


   TGraphErrors *gr3 = new TGraphErrors(numbins,cut,fom,e_cut,e_fom);
//   TGraphErrors *gr3 = new TGraphErrors(200,cut,fom,e_cut,e_fom);
   gr3->SetTitle("Figure of Merit");
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(0.7);
   gr3->SetMarkerColor(kBlue);
   cnv1->cd(4); gr3->Draw("AP");

   TGraphErrors *gr3p = new TGraphErrors(50,cut,fom,e_cut,e_fom);
   gr3p->SetTitle("Figure of Merit");
   gr3p->SetMarkerStyle(20);
   gr3p->SetMarkerSize(0.7);
   gr3p->SetMarkerColor(kBlue);
   cnv4->cd(); gr3p->Draw("AP");




   TGraphErrors *gr4 = new TGraphErrors(numbins,cut,sig,e_cut,e_sig);
//   TGraphErrors *gr4 = new TGraphErrors(200,cut,sig,e_cut,e_sig);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(0.7);
   gr4->SetMarkerColor(kBlue);

   TGraphErrors *gr4p = new TGraphErrors(50,cut,sig,e_cut,e_sig);
   gr4p->SetMarkerStyle(20);
   gr4p->SetMarkerSize(0.7);
   gr4p->SetMarkerColor(kBlue);
   cnv2->cd(); gr4p->Draw("AP");

   TGraphErrors *gr5 = new TGraphErrors(numbins,cut,bkg,e_cut,e_bkg);
//   TGraphErrors *gr5 = new TGraphErrors(200,cut,bkg,e_cut,e_bkg);
   gr5->SetMarkerStyle(20);
   gr5->SetMarkerSize(0.7);
   gr5->SetMarkerColor(kRed);


   TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Signal and Background levels");
     mg->Add(gr4);
     mg->Add(gr5);

   cnv1->cd(1); mg->Draw("AP");

//      	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;
//   	cnv2->cd(4); gr3->Draw("AP");


cout<<"best f.o.m. ="<<bestfom<<endl;
cout<<"best cut = +/- "<<cut[besti]<<endl;





}



