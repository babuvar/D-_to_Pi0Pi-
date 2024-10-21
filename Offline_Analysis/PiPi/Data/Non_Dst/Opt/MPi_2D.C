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
void MPi_2D(void)
{


float fom[400][400], eff[400][400], cut[400], sig[400][400]={0}, bkg[400][400]={0}, siba_MC[400][400]={0}, siba_data[400][400]={0}, ratio[400][400];  //cut in terms of sigma
float numbins=32;
float width=(0.151 - 0.119)/numbins;


cout<<"width = "<<width<<endl;


int stopbin, stopbin;
float fullarea=0;



for(int i=0;i<numbins;i++){cut[i]=0.119+(i*width); }

  //LOAD DATA FILE
  TChain* chain=new TChain("h2");


chain->Add("pipi_MC_4s.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

chain->Add("pipi_data_4s.root");
  Int_t nevt2=(int)chain->GetEntries();



cout<<"DONE 1"<<endl;
 Float_t Deltam, f_PD, f_dzero, f_Pizmass, f_Dstf, f_Df, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet,
 f_Gam2thet, f_Pizmass;



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



  for(int i=0; i < nevt ;i++) 
//  for(int i=0;i<100000;i++)
//    for(int i=0;i<1000;i++)
   {


      chain->GetEntry(i);
stopbin=floor((f_Pizmass-0.119)/width);



if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06){ 
if(f_Pimom > 1.31){
if(f_PD > 2.72){


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




//cout<<"f_Pizmom = "<<f_Pizmom<<endl;
//cout<<"stopbin1 = "<<stopbin1<<endl;


if(f_dzero >1.80  && f_dzero < 1.90){       //~3sigma range to estimate F.O.M.
if(f_Df != 1){
for(int j=0;j<stopbin;j++){
for(int k=stopbin;k<numbins;k++){bkg[j][k]++; }}

}

}       //~3sigma range to estimate F.O.M.

//Fill sideband events
if( (f_dzero >1.70  && f_dzero < 1.76) || (f_dzero >1.92  && f_dzero < 2.0) ){//MD-sideband
for(int j=0;j<stopbin;j++){
for(int k=stopbin;k<numbins;k++){siba_MC[j][k]++; }}
}//MD-sideband

}}}}

}
    }

//Generate Signal
TRandom *MyRandom= new TRandom(12345);
Float_t xrand;
for (Int_t i = 0; i < 326280; i++) {//10 x sig
//xrand = MyRandom->Gaus(0.13498,0.005761);
xrand = MyRandom->Gaus(0.13498,0.005918);
if(xrand > 0.119  && xrand < 0.151){
stopbin=floor((xrand-0.119)/width);
for(int j=0;j<stopbin;j++){
for(int k=stopbin;k<numbins;k++){sig[j][k]++; }}
}
}
fullarea=326280;//10 x sig



cout<<"DONE 2"<<endl;

//------------------------------DATA-------------------------------

 //Fill sideband events for Data
  for(int i=nevt;i<nevt2;i++)
//   for(int i=0;i<50000;i++)
   {
      chain->GetEntry(i);
stopbin=floor((f_Pizmass-0.119)/width);


if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06){ 
if(f_Pimom > 0.84){
if(f_PD > 2.72){

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
if( (f_dzero >1.70  && f_dzero < 1.76) || (f_dzero >1.92  && f_dzero < 2.0) ){//MD-sideband
for(int j=0;j<stopbin;j++){
for(int k=stopbin;k<numbins;k++){siba_data[j][k]++; }}
}//MD-sideband

}}}}
}

    }

cout<<"DONE 3"<<endl;


float bestfom=0, bestcut_sigma, bestcut, besti, bestj; 

TGraph2D *gr1 = new TGraph2D(); TGraph2D *gr2 = new TGraph2D();
TGraph2D *gr3 = new TGraph2D(); 

int count=-1; float x,y,z;

//Calculating FOM and Efficiency
for(int i=0;i<numbins;i++){
for(int j=i+1;j<numbins;j++){
eff[i][j] =sig[i][j]/fullarea;
if(siba_MC[i][j] != 0){ratio[i][j] = (siba_data[i][j]*2)/siba_MC[i][j];} else{ratio[i][j] =1;}

sig[i][j]=sig[i][j]*0.477*0.5*0.1;
bkg[i][j]=bkg[i][j]*ratio[i][j]*0.5;


fom[i][j] = sig[i][j]/sqrt(sig[i][j]+bkg[i][j]); //cout<<sig[i]<<endl;
 if(fom[i][j] > bestfom){bestfom =fom[i][j]; besti=i; bestj=j;}

}}
// cout<<"-------------------------DONE-----------------------"<<endl;

count=;
//Plotting
for(int i=0;i<numbins;i++){
for(int j=0;j<numbins;j++){
//for(int i=0;i<59;i++){
//for(int j=0;j<59;j++){
 count++;
 x=cut[i];
 y=cut[j];
 z=fom[i][j];
 gr1->SetPoint(count,x,y,z);
 z=eff[i][j];
 gr2->SetPoint(count,x,y,z);
 z=ratio[i][j];
 gr3->SetPoint(count,x,y,z);

}}




      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;
	cnv1->Divide(2,2);

   gr1->SetTitle("Figure of Merit");
//   cnv1->cd(1); gr1->Draw("p");
   cnv1->cd(1); gr1->Draw("colz0");


   gr2->SetTitle("Efficiency");
//   cnv1->cd(2); gr2->Draw("p");
   cnv1->cd(2); gr2->Draw("colz0");



   gr3->SetTitle("Data/MC Background Ratio");
//   cnv1->cd(3); gr3->Draw("p");
   cnv1->cd(3); gr3->Draw("colz0");



      	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;
   	cnv2->cd(1); gr1->Draw("colz0");


cout<<"best f.o.m. ="<<bestfom<<endl;
cout<<"best down cut =  "<<cut[besti]<<"\t best up cut =  "<<cut[bestj]<<endl;
//cout<<"best cut = +/- "<<cut1[besti]*0.0055<<" ("<<cut1[besti]<<" sigma cut)"<<endl;




}



