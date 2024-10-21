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
void MomPiH_opt_2D(void)
{


float fom[400][400], eff[400][400], cut1[400], cut2[400], sig[400][400]={0}, bkg[400][400]={0}, siba_MC[400][400]={0}, siba_data[400][400]={0}, ratio[400][400];  //cut in terms of sigma
float numbins=400;
float width1=(4.75 -1.06)/numbins;
float width2=(4.75 -0.84)/numbins;

cout<<"width1 = "<<width1<<endl;
cout<<"width2 = "<<width2<<endl;

int stopbin1, stopbin2;
float fullarea=0;



for(int i=0;i<numbins;i++){cut1[i]=1.06+(i*width1); cut2[i]=0.84+(i*width2);}


  //LOAD DATA FILE
  TChain* chain=new TChain("h2");


chain->Add("pipi_MC_4s.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

chain->Add("pipi_data_4s.root");
  Int_t nevt2=(int)chain->GetEntries();




cout<<"DONE 1"<<endl;
 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Df, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet,
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
stopbin1=floor((f_Pizmom-1.06)/width1);
stopbin2=floor((f_Pimom-0.84)/width2);




if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 && f_Pizmom < 4.75){
if(f_Pimom > 0.84 && f_Pimom < 4.75){
if(f_PD > 2.62){


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
if(f_Df == 1){
fullarea++;
for(int j=0;j<stopbin1;j++){
for(int k=0;k<stopbin2;k++){sig[j][k]++; }}

}
else{
for(int j=0;j<stopbin1;j++){
for(int k=0;k<stopbin2;k++){bkg[j][k]++; }}

}

}       //~3sigma range to estimate F.O.M.

//Fill sideband events
if( (f_dzero >1.70  && f_dzero < 1.76) || (f_dzero >1.92  && f_dzero < 2.0) ){//MD-sideband
for(int j=0;j<stopbin1;j++){
for(int k=0;k<stopbin2;k++){siba_MC[j][k]++; }}
}//MD-sideband

}}}}

}
    }


cout<<"DONE 2"<<endl;

//------------------------------DATA-------------------------------

 //Fill sideband events for Data
  for(int i=nevt;i<nevt2;i++)
//   for(int i=0;i<50000;i++)
   {


      chain->GetEntry(i);
stopbin1=floor((f_Pizmom-1.84)/width1);
stopbin2=floor((f_Pimom-0.86)/width2);




if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 && f_Pizmom < 4.75){
if(f_Pimom > 0.84 && f_Pimom < 4.75){
if(f_PD > 2.62){


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
for(int j=0;j<stopbin1;j++){
for(int k=0;k<stopbin2;k++){siba_data[j][k]++; }}
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
for(int j=0;j<numbins;j++){
eff[i][j] =sig[i][j]/fullarea;
if(siba_MC[i][j] != 0){ratio[i][j] = (siba_data[i][j])/siba_MC[i][j];} else{ratio[i][j] =1;}

sig[i][j]=sig[i][j]*0.477;
bkg[i][j]=bkg[i][j]*ratio[i][j];


fom[i][j] = sig[i][j]/sqrt(sig[i][j]+bkg[i][j]); //cout<<sig[i]<<endl;
 if(fom[i][j] > bestfom){bestfom =fom[i][j]; besti=i; bestj=j;}

}}
// cout<<"-------------------------DONE-----------------------"<<endl;

count=;
//Plotting
//for(int i=0;i<numbins;i++){
//for(int j=0;j<numbins;j++){
for(int i=0;i<200;i++){
for(int j=0;j<200;j++){
 count++;
 x=cut1[i];
 y=cut2[j];
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






cout<<"best f.o.m. ="<<bestfom<<endl;
cout<<"best cut(Mom_Pi0) =  "<<cut1[besti]<<"\t best cut(Mom_Pi+) =  "<<cut2[bestj]<<endl;
//cout<<"best cut = +/- "<<cut1[besti]*0.0055<<" ("<<cut1[besti]<<" sigma cut)"<<endl;




}



