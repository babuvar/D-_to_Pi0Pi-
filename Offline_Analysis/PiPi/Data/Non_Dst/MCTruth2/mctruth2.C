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
void mctruth2(void)
{

float Nsignal_p=0, Nsignal_n=0, A_raw, e_Araw;

float Nbkg_p=0, Nbkg_n=0, A_bkg, e_Abkg;

float Nb1_p=0, Nb1_n=0, A_b1, e_Ab1;

float Nb2_p=0, Nb2_n=0, A_b2, e_Ab2;



    TH1F *h15 = new TH1F("my_hist15","hist15",50,1.68,2.06);  //d0 mass
    h15->SetFillColor(5);



  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  chain->Add("pipi_mc.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df, f_Did; 
    Float_t   f_Dcharge, f_Truevt;

/*
  h1->SetBranchAddress("Pstd",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Did",&f_Did);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);
  h1->SetBranchAddress("Truevt",f_Truevt);
*/
  h1->SetBranchAddress("Df",&f_Df);


  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {
/*

//      chain->GetEntry(i);
      chain->GetEntry(i);
//if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.52){
//if(f_PD > 2.5 && f_PD < 5.0){
//Hard Photon 1
if(f_Egamma1 > 0.150 && f_Egamma2 > 0.150){
*/

//cout<<"f_Truevt = "<<f_Truevt<<endl;
//cout<<"f_Did = "<<f_Did<<endl;
cout<<"f_Df = "<<f_Df<<endl;

// if(f_Df == 1 ){  h15->Fill(f_dzero); }

//if( f_Truevt == 1){  h15->Fill(f_dzero); }

/*
 if(f_Df == 1 || f_Df == 10){  h5->Fill(f_dzero); }
//else if(f_Did == 411 || f_Did == -411  || f_Did == 421 || f_Did == -421 || f_Did == 431  || f_Did == -431 ){  h6->Fill(f_dzero); }
else if(f_Did == 411 || f_Did == -411){  h8->Fill(f_dzero); }
else if(f_Did == 421 || f_Did == -421){  h9->Fill(f_dzero); }
else if(f_Did == 431 || f_Did == -431){  h10->Fill(f_dzero); }
*/


//}}}}}}//}
    }
 



 




 


    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
 

c1->cd(); h15->Draw();  







}





















