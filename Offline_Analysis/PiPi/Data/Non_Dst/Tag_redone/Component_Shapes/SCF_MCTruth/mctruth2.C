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




    TH1F *h8 = new TH1F("my_hist8","hist8",10,1.68,2.06);  //d0 mass

    TH1F *h5 = new TH1F("my_hist5","hist5",38,1.68,2.06);  //d0 mass
    TH1F *h6 = new TH1F("my_hist6","hist6",38,1.68,2.06);  //d0 mass



  h8->SetFillColor(3);

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  chain->Add("mc.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df, f_Did; 
    Float_t   f_Dcharge, f_Genpizid, f_Genpiid , f_Piz1f, f_Pif ;


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
  h1->SetBranchAddress("Genpizid",&f_Genpizid);
  h1->SetBranchAddress("Genpiid",&f_Genpiid);
  h1->SetBranchAddress("Piz1f",&f_Piz1f);
  h1->SetBranchAddress("Pif",&f_Pif);


  THStack *hs2 = new THStack("hs2","D Mass");

   THStack *hs3 = new THStack("hs3","D Mass");




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {


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


if(f_Df == -1){h5->Fill(f_dzero);}
else{if(f_Df != 1 && f_Df != 10 ){h6->Fill(f_dzero);
cout<<"Df = "<<f_Df<<"\t\t Did = "<<f_Did<<endl;
cout<<"-----------------------------------------------------"<<endl;
}}


if((f_Df != 1 && f_Df != 10) && (f_Did == 411 || f_Did == -411)){
if(f_Genpizid ==111 && (f_Genpiid == 211 || f_Genpiid == -211)){

h8->Fill(f_dzero); 
/*
cout<<"**********************************"<<endl;
cout<<"Df = "<<f_Df<<endl;
cout<<"Piz1f = "<<f_Piz1f<<endl;
cout<<"Pif = "<<f_Pif<<endl;
cout<<"**********************************"<<endl;
*/
}}




//if(f_Df != 1 && f_Df != 10 && f_Df != -1 && f_Df != 0){  cout<<f_Did<<endl; }

}}}}}}
    }
 







    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
c2->cd(); h8->Draw();

    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");
c3->cd(); h5->Draw();


    TCanvas *c4 = new TCanvas("myCanvas4","My Canvas4");
c4->cd(); h6->Draw();

}





















