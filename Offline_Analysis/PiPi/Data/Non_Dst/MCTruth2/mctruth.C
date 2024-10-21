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
void mctruth(void)
{

float Nsignal_p=0, Nsignal_n=0, A_raw, e_Araw;

float Nbkg_p=0, Nbkg_n=0, A_bkg, e_Abkg;

float Nb1_p=0, Nb1_n=0, A_b1, e_Ab1;

float Nb2_p=0, Nb2_n=0, A_b2, e_Ab2;





    TH1F *h5 = new TH1F("my_hist5","hist5",50,1.68,2.06);  //d0 mass
    TH1F *h6 = new TH1F("my_hist6","hist6",50,1.68,2.06);  //d0 mass
    TH1F *h7 = new TH1F("my_hist7","hist7",50,1.68,2.06);  //d0 mass


    TH1F *h8 = new TH1F("my_hist8","hist8",50,1.68,2.06);  //d0 mass
    TH1F *h9 = new TH1F("my_hist9","hist9",50,1.68,2.06);  //d0 mass
    TH1F *h10 = new TH1F("my_hist10","hist10",50,1.68,2.06);  //d0 mass




  h5->SetFillColor(4);
  h6->SetFillColor(8);
  h7->SetFillColor(2);


  h8->SetFillColor(3);
  h9->SetFillColor(7);
  h10->SetFillColor(5);



  //LOAD DATA FILE
  TChain* chain=new TChain("h2");
  chain->Add("pipi_mc.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df, f_Did; 
    Float_t   f_Dcharge;


  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Df",&f_Df);
  h2->SetBranchAddress("Did",&f_Did);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);


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



if(f_Df == 1 || f_Df == 10){  h5->Fill(f_dzero); }
//else if(f_Did == 411 || f_Did == -411  || f_Did == 421 || f_Did == -421 || f_Did == 431  || f_Did == -431 ){  h6->Fill(f_dzero); }
else if(f_Did == 411 || f_Did == -411){  h8->Fill(f_dzero); }
else if(f_Did == 421 || f_Did == -421){  h9->Fill(f_dzero); }
else if(f_Did == 431 || f_Did == -431){  h10->Fill(f_dzero); }

else{  h7->Fill(f_dzero); }



//if(f_Df != 1 && f_Df != 10 && f_Df != -1 && f_Df != 0){  cout<<f_Did<<endl; }

}}}}}}
    }
 




  hs2->Add(h5);
  hs2->Add(h7);
  hs2->Add(h6);
  hs2->Add(h5);


  hs3->Add(h10);
  hs3->Add(h8);
  hs3->Add(h9);





TLegend * leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->SetHeader("Components");
  leg->AddEntry(h6,"Background from D^{+}, D^{0} and D_s","f");
  leg->AddEntry(h7,"Combinatorial","f");
  leg->AddEntry(h5,"Signal","f");

TLegend * leg2 = new TLegend(0.6,0.7,0.89,0.89);
  leg2->SetHeader("Components");
  leg2->AddEntry(h8,"Background from D^{+} ","f");
  leg2->AddEntry(h9,"Background from  D^{0} ","f");
  leg2->AddEntry(h10,"Background from  D_s ","f");




    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
  c1->Divide(2,2) ;

c1->cd(1); hs2->Draw();  leg->Draw();
c1->cd(2); h5->Draw();
c1->cd(3); h6->Draw();
c1->cd(4); h7->Draw();

    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
  c2->Divide(2,2) ;
c2->cd(1); hs3->Draw();  leg2->Draw();
c2->cd(2); h8->Draw();
c2->cd(3); h9->Draw();
c2->cd(4); h10->Draw();

}





















