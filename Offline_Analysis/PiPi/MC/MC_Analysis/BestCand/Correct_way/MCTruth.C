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
void MCTruth(void)
{

float Nsignal_p=0, Nsignal_n=0, A_raw, e_Araw;

float Nbkg_p=0, Nbkg_n=0, A_bkg, e_Abkg;

float Nb1_p=0, Nb1_n=0, A_b1, e_Ab1;

float Nb2_p=0, Nb2_n=0, A_b2, e_Ab2;


    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");


    TH1F *h11 = new TH1F("my_hist11","hist11",50,0.135,0.155);  //deltam
    TH1F *h2 = new TH1F("my_hist2","hist2",50,0.135,0.155);  //deltam
    TH1F *h3 = new TH1F("my_hist3","hist3",50,0.135,0.155);  //deltam
    TH1F *ha = new TH1F("my_hista","hista",50,0.135,0.155);  //deltam
    TH1F *hd = new TH1F("my_histd","histd",50,0.135,0.155);  //deltam
//    TH1F *h4 = new TH1F("my_hist4","hist4",50,0.135,0.155);  //deltam

    TH1F *h5 = new TH1F("my_hist5","hist5",50,1.70,2.0);  //d0 mass
    TH1F *h6 = new TH1F("my_hist6","hist6",50,1.70,2.0);  //d0 mass
    TH1F *h7 = new TH1F("my_hist7","hist7",50,1.70,2.0);  //d0 mass
    TH1F *hb = new TH1F("my_histb","histb",50,1.70,2.0);  //d0 mass
    TH1F *he = new TH1F("my_histe","histe",50,1.70,2.0);  //d0 mass
//    TH1F *h8 = new TH1F("my_hist8","hist8",50,1.70,2.0);  //d0 mass


    TH1F *h8 = new TH1F("my_hist8","hist5",50,2.5,5.0);  //P*(D*)
    TH1F *h9 = new TH1F("my_hist9","hist6",50,2.5,5.0);  //P*(D*)
    TH1F *h10 = new TH1F("my_hist10","hist7",50,2.5,5.0);  //P*(D*)
    TH1F *hc = new TH1F("my_histc","histc",50,2.5,5.0);  //P*(D*)
    TH1F *hf = new TH1F("my_histf","histf",50,2.5,5.0);  //P*(D*)


  h11->SetFillColor(4);
  h2->SetFillColor(8);
  h3->SetFillColor(2);
//  h4->SetFillColor(6);
  ha->SetFillColor(6);
  hd->SetFillColor(5);

  h5->SetFillColor(4);
  h6->SetFillColor(8);
  h7->SetFillColor(2);
//  h8->SetFillColor(6);
  hb->SetFillColor(6);
  he->SetFillColor(5);


  h8->SetFillColor(4);
  h9->SetFillColor(8);
  h10->SetFillColor(2);
  hc->SetFillColor(6);
  hf->SetFillColor(5);




  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  chain->Add("pipi_100K_after.root");
//  chain->Add("pipi_100K_before.root");
//  chain->Add("DtoPiZPi_50000x2_Y4s_old.root");
//  chain->Add("pipi_MC1.root");
  chain->Add("pipi_MC0_Y4s_BCS.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df, f_Did;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Did",&f_Did);

  THStack *hs1 = new THStack("hs1","M_{D*} - M_{D}");
  THStack *hs2 = new THStack("hs2","D Mass");
  THStack *hs3 = new THStack("hs3","P*_{D*}");
 







  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {


//      chain->GetEntry(i);
      chain->GetEntry(i);


if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
//if(f_dzero >1.70  && f_dzero < 2.0){       
if(f_dzero >1.68  && f_dzero < 2.0){

//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD < 2.95 ){

//if(f_Categ == 5){
//if(f_Categ == 1 || f_Categ == 2 || f_Categ == 6 ){
//if(f_Categ == 1 || f_Categ == 2 || f_Categ == 5 || f_Categ == 6 ){

//if(f_Categ != 1 && f_Categ != 2 && f_Categ != 6 && f_Categ != 5){
/*

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
if(photon1cutflag == 1 && photon2cutflag == 1){
if(sphoton1cutflag == 1 && sphoton2cutflag == 1){
*/


if(f_Df == 1){	if(f_Dcharge == 1){Nsignal_p++;}	if(f_Dcharge == -1){Nsignal_n++;}	}
else{	if(f_Dcharge == 1){Nbkg_p++;}	if(f_Dcharge == -1){Nbkg_n++;}	

if(f_Did == 411 || f_Did == -411  ){	if(f_Dcharge == 1){Nb1_p++;}	if(f_Dcharge == -1){Nb1_n++;}		}
if(f_Did == 421 || f_Did == -421  ){	if(f_Dcharge == 1){Nb2_p++;}	if(f_Dcharge == -1){Nb2_n++;}		}
}



if(f_Dstf == 1){h11->Fill(Deltam); h5->Fill(f_dzero); h8->Fill(f_PD);}
//else{
//if(f_Df == 1){h2->Fill(Deltam); h6->Fill(f_dzero); h9->Fill(f_PD);}
else if(f_Df == 1){h2->Fill(Deltam); h6->Fill(f_dzero); h9->Fill(f_PD);}
else if(f_Did == 411 || f_Did == -411  ){ ha->Fill(Deltam); hb->Fill(f_dzero); hc->Fill(f_PD); }
else if(f_Did == 421 || f_Did == -421  ){ hd->Fill(Deltam); he->Fill(f_dzero); hf->Fill(f_PD); }
else{ h3->Fill(Deltam); h7->Fill(f_dzero); h10->Fill(f_PD); }
//else{if(f_Did == 411 || f_Did == -411  ){ h3->Fill(Deltam); h7->Fill(f_dzero); h10->Fill(f_PD); }}
//}


//}}//Photon cuts

}}}}}}}//}
    }
 



  hs1->Add(h2);
  hs1->Add(h3);
  hs1->Add(h11);
  hs1->Add(ha);
  hs1->Add(hd);


  hs2->Add(h7);
  hs2->Add(hb);
  hs2->Add(he);
  hs2->Add(h6);
  hs2->Add(h5);


  hs3->Add(h10);
  hs3->Add(h9);
  hs3->Add(h8);
  hs3->Add(hc);
  hs3->Add(hf);


TLegend * leg = new TLegend(0.6,0.7,0.89,0.89);
//  leg = new TLegend();
  leg->SetHeader("Reconstructed Object");
  leg->AddEntry(h11,"Signal","f");
  leg->AddEntry(h2,"Wrong slow #pi^{0}","f");
  leg->AddEntry(h3,"combinatorial","f");
  leg->AddEntry(ha,"D+ contamination","f");
  leg->AddEntry(hd,"D0 contamination","f");


c1->cd(); hs1->Draw();  leg->Draw();
c2->cd(); hs2->Draw();  leg->Draw();
c3->cd(); hs3->Draw();  leg->Draw();


//Asymmetry calculationsa
A_raw 	= (Nsignal_p-Nsignal_n)/(Nsignal_p+Nsignal_n);
//e_Araw 	= 
A_bkg	= (Nbkg_p-Nbkg_n)/(Nbkg_p+Nbkg_n);
A_b1	= (Nb1_p-Nb1_n)/(Nb1_p+Nb1_n);
A_b2	= (Nb2_p-Nb2_n)/(Nb2_p+Nb2_n);




cout<<"Nsig_p = "<<Nsignal_p<<endl;
cout<<"Nsig_n = "<<Nsignal_n<<endl;

cout<<"Nbkg_p = "<<Nbkg_p<<endl;
cout<<"Nbkg_n = "<<Nbkg_n<<endl;

cout<<"A_raw = "<<A_raw<<endl;
cout<<"A_bkg = "<<A_bkg<<endl;
cout<<"A_b1 = "<<A_b1<<endl;
cout<<"A_b2 = "<<A_b2<<endl;




}















