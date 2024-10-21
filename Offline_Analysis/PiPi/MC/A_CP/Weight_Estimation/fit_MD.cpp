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
  RooRealVar dzero("dzero","dzero",1.70,2.00,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
void fit_MD(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");


  chain->Add("pipi_MC6.root");
//  chain->Add("pipi_data.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
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

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.70  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD < 5.0){

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

data->add(RooArgSet(dzero));

num_pass++;
}}//Photon cuts

}}}}}}}
    }
 
cout<<"num_pass"<<num_pass<<endl;



  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,1000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,1000000);
//Signal
  RooRealVar m("#mu","m",1.8690531,1.85,1.88);//,1.86,1.80,1.90);//
  RooRealVar s("#sigma","s",0.0160014,0,0.04);//,0.01,0,0.04);//
  RooRealVar a("#alpha","a",0.668437,0,2);
  RooRealVar n("n","n",14.59,0,200);


  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n); 
//Background
  RooRealVar p("p","p",0.0);   
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooRealVar p1("p1","p1",0.5,0,1);   
  RooChebychev bkg3("bkg3", "bkg3", dzero, RooArgList(p1));

  RooRealVar c("c","c",-9,-30,10);
  RooExponential bkg2("bkg2", "bkg2", dzero, c);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
  RooAddPdf Bkg("Bkg","Bkg",RooArgList(bkg1,bkg2),RooArgList(f));


//Full Model
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));

//Fit
   RooFitResult* fitRes = Model.fitTo(*data,Save());


  TCanvas* can = new TCanvas("c","c",700,800) ;
           can->cd();
  RooPlot* dzero_frame = dzero.frame(Bins(30),Title("Mass of D candidate"));
  dzero_frame->SetTitle("Mass of D candidate");


  data->plotOn(dzero_frame);

  Model.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
  Model.plotOn(dzero_frame,Components("Bkg"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
  Model.paramOn(dzero_frame);



dzero_frame->Draw();




}

