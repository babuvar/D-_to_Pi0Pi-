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


  chain->Add("DtoPiZPi_50000x2_Y4s_conv.root");
//  chain->Add("DtoPiZPi_50000x2_Y4s.root");
//  chain->Add("DtoPiZPi_exp55_Data_conv.root");
//  chain->Add("DtoPiZPi_exp55_Data.root");



  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t f_dzero, Deltam, f_PD, f_Egam1s, f_Egam2s, f_Egam1, f_Egam2, f_Gam1thet, f_Gam2thet, f_MPi0, f_MPi0s, f_type, f_Dstf, f_Df, f_Did;

  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egam1);
  h1->SetBranchAddress("Egamma2",&f_Egam2);
  h1->SetBranchAddress("Gam1thet",&f_Gam1thet);
  h1->SetBranchAddress("Gam2thet",&f_Gam2thet);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizmass",&f_MPi0);
  h1->SetBranchAddress("Pizsmass",&f_MPi0s);
  h1->SetBranchAddress("Type",&f_type);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Did",&f_Did);


int photon1cutflag=0, photon2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

//if(Deltam > 0.139 && Deltam < 0.142){
if(f_dzero > 1.7 && f_dzero < 2.0){
//if(f_type == 2 || f_type == 3|| f_type == 4){
if(f_type == 1){
//if(fabs(f_MPi0-0.135) < 0.016){
//if(fabs(f_MPi0s-0.135) < 0.016){
//if(f_PD > 2.9){

//Hard Photon 1
if(f_Gam1thet < -60 && f_Egam1 > 0.150){photon1cutflag=1;}  //BWD
else if(f_Gam1thet > 73 && f_Egam1 > 0.100){photon1cutflag=1;}  //FWD
else if(f_Gam1thet > -60 && f_Gam1thet < 73 && f_Egam1 > 0.050){photon1cutflag=1;}  //BARREL
else{photon1cutflag=0;}

//Hard Photon 2
if(f_Gam2thet < -60 && f_Egam2 > 0.150){photon2cutflag=1;}  //BWD
else if(f_Gam2thet > 73 && f_Egam2 > 0.100){photon2cutflag=1;}  //FWD
else if(f_Gam2thet > -60 && f_Gam2thet < 73 && f_Egam2 > 0.050){photon2cutflag=1;}  //BARREL
else{photon2cutflag=0;}



//Photon cuts
//if(photon1cutflag == 1 && photon2cutflag == 1){
//if(f_Egam1s > 0.030 && f_Egam2s > 0.030){

//if(f_Df != 1){
//if(fabs(f_Did) == 411){

data->add(RooArgSet(dzero));

//}


//}}//Photon cuts

}}//}//}}
    }
 




  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,1000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,1000000);
//Signal
  RooRealVar m("#mu","m",1.8690531,1.80,1.90);//,1.86,1.80,1.90);//
  RooRealVar s("#sigma","s",0.0160014,0,0.04);//,0.01,0,0.04);//
  RooRealVar a("#alpha","a",0.668437,0,2);//,0.5,0,2);
  RooRealVar n("n","n",14.59,0,150);//2,0,150);


  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n); 
//Background
  RooRealVar p("p","p",0.0);   
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooRealVar p1("p1","p1",0.5,0,1);   
  RooChebychev bkg3("bkg3", "bkg3", dzero, RooArgList(p1));

//  RooRealVar c("c","c",-9,-30,10);
//  RooExponential bkg2("bkg2", "bkg2", dzero, c);
//  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.0);//,0.5,0,1); 
//  RooAddPdf Bkg("Bkg","Bkg",RooArgList(bkg1,bkg2),RooArgList(f));


RooRealVar coef_x1("coef_x1", "coef_x0", 1.0,-1000,1000); 
// RooRealVar coef_x1("coef_x1", "coef_x0", 0.509);
RooRealVar coef_x2("coef_x2", "coef_x2", 1.0,-1000,1000);
RooRealVar coef_x3("coef_x3", "coef_x3", 1.0,-1000,1000);
RooPolynomial Bkg("Bkg","c1 dzero + c2 dzero2 + c3 dzero3 ", dzero ,RooArgList(coef_x1,coef_x2,coef_x3));
//RooChebychev Bkg("Bkg","1+c1 dzero + c2 (dzero2-1)", dzero ,RooArgList(coef_x1, coef_x2));



//Full Model
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));

//Fit
   RooFitResult* fitRes = Model.fitTo(*data,Save());
//   RooFitResult* fitRes = Sig.fitTo(*data);

  TCanvas* can = new TCanvas("c","c",800,800) ;
           can->cd();
  RooPlot* dzero_frame = dzero.frame(Bins(30),Title("Mass of D candidate"));
  dzero_frame->SetTitle("Mass of D candidate");


  data->plotOn(dzero_frame);

  Model.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
  Model.plotOn(dzero_frame,Components("Bkg"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
  Model.paramOn(dzero_frame);
//  Sig.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2)); 
//  Sig.paramOn(dzero_frame);


dzero_frame->Draw();


/*
  RooPlot *yframe1 = data->plotOn(dzero.frame(100),MarkerColor(kBlue));

  TGaxis::SetMaxDigits(3);
  yframe1->Draw();
*/

}

