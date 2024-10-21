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
  RooRealVar dzero("dzero","dzero",1.68,2.0,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
void fit_MD(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");


  chain->Add("pipi_data_BkgShp.root");
//  chain->Add("Merged_GMC.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df;


  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Deltam",&Deltam);

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
if(Deltam > 0.145 && Deltam < 0.148){       //optimized for M_D fitting
if(f_dzero > 1.68 && f_dzero < 2.0 ){
data->add(RooArgSet(dzero));
}}

}


cout<<"num_pass"<<num_pass<<endl;



  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,10000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,10000000);

  RooRealVar mean("mean","mean",1.65,1.5,1.9);
  RooRealVar sig("sig","sig",0.009,0.00,0.25);
  RooRealVar sig1("sig1","sig1",0.005,0.00,0.2);
  RooRealVar sig2("sig2","sig2",0.005,0.00,0.2);
  RooRealVar f_Sig("f_Sig","f_Sig",0.3,0.0,1.0);


//Signal
  RooGaussian sig_p1("sig_p1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_p2("sig_p2", "signal Gaussian 2",dzero,mean,sig1,sig2);
  RooAddPdf Sig("Sig_p","Sig",RooArgList(sig_p1,sig_p2),f_Sig);

  RooRealVar c("c","c",-9,-30,10);
  RooExponential SigP("SigP", "SigP", dzero, c);

//Background
  RooRealVar p("p","p",0.0,-5,0.5);   
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooRealVar p1("p1","p1",0.5,0,1);   
  RooChebychev bkg3("bkg3", "bkg3", dzero, RooArgList(p1));


  RooExponential bkg2("bkg2", "bkg2", dzero, c);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
  RooAddPdf Bkg("Bkg","Bkg",RooArgList(bkg1,bkg2),RooArgList(f));


//Full Model
//  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));
//  RooAddPdf Model("Model","Model",RooArgList(Sig,bkg1),RooArgList(N_sig,N_bkg));
  RooAddPdf Model("Model","Model",RooArgList(sig_p1,bkg1),RooArgList(N_sig,N_bkg));
//  RooAddPdf Model("Model","Model",RooArgList(SigP,bkg1),RooArgList(N_sig,N_bkg));

//Fit
   RooFitResult* fitRes = Model.fitTo(*data,Save());


  TCanvas* can = new TCanvas("c","c",700,800) ;
           can->cd();
  RooPlot* dzero_frame = dzero.frame(Bins(30),Title("Mass of D candidate"));
  dzero_frame->SetTitle("Mass of D candidate");


  data->plotOn(dzero_frame);

  Model.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
//  Model.plotOn(dzero_frame,Components("Bkg"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
  Model.plotOn(dzero_frame,Components("bkg1"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
 Model.paramOn(dzero_frame);



dzero_frame->Draw();




}

