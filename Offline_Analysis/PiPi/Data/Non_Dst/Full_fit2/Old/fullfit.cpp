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
//  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");


TH1D* data=new TH1D("pdat", "pdat", 100, 1.68, 2.06);


void fullfit(void)
{
  //LOAD DATA FILE
//  TChain* chain=new TChain("h2");
  TChain* chain=new TChain("h2");


//  chain->Add("pipi_d4s.root");
  chain->Add("pipi_mc.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge;


//  h2->SetBranchAddress("Deltam",&Deltam);
  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
//  h2->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
//  h2->SetBranchAddress("Dstf",&f_Dstf);
  h2->SetBranchAddress("Df",&f_Df);
//  h2->SetBranchAddress("Egam1s",&f_Egam1s);
//  h2->SetBranchAddress("Egam2s",&f_Egam2s);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
//  h2->SetBranchAddress("Categ",&f_Categ);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


  for(int i=0;i<nevt;i++) 
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);


if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 3.5){


//if(f_Egamma1 > 0.150 && f_Egamma2 > 0.150){

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




//if(f_Df==1 || f_Df==10){
//if(f_Df != 1 && f_Df != 10){
data->Fill(f_dzero);
//}


}
}}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;


RooDataHist("Data","Data",dzero,data) ; 



  //DEFINE PDF
//Common
  //____________________________________________________

  RooRealVar N_t("N_{Sig}","N_t",100,0,1000000);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,200000000);


  //_____________________________________________________

//Common Pars
  RooRealVar m("m","m", 1.86805);//,1.86,1.88);
  RooRealVar s("s","s", 0.015559);//,0.005,0.03);
  RooRealVar a("#alpha","a",0.901);//,0.65,1.15);
  RooRealVar n("n","n",2.337);//,0,1000);
  RooRealVar m_g("m2_g","m2_g", 1.8942);//,1.780,1.98);
  RooRealVar s_g("s2_g","s2_g", 0.0441);//,0.0,0.15);



  RooRealVar p("p","p",0.0);
  RooRealVar c("c","c",-20,-50,0.1);
  RooRealVar mean_bkg("#mu_{bkg}","mean_bkg",1.6799,1.5,1.70);
  RooRealVar sig_bkg("#sigma_{bkg}","sig_bkg",0.038,0.0,0.1);
 
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.9,0.0,1);
 RooRealVar f2("frac_{Bkg/Bkg3}","f3",0.9,0.0,1);



//Signal 

  RooCBShape Sig1("Sig1", "Cystal Ball Function", dzero, m, s, a, n); 
  RooGaussian Sig2("Sig2", "Sig2",dzero,m_g,s_g);
  RooRealVar f_sig("f_sig","f_sig",0.9566);//,0.8,1.0);
  RooAddPdf Sig("Sig","Sig",RooArgList(Sig1,Sig2),RooArgList(f_sig));



//Background 
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooExponential bkg2("bkg2", "bkg2", dzero, c);
  RooGaussian bkg3("bkg3", " Gaussian bkg",dzero,mean_bkg,sig_bkg);
  RooAddPdf Bkg1("Bkg1","Bkg1",RooArgList(bkg1,bkg2),RooArgList(f));
  RooAddPdf Bkg2("Bkg2","Bkg2",RooArgList(Bkg1,bkg3),RooArgList(f2));

//Full Model 
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg2),RooArgList(N_t,N_tb));





//Fit
   RooFitResult* fitRes = Model.fitTo(*Data,Save());



  TCanvas* can = new TCanvas("c","c",700,700) ;


  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  Data->plotOn(xframe_1);
  Model.plotOn(xframe_1);
  Model.plotOn(xframe_1,Components("Sig"),LineStyle(kDashed),LineColor(kRed));
  Model.plotOn(xframe_1,Components("Bkg2"),LineStyle(kDashed));
//  Model.plotOn(xframe_1,Components("bkg1"),LineStyle(kDashed),LineColor(kGreen));
//  Model.plotOn(xframe_1,Components("bkg2"),LineStyle(kDashed),LineColor(kGreen));
//  Model.plotOn(xframe_1,Components("bkg3"),LineStyle(kDashed),LineColor(kGreen));
   Model.paramOn(xframe_1);



  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;

cout<<"num_pass"<<num_pass<<endl;

}










