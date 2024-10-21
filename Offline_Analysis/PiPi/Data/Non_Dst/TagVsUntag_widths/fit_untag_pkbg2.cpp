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
//  RooRealVar dzero("dzero","dzero",1.68,1.95,"GeV");
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");

  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero));
//TH1D* data2=new TH1D("pdat", "pdat", 100, 1.68, 2.06);


void fit_untag_pkbg2(void)
{



  //LOAD DATA FILE
  TChain* chain2=new TChain("h2");

  chain2->Add("pipi_4smc_0F.root");  
  chain2->Add("pipi_4smc_2F.root");  
  chain2->Add("pipi_4smc_4F.root"); 
  chain2->Add("pipi_4smc_1F.root"); 
  chain2->Add("pipi_4smc_3F.root");  
  chain2->Add("pipi_4smc_5F.root"); 


  Int_t nevt2=(int)chain2->GetEntries();
  cout<<"nevt2\t"<<nevt2 <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge, f_Did;





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
  h2->SetBranchAddress("Dcharge",&f_Dcharge);
  h2->SetBranchAddress("Did",&f_Did);

int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;




  for(int i=0;i<nevt2;i++) 
//  for(int i=0;i<500000;i++) 
    {


      chain2->GetEntry(i);
		  dzero.setVal(f_dzero);
 
//if(f_dzero >1.68 && f_dzero < 2.06){      
if(f_dzero >1.68 && f_dzero < 2.06){ 
//if(f_dzero >1.68 && f_dzero < 1.95){     
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){


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

if(f_Df != 1 && f_Df != 10){
//if(f_Did == 411 || f_Did == -411  || f_Did == 421 || f_Did == -421 || f_Did == 431  || f_Did == -431 ){

data2->add(RooArgSet(dzero));}//}
//data2->Fill(f_dzero);}//}
}

}}}}}
    }
 
//RooDataHist("Data","Data",dzero,data2) ; 


//Common Pars - II
  RooRealVar m2("m2","m2", 1.5,1.1,1.7);
  RooRealVar s2("s2","s2", 0.01532,0.001,0.3);
  RooRealVar a2("#alpha2","a2",-0.7937,-3,3);
  RooRealVar n2("n2","n2",16.44,-1000,1000);
  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m2, s2, a2, n2); 

//Fit
   RooFitResult* fitRes2 = Sig.fitTo(*data2,Save());
//   RooFitResult* fitRes2 = Sig.fitTo(*Data,Save());


  TCanvas* can = new TCanvas("c","c",700,700) ;

 
  RooPlot *xframe_2 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
//  Data->plotOn(xframe_2);
  data2->plotOn(xframe_2);
  Sig.plotOn(xframe_2, LineStyle(kSolid),LineWidth(2));
//  Bkg.plotOn(xframe_2,Components("Bkg1"),LineStyle(kDashed));
  Sig.paramOn(xframe_2);

  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;




}









