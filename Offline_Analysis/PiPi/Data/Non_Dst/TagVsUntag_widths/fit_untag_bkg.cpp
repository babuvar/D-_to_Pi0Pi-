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

void fit_untag_bkg(void)
{
  //LOAD DATA FILE
  TChain* chain2=new TChain("h2");


  chain2->Add("pipi_4smc_0F.root");



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
//  for(int i=0;i<1000;i++) 
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
if(f_Did != 411 && f_Did != -411  && f_Did != 421 && f_Did != -421 && f_Did != 431  && f_Did != -431 ){

data2->add(RooArgSet(dzero));}}

}

}}}}}
    }
 


//PDF
  RooRealVar p("p","p",0.0,-1,1);
  RooRealVar p1("p1","p1",0.0,-5,5);
  RooRealVar c("c","c",-20,-50,0.1);
//  RooChebychev Bkg("Bkg", "Bkg", dzero, RooArgList(p,p1));
//  RooExponential Bkg2("Bkg2", "Bkg2", dzero, c);
  RooExponential Bkg("Bkg", "Bkg", dzero, c);

  RooRealVar f("f","f",0.5,0.0,1.0);
//  RooAddPdf Bkg("Sig","Sig",RooArgList(Bkg1,Bkg2),RooArgList(f));



//Fit
   RooFitResult* fitRes2 = Bkg.fitTo(*data2,Save());

  RooPlot *xframe_2 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  data2->plotOn(xframe_2);
  Bkg.plotOn(xframe_2, LineStyle(kSolid),LineWidth(2));
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(15);
    frame2->SetMinimum(-15);
  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;
  Bkg.paramOn(xframe_2);



  TCanvas* can = new TCanvas("c","c",700,700) ;
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();
 
  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;


}










