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
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
void sig(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");


  chain->Add("pipi_MC6.root");
//  chain->Add("pipi_MC7.root");
//  chain->Add("pipi_1M.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
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

//if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
//if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
//if(f_Pizmom > 1.06 ){
//if(f_Pimom > 0.84 ){
//if(f_PD > 2.95){
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
//if(photon1cutflag == 1 && photon2cutflag == 1){
//if(sphoton1cutflag == 1 && sphoton2cutflag == 1){

//if(f_Df!=1 && f_Df!=10){
if(f_Df==1 || f_Df==10){
data->add(RooArgSet(dzero));
cout<<"Df = "<<f_Df<<endl;
}



num_pass++;
//}}//Photon cuts

}//}}}}}}
    }
 
cout<<"num_pass"<<num_pass<<endl;




//Common Pars - II
  RooRealVar m2("m2","m2", 1.86861,1.86,1.88);
  RooRealVar s2("s2","s2", 0.01532,0.005,0.03);
  RooRealVar a2("#alpha2","a2",0.7937,0.65,2.15);
  RooRealVar n2("n2","n2",16.44,0,1000);
  RooRealVar m2_g("m2_g","m2_g", 1.86861,1.780,1.98);
  RooRealVar s2_g1("s2_g1","s2_g1", 0.01532,0.0,0.15);
  RooRealVar s2_g2("s2_g2","s2_g2", 0.01532,0.0,0.15);

//Signal  - II
  RooCBShape Sig2("Sig2", "Cystal Ball Function", dzero, m2, s2, a2, n2); 
//Background 
//  RooGaussian Bkg2("Bkg2", "Bkg2",dzero,m2_g,s2_g1);
//  RooGaussian Bkg2("Bkg2", "Bkg2",dzero,m2,s2_g1);

//  RooBifurGauss Bkg2("Bkg2", "Bkg2",dzero,m2_g,s2_g1,s2_g2);
  RooBifurGauss Bkg2("Bkg2", "Bkg2",dzero,m2,s2_g1,s2_g2);

  RooRealVar f("f","f",0.99,0.8,1.0);


//Full Model  - II
//  RooAddPdf Model2("Model2","Model2",RooArgList(Sig2,Bkg2),RooArgList(N_t2,N_tb2));
  RooAddPdf Model2("Model2","Model2",RooArgList(Sig2,Bkg2),RooArgList(f));


//Fit
   RooFitResult* fitRes2 = Sig2.fitTo(*data,Save());
//   RooFitResult* fitRes2 = Model2.fitTo(*data,Save());




 
  RooPlot *xframe_2 =dzero.frame(Bins(25));
  data->plotOn(xframe_2);
//  Model2.plotOn(xframe_2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
//  Model2.paramOn(xframe_2);
    Sig2.plotOn(xframe_2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
    Sig2.paramOn(xframe_2);

    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame() ;
    frame2->addPlotable(hpull2,"P") ;
   frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
//    frame2->SetMaximum(10);
//    frame2->SetMinimum(-10);

  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;

//  Model2.plotOn(xframe_2,Components("Bkg2"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));

  TCanvas* can = new TCanvas("c","c",700,700) ;
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();

  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;


}














