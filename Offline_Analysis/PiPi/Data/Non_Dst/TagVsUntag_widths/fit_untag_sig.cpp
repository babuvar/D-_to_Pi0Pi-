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

double low_lim = 1.68, hi_lim = 2.06;
//double low_lim = 1.7, hi_lim = 1.95;


  RooRealVar dzero("dzero","dzero",low_lim,hi_lim,"GeV");

  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero));

void fit_untag_sig(void)
{
  //LOAD DATA FILE
  TChain* chain2=new TChain("h2");


  chain2->Add("pipi_4smc_0F.root");
  chain2->Add("pipi_5s_gmc0.root");

  Int_t nevt2=(int)chain2->GetEntries();
  cout<<"nevt2\t"<<nevt2 <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge;





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


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;




  for(int i=0;i<nevt2;i++)
//   for(int i=0;i<10000;i++)
    {


      chain2->GetEntry(i);
		  dzero.setVal(f_dzero);
 
if(f_dzero >low_lim && f_dzero < hi_lim){          
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
//if(f_PD > 2.5){
if(f_PD > 2.65){
//if(f_PD > 3.0){


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
if(f_Df==1 || f_Df==10){data2->add(RooArgSet(dzero));}
}

}}}}}
    }
 



//Common Pars - II
  RooRealVar m_diff("#mu_{diff}","m_diff",-0.000635,-0.1,0.1);
  RooRealVar s_ratio("#sigma_{ratio}","s_ratio",1.000,0.5,10.0);

  RooRealVar m("#mu","m", 1.86861,1.86,1.88);
  RooRealVar s("#sigma","s", 0.01532,0.005,0.03);
  RooRealVar a("#alpha","a",0.7937,0.65,2.15);
  RooRealVar n("n","n",16.44,0,1000);
  RooFormulaVar m_g("#mu_{G}","#mu+#mu_{diff}",RooArgList(m,m_diff));
  RooFormulaVar s_g("#sigma_{G}","#sigma*#sigma_{ratio}",RooArgList(s,s_ratio));
//  RooRealVar m_g("#mu_{G}","m_g", 1.86861,1.780,1.98);
//  RooRealVar s_g("#sigma_{G}","s_g", 0.01532,0.0,0.15);


//Signal  - II
  RooCBShape Sig2("Sig2", "Cystal Ball Function", dzero, m, s, a, n); 
//Background 
  RooGaussian Bkg2("Bkg2", "Bkg2",dzero,m_g,s_g);

  RooRealVar f("frac","f",0.99,0.0,1.0);
//Full Model  - II
  RooAddPdf Model2("Model2","Model2",RooArgList(Sig2,Bkg2),RooArgList(f));


//Fit
//   RooFitResult* fitRes2 = Sig2.fitTo(*data2,Save());
   RooFitResult* fitRes2 = Model2.fitTo(*data2,Save());


 
  RooPlot *xframe_2 =dzero.frame(Bins(25));
  data2->plotOn(xframe_2);
  Model2.plotOn(xframe_2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
  Model2.paramOn(xframe_2);
//    Sig2.plotOn(xframe_2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
//    Sig2.paramOn(xframe_2);

    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame() ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
//    frame2->SetMaximum(30);
//    frame2->SetMinimum(-30);

  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;

  Model2.plotOn(xframe_2,Components("Bkg2"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));

  TCanvas* can = new TCanvas("c","c",700,700) ;
  can->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;

/*
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();

  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;
*/



}










