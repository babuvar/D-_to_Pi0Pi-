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
#include "TPaveLabel.h"
using namespace RooFit;
  RooRealVar dzero("dzero","M_{D}   ",1.68,2.06,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));



int nsig1=0,nsig2=0,psig1=0,psig2=0;

void pkbkg3()
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");


//  chain->Add("../Root_Files/pipi_4smc_0F.root");
//  chain->Add("../Root_Files/pipi_4smc_1F.root"); 
//  chain->Add("../Root_Files/pipi_4smc_2F.root");
//  chain->Add("../Root_Files/pipi_4smc_3F.root"); 
//  chain->Add("../Root_Files/pipi_4smc_4F.root"); 
  chain->Add("../Root_Files/pipi_4smc_5F.root");


//  chain->Add("../Root_Files/pipi_5s_gmc0.root"); 
//  chain->Add("../Root_Files/pipi_5s_gmc1.root"); 
//  chain->Add("../Root_Files/pipi_5s_gmc2.root");
//  chain->Add("../Root_Files/pipi_5s_gmc3.root"); 
//  chain->Add("../Root_Files/pipi_5s_gmc4.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc5.root");




  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df,f_Did;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
//  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Did",&f_Did);
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



int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<5000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){


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


if( f_Df != 1 && f_Df != 10){
//Optimized P*D* bin
if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD < 2.95){

//if(f_Dcharge == 1){
//if(f_Dcharge == -1){
data->add(RooArgSet(dzero));
//}

}
}


}}//Photon cuts

}}}}}}//}//}
    }
 


//Common Pars - II
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.0,-0.04,0.04);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0,0.5,2.5);


  RooFormulaVar m("m_sig","1.6802+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.03656*#sigma_{fudge}",RooArgList(s_fudge));
//  RooRealVar m("m","m",1.6802);
//  RooRealVar s("s","s",0.03656);
  RooRealVar a("#alpha","a",-1.8526);
  RooRealVar n("n","n",1.575);
//  RooRealVar p("p","p",-0.4235);

/*
  RooRealVar m("m","m",1.67885);
  RooRealVar s("s","s",0.03652);
  RooRealVar a("#alpha","a2",-1.8661);
  RooRealVar n("n","n",1.291);
//  RooRealVar p("p","p",-0.28888);
*/

  RooRealVar p("p","p",-0.4235,-1.0,1.0);


  RooRealVar f("f","f",0.5,0.0,1.0);
  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n);
   RooChebychev Bkg("Bkg", "Bkg", dzero, RooArgList(p));
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(f));


//Fit
   RooFitResult* fitRes2 = Model.fitTo(*data,Save());


  TCanvas* can = new TCanvas("c","c",700,700) ;
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();
 
  RooPlot *xframe_2 =dzero.frame(Bins(38));
  data->plotOn(xframe_2);
  Model.plotOn(xframe_2, LineStyle(kSolid),LineWidth(2));
  Model.paramOn(xframe_2);
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame() ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;
  Model.plotOn(xframe_2,Components("Sig"),LineStyle(kDashed),LineColor(kRed));
  Model.plotOn(xframe_2,Components("Bkg"),LineStyle(kDashed));



  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;

}




   






