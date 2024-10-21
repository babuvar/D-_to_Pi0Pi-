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


//  chain->Add("pipi_data_4s.root");
  chain->Add("pipi_MC_4s.root");


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
//  for(int i=0;i<2000000;i++) 
//  for(int i=0;i<100000;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

//if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
//if(f_dzero >1.7 && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
//if(f_PD > 2.62){
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



if(f_Dcharge == 1){
if(f_Df==1 || f_Df==10){num_pass++;}
data->Fill(f_dzero);}


}
}}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;


RooDataHist("Data","Data",dzero,data) ; 



  //DEFINE PDF
//Common
  //____________________________________________________

  RooRealVar N_t("N_{Sig}","N_t",100,0,1000000);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,2000000);


  //_____________________________________________________

//Common Pars
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.0004119,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.072552,0.0,1.9);

  RooFormulaVar m("m","1.86861+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01532*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.7937);//,0.75,0.85);
  RooRealVar n("n","n",16.44);//,13,30);

//  RooRealVar p("p","p",-0.20338,-10,0.2);
  RooRealVar p("p","p",0.0);
  RooRealVar c("c","c",-20,-50,0.1);
//  RooRealVar mean_bkg("#mu_{bkg}","mean_bkg",1.6799,1.5,1.70);
//  RooRealVar sig_bkg("#sigma_{bkg}","sig_bkg",0.038,0.02,0.2);
  RooRealVar mean_bkg("#mu_{bkg}","mean_bkg",1.6799,1.5,1.90);
  RooRealVar sig_bkg("#sigma_{bkg}","sig_bkg",0.038,0.0,0.5);
 
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.9,0.68,1);
 RooRealVar f2("frac_{Bkg/Bkg3}","f3",0.9,0.56,1);

//Signal 
  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n); 
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
   Model.paramOn(xframe_1);



  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;

cout<<"num_pass"<<num_pass<<endl;

}










