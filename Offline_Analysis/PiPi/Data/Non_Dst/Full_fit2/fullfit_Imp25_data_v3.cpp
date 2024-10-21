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
//  RooRealVar dzero("dzero","dzero",1.68,2.00,"GeV");
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");


TH1D* data=new TH1D("pdat", "pdat", 100, 1.68, 2.06);
//TH1D* data=new TH1D("pdat", "pdat", 100, 1.68, 2.0);


void fullfit_Imp25_data_v3(void)
{
  //LOAD DATA FILE
//  TChain* chain=new TChain("h2");
  TChain* chain=new TChain("h2");


//  chain->Add("pipi_d4s.root");
  chain->Add("pipi_d4s.root");


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
//  for(int i=0;i<1000;i++) 
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);


if(f_dzero >1.68 && f_dzero < 2.06){       
//if(f_dzero >1.68 && f_dzero < 2.0){    

if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){
//if(f_PD > 3.5){



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
//if(f_Df !=1 && f_Df !=10){

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

  RooRealVar N_t("N_{Sig}","N_t",35000,25000,155000);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,200000000);


  //_____________________________________________________

//Signal

//Fudge Factors
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.000,-0.005,0.005);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.066);
  RooRealVar m_diff("#mu_{diff}","m_diff",0.0313);
  RooRealVar s_ratio("#sigma_{ratio}","s_ratio",3.747);
//  RooRealVar m_fudgeBkg("#mu_{fudgeBkg}","m_fudgeBkg",0.0,-0.05,0.05);
//  RooRealVar s_fudgeBkg("#sigma_{fudgeBkg}","s_fudgeBkg",1.0,0.5,2.5);


//Signal
  RooFormulaVar m("m_sig","1.86771+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.015873*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.905);  
  RooRealVar n("n","n",2.469); 
  RooFormulaVar m_g("#mu_{G}","1.86771+#mu_{fudge}+#mu_{diff}",RooArgList(m_fudge,m_diff));
  RooFormulaVar s_g("#sigma_{G}","0.015835*#sigma_{fudge}*#sigma_{ratio}",RooArgList(s_fudge,s_ratio));
  RooRealVar f_sig("f_sig","f_sig",0.9619);



  RooCBShape Sig1("Sig1", "Cystal Ball Function", dzero, m, s, a, n); 
  RooGaussian Sig2("Sig2", "Sig2",dzero,m_g,s_g);
  RooAddPdf Sig("Sig","Sig",RooArgList(Sig1,Sig2),RooArgList(f_sig));


//Background 
  RooRealVar p1("p1","p1",-0.46,-0.7,0.01);
  RooRealVar p2("p2","p2",0.08,0.0,0.15);


  RooChebychev Bkg1("Bkg1", "Bkg1", dzero, RooArgList(p1,p2));

  RooRealVar m_b("#mu_{bkg}","m_b",1.66838,1.59,1.71);
  RooRealVar s_b("#sigma_{bkg}","s_b",0.045,0.025,0.095);
  RooRealVar a_b("#alpha_bkg","a_bkg",-1.9750);
  RooRealVar n_b("n_bkg","n_bkg",1.167);



  RooCBShape Bkg2("Bkg2", "Bkg2", dzero, m_b, s_b, a_b, n_b); 
  RooRealVar f_bkg("f_bkg","f_bkg",0.93,0.92,0.995);


//Full Background 
  RooAddPdf Bkg("Bkg","Bkg",RooArgList(Bkg1,Bkg2),RooArgList(f_bkg));


//Full Model 
 RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_t,N_tb));





//Fit
   RooFitResult* fitRes = Model.fitTo(*Data,Save());




  //DeltaM PLOTING
//  RooPlot *xframe_1 =dzero.frame(Bins(38),Title("D #rightarrow #pi #pi"));
  RooPlot *xframe_1 =dzero.frame(Bins(100),Title("D #rightarrow #pi #pi"));
  Data->plotOn(xframe_1);
  Model.plotOn(xframe_1);
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame(Title("Pull Distribution")) ;
    frame1->addPlotable(hpull1,"P") ;
    frame1->SetMaximum(6);
    frame1->SetMinimum(-6);
  cout<<" signalchi2 = "<<xframe_1->chiSquare()<<endl;
 Model.plotOn(xframe_1,Components(Bkg),LineStyle(kDashed));
 Model.plotOn(xframe_1,Components(Sig),LineStyle(kDashed),LineColor(kGreen));
 Model.plotOn(xframe_1,Components(Bkg2),LineStyle(kDashed),LineColor(kRed));
   Model.paramOn(xframe_1);



  TCanvas* can3 = new TCanvas("c3","c3",600,800) ;
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();

  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;


}










