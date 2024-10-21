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


TH1D* data_p=new TH1D("pdat", "pdat", 100, 1.68, 2.06);
TH1D* data_n=new TH1D("ndat", "ndat", 100, 1.68, 2.06);

void check4(void)
{
  //LOAD DATA FILE
//  TChain* chain=new TChain("h2");
  TChain* chain=new TChain("h2");


//  chain->Add("pipi_data_4s.root");
//  chain->Add("pipi_data_5s.root");  
//  chain->Add("pipi_data_cont.root");
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
//  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
//  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


  for(int i=0;i<nevt;i++)  
//  for(int i=0;i<100000;i++)
    {
photon1cutflag=0, photon2cutflag=0;

      chain->GetEntry(i);
		  dzero.setVal(f_dzero);


if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmass > 0.122  && f_Pizmass < 0.148){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 1.31 ){
if(f_PD > 2.72){


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

//if(f_Df!=1){//for MC

num_pass++;
if(f_Dcharge == 1){data_p->Fill(f_dzero);}
if(f_Dcharge == -1){data_n->Fill(f_dzero);}

//}//for MC


}//Photon cuts


}}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;


RooDataHist("data_P","data_P",dzero,data_p) ; 
RooDataHist("data_N","data_N",dzero,data_n) ; 

  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",100000,0,1000000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",1000000,0,10000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________

//Common Pars
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.0004119);//0.0005,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.072552);//,0.0,1.7);





  RooFormulaVar m("m","1.86861+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01532*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.7937);//,0.05,0.9);
  RooRealVar n("n","n",16.44);//,13,30);

  RooRealVar p("p","p",0.0,-1,0.01); 
  RooRealVar p2("p2","p2",0.0,-1,0.01); 

  RooRealVar mean_bkg("#mu_{bkg}","mean_bkg",1.6799,1.5,1.70);
  RooRealVar mean_bkg2("#mu_{bkg2}","mean_bkg2",1.6799,1.5,1.70);

  RooRealVar sig_bkg("#sigma_{bkg}","sig_bkg",0.038,0.02,0.2);
  RooRealVar sig_bkg2("#sigma_{bkg2}","sig_bkg2",0.038,0.02,0.2);

  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.85,0.7,1.0); 
  RooRealVar f2("frac2_{Bkg1/Bkg2}","f",0.85,0.7,1.0); 

//Signal - P
  RooCBShape Sig_p("Sig_p", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - P  
  RooChebychev bkg1_p("bkg1_p", "bkg1_p", dzero, RooArgList(p));
  RooGaussian bkg2_p("bkg2_p", " Gaussian 1",dzero,mean_bkg,sig_bkg);
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(bkg1_p,bkg2_p),RooArgList(f));
//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig_n("Sig_n", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - N  
//Data
//  RooChebychev bkg1_n("bkg1_n", "bkg1_n", dzero, RooArgList(p));
//  RooGaussian bkg2_n("bkg2_n", " Gaussian 2",dzero,mean_bkg,sig_bkg);
//  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(bkg1_n,bkg2_n),RooArgList(f));
//MC
  RooChebychev bkg1_n("bkg1_n", "bkg1_n", dzero, RooArgList(p2));
  RooGaussian bkg2_n("bkg2_n", " Gaussian 2",dzero,mean_bkg2,sig_bkg2);
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(bkg1_n,bkg2_n),RooArgList(f2));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
    RooDataHist* combData =new RooDataHist("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_P),Import("neg",*data_N));


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  //_____________________________________________________




//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());


  TCanvas* can = new TCanvas("c","c") ;
  can->Divide(2,1) ;

  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
   Model_p.paramOn(xframe_1,data_P);

  RooPlot *xframe_2 = dzero.frame(Bins(38),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
//     Model_n.paramOn(xframe_2,data_N); 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;

}










