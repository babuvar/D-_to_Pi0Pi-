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
  RooRealVar dzero("dzero","dzero",1.7,2.0,"GeV");

  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));
void fit_MD3(void)
{
  //LOAD DATA FILE
//  TChain* chain=new TChain("h1");
  TChain* chain=new TChain("h1");


//  chain->Add("DtoPiZPi_exp55.root");
//  chain->Add("pipi_nonDst_GMC0_4s.root");

  chain->Add("pipi_nonDst_Data_4s.root");
  chain->Add("pipi_nonDst_Data_5s.root");
  chain->Add("pipi_nonDst_Data_cont.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge;


//  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstd",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
//  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
//  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
//  h1->SetBranchAddress("Egam1s",&f_Egam1s);
//  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
//  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
//  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
//  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<2000000;i++) 
//  for(int i=0;i<100000;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

//if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.7 && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
//if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 3.5){
//if(f_PD > 2.5 && f_PD < 5.0){
//Hard Photon 1
if(f_Egamma1 > 0.150 && f_Egamma2 > 0.150){
//if(f_Df==1){

num_pass++;
if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));}

}//}
}}}}}
    }
 
cout<<"num_pass"<<num_pass<<endl;





  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",100,0,100000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,2000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________

//Common Pars
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0005,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.35,1.3,1.4);


  RooFormulaVar m("m","1.86861+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01532*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.81,0.0,1.5);
  RooRealVar n("n","n",14,25,50);


  RooRealVar p("p","p",-0.20338,-0.5,0.1); 
  RooRealVar c("c","c",-25.92,-40,5);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 

//Signal - P
  RooCBShape Sig_p("Sig_p", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - P  
  RooChebychev bkg1_p("bkg1_p", "bkg1_p", dzero, RooArgList(p));
  RooExponential bkg2_p("bkg2_p", "bkg2_p", dzero, c);
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(bkg1_p,bkg2_p),RooArgList(f));
//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig_n("Sig_n", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - N  
  RooChebychev bkg1_n("bkg1_n", "bkg1_n", dzero, RooArgList(p));
  RooExponential bkg2_n("bkg2_n", "bkg2_n", dzero, c);
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(bkg1_n,bkg2_n),RooArgList(f));
//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_p),Import("neg",*data_n));
  


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
   Model_p.paramOn(xframe_1,data_p);

  RooPlot *xframe_2 = dzero.frame(Bins(38),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
 //    Model_n.paramOn(xframe_2,data_n); 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;

}










