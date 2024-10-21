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
#include "TChain.h"
using namespace RooFit;
  RooRealVar dzero("dzero","dzero",1.80,1.94,"GeV");

TH1D* data_p=new TH1D("pdat", "pdat", 100, 1.8, 1.94);
TH1D* data_n=new TH1D("ndat", "ndat", 100, 1.8, 1.94);

void acp_MD(void)
{

//LOAD DATA FILE
  TChain* chain=new TChain("h2");
  chain->Add("kspi_4sd.root");
  chain->Add("kspi_5sd.root");
  chain->Add("kspi_contd.root");

//  chain->Add("kspi_4smc.root");
//  chain->Add("KsPi_MC0_5s.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Ksmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;

cout<<"done 1"<<endl;


  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Ksmom",&f_Ksmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);



int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
    {
      h2->GetEntry(i);
		  dzero.setVal(f_dzero);

if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
//if(f_PD > 2.54){
if(f_PD > 2.65){


if(f_Dcharge == 1){data_p->Fill(f_dzero);}
if(f_Dcharge == -1){data_n->Fill(f_dzero);}


}}}}

    }

RooDataHist("data_P","data_P",dzero,data_p) ; 
RooDataHist("data_N","data_N",dzero,data_n) ; 


  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",100,0,10000000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,10000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________


//Common Pars
//Signal

  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0005,0.0,3.0);

  RooRealVar mean("mean","mean",1.869949, 1.85, 1.89);
  RooRealVar mean2("mean2","mean2",1.869949, 1.85, 1.89);
  RooRealVar mean3("mean3","mean3",1.869949, 1.85, 1.89);
  RooRealVar mean4("mean4","mean4",1.869949, 1.85, 1.89);


  RooRealVar sig("sig","sig",0.01307,0.006,0.020);
  RooRealVar sigL("sigL","sigL",0.004818,0.002,0.007);
  RooRealVar sigR("sigR","sigR",0.005236,0.002,0.007);
  RooRealVar f_Sig("f_Sig","f_Sig",0.247,0.0,0.5);


  RooRealVar sig2("sig2","sig2",0.01307,0.006,0.020);
  RooRealVar sigL2("sigL2","sigL2",0.004818,0.002,0.007);
  RooRealVar sigR2("sigR2","sigR2",0.005236,0.002,0.007);
  RooRealVar f_Sig2("f_Sig2","f_Sig2",0.247,0.0,0.5);




  RooRealVar c("c","c",-9,-30,1);



//Signal - P
  RooGaussian sig_p1("sig_p1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_p2("sig_p2", "signal Gaussian 2",dzero,mean2,sigL,sigR);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(sig_p1,sig_p2),f_Sig);


 
//Background - P  
  RooExponential Bkg_p("Bkg_p", "Bkg_p", dzero, c);
//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//------------------------------------------------------------------------------------------------------------

//Signal - N
  RooGaussian sig_n1("sig_n1", "signal Gaussian 1",dzero,mean3,sig2);
  RooBifurGauss sig_n2("sig_n2", "signal Gaussian 2",dzero,mean4,sigL2,sigR2);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(sig_n1,sig_n2),f_Sig);



//Background - N  
  RooExponential Bkg_n("Bkg_n", "Bkg_n", dzero, c);
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

  //MD PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(30),Title("D^{+} #rightarrow K^{0}_{S} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
//   Model_p.paramOn(xframe_1);

  RooPlot *xframe_2 = dzero.frame(Bins(30),Title("D^{-} #rightarrow K^{0}_{S} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
//     Model_n.paramOn(xframe_2); 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; TGaxis::SetMaxDigits(3); xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; TGaxis::SetMaxDigits(3); xframe_2->Draw() ;

Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;

}










