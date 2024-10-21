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
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));


void acp_MD(void)
{


  TFile *file = new TFile("KsPi_data_4s.root");


  TTree *t1 = (TTree*)file->Get("h1");
  Int_t nevt = (Int_t)t1->GetEntries();


 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Ksmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;

cout<<"done 1"<<endl;

  t1->SetBranchAddress("Deltam",&Deltam);
  t1->SetBranchAddress("Pstdst",&f_PD);
  t1->SetBranchAddress("Dmass",&f_dzero);
  t1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  t1->SetBranchAddress("Egam1s",&f_Egam1s);
  t1->SetBranchAddress("Egam2s",&f_Egam2s);
  t1->SetBranchAddress("Categ",&f_Categ);
  t1->SetBranchAddress("Ksmom",&f_Ksmom);
  t1->SetBranchAddress("Pimom",&f_Pimom);
  t1->SetBranchAddress("Dcharge",&f_Dcharge);



int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {


      t1->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
//if(f_Pizmass > 0.119  && f_Pizmass < 0.151){	//***********
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD <= 2.95){


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
if(sphoton1cutflag == 1 && sphoton2cutflag == 1){




if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));}

//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}//Photon cuts

}}}}}}//}//}

    }
 
/*
  for(int i=0;i<nevt;i++) 
    {


      t1->GetEntry(i);
		  dzero.setVal(f_dzero);

if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.


if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));}

}
    }
*/

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
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0001,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0005,0.0,3.0);
  RooFormulaVar mean("mean","1.869949+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar sig("sig","0.01307*#sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar sig1("sig1","0.004818*#sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar sig2("sig2","0.005236*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar f_Sig("f_Sig","f_Sig",0.247);



  RooRealVar c("c","c",-9,-30,10);



//Signal - P
  RooGaussian sig_p1("sig_p1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_p2("sig_p2", "signal Gaussian 2",dzero,mean,sig1,sig2);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(sig_p1,sig_p2),f_Sig);


 
//Background - P  
  RooExponential Bkg_p("Bkg_p", "Bkg_p", dzero, c);
//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//------------------------------------------------------------------------------------------------------------

//Signal - N
  RooGaussian sig_n1("sig_n1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_n2("sig_n2", "signal Gaussian 2",dzero,mean,sig1,sig2);
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
  RooPlot *xframe_1 =dzero.frame(Bins(30),Title("D^{+} #rightarrow K^{0}_{S} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
   Model_p.paramOn(xframe_1,data_p);

  RooPlot *xframe_2 = dzero.frame(Bins(30),Title("D^{-} #rightarrow K^{0}_{S} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
 //    Model_n.paramOn(xframe_2,data_n); 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;

}

