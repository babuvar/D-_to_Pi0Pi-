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
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");
  RooRealVar dzero2("dzero2","dzero2",1.7,2.0,"GeV");

  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(dzero));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(dzero));

  RooDataSet* data_P=new RooDataSet("data_P","data_P",RooArgSet(dzero2));
  RooDataSet* data_N=new RooDataSet("data_N","data_P",RooArgSet(dzero2));


void sixDfit(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  TChain* chain2=new TChain("h2");



  chain->Add("pipi_data_4s.root");
  chain->Add("pipi_data_5s.root");  
  chain->Add("pipi_data_cont.root");

  chain2->Add("pipi_data_4s.root");
  chain2->Add("pipi_data_5s.root");  
  chain2->Add("pipi_data_cont.root");




  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

  Int_t nevt2=(int)chain2->GetEntries();
  cout<<"nevt2 =\t"<<nevt2<<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
//  h1->SetBranchAddress("Dstf",&f_Dstf);
//  h1->SetBranchAddress("Df",&f_Df);
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
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){

photon1cutflag=0; photon2cutflag=0; sphoton1cutflag=0; sphoton2cutflag=0;

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



//Optimized P*D* bin
if(f_PD > 2.95){

if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));}

}


//2nd P*D* bin
if(f_PD > 2.5 && f_PD <= 2.95){

if(f_Dcharge == 1){data_p2->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n2->add(RooArgSet(dzero));}

}



//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}}//Photon cuts

}}}}}}//}//}
    }
 


  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);

  for(int i=0;i<nevt2;i++) 
    {

      chain2->GetEntry(i);
		  dzero2.setVal(f_dzero);


photon1cutflag=0; photon2cutflag=0;


if(f_dzero >1.7 && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 3.5){

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


num_pass++;
if(f_Dcharge == 1){data_P->add(RooArgSet(dzero2));}
if(f_Dcharge == -1){data_N->add(RooArgSet(dzero2));}

}//}
}}}}}
    }
 

  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t1("N_{Sig1}","N_t1",100,0,1000000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));


  RooRealVar N_t2("N_{Sig2}","N_t2",100,0,1000000);
  RooFormulaVar N_n2("N_n2","(0.5)*(1-A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));
  RooFormulaVar N_p2("N_p2","(0.5)*(1+A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));



  RooRealVar Abkg1("A_{Bkg1}","ABkg1",0,-1,1);
  RooRealVar N_tb1("N_{Bkg1}","N_tb1",100,0,1000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));

  RooRealVar Abkg2("A_{Bkg2}","ABkg2",0,-1,1);
  RooRealVar N_tb2("N_{Bkg2}","N_tb2",100,0,1000000);
  RooFormulaVar N_nb2("N_nb2","(0.5)*(1-A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));
  RooFormulaVar N_pb2("N_pb2","(0.5)*(1+A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));



  //_____________________________________________________


//Common Pars
//Signal


  RooRealVar m_fudge("mu_{fudge}","m_fudge",0.0005,-0.01,0.01);
  RooRealVar s_fudge("sigma_{fudge}","s_fudge",1.0005,0.0,2.0);
  RooFormulaVar m("m","1.86861+mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01532*sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar s2("s2","0.01558*sigma_{fudge}",RooArgList(s_fudge));

  RooRealVar a("#alpha","a",0.743);//,0,2);
  RooRealVar n("n","n",12.8);//,0,150);

//  RooRealVar a("#alpha","a",0.743,0.7,0.8);
//  RooRealVar n("n","n",12.8,5,20);

//Background
//Bin-1
  RooRealVar p("p","p",-0.5,-5,0); 
  RooRealVar mean_bkg("mean_bkg","mean_bkg",1.6799,1.5,1.70);
  RooRealVar sig_bkg("sig_bkg","sig_bkg",0.03381,0.02,0.06);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
//Bin-2
  RooRealVar p2("p2","p2",-0.5,-5,0); 
  RooRealVar sig_bkg2("sig_bkg2","sig_bkg2",0.03381,0.02,0.06);
  RooRealVar f2("frac2_{Bkg1/Bkg2}","f2",0.5,0,1); 



//--------------------------------------------------------------------------------------
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
  RooChebychev bkg1_n("bkg1_n", "bkg1_n", dzero, RooArgList(p));
  RooGaussian bkg2_n("bkg2_n", " Gaussian 2",dzero,mean_bkg,sig_bkg);
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(bkg1_n,bkg2_n),RooArgList(f));
//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));
//--------------------------------------------------------------------------------------
//Signal - P2
  RooCBShape Sig_p2("Sig_p2", "Cystal Ball Function", dzero, m, s2, a, n); 
//Background - P2 
  RooChebychev bkg1_p2("bkg1_p2", "bkg1_p2", dzero, RooArgList(p2));
  RooGaussian bkg2_p2("bkg2_p2", " Gaussian 12",dzero,mean_bkg,sig_bkg2);
//  RooGaussian bkg2_p2("bkg2_p2", " Gaussian 12",dzero,mean_bkg,sig_bkg);

  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(bkg1_p2,bkg2_p2),RooArgList(f2));
//Full Model - P2
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N2
  RooCBShape Sig_n2("Sig_n2", "Cystal Ball Function", dzero, m, s2, a, n); 
//Background - N2
  RooChebychev bkg1_n2("bkg1_n2", "bkg1_n2", dzero, RooArgList(p2));
  RooGaussian bkg2_n2("bkg2_n2", " Gaussian 22",dzero,mean_bkg,sig_bkg2);
//  RooGaussian bkg2_n2("bkg2_n2", " Gaussian 22",dzero,mean_bkg,sig_bkg);

  RooAddPdf Bkg_n2("Bkg_n2","Bkg_n2",RooArgList(bkg1_n2,bkg2_n2),RooArgList(f2));
//Full Model - N2
  RooAddPdf Model_n2("Model_n2","Model_n2",RooArgList(Sig_n2,Bkg_n2),RooArgList(N_n2,N_nb2));
//--------------------------------------------------------------------------------------



//***********************************************************************************************************************************************
  //DEFINE PDF
//Common
  RooRealVar N_T("N_{Sig3}","N_T",100,0,1000000);
  RooFormulaVar N_N("N_N","(0.5)*(1-A_{Raw})*N_{Sig3}",RooArgList(Araw,N_T));
  RooFormulaVar N_P("N_P","(0.5)*(1+A_{Raw})*N_{Sig3}",RooArgList(Araw,N_T));

  //____________________________________________________
  RooRealVar Abkg3("A_{Bkg3}","ABkg3",0,-1,1);
  RooRealVar N_TB("N_{Bkg3}","N_TB",100,0,2000000);
  RooFormulaVar N_NB("N_NB","(0.5)*(1-A_{Bkg3})*N_{Bkg3}",RooArgList(Abkg3,N_TB));
  RooFormulaVar N_PB("N_PB","(0.5)*(1+A_{Bkg3})*N_{Bkg3}",RooArgList(Abkg3,N_TB));

  //_____________________________________________________



  RooRealVar p3("p3","p3",-0.20338,-0.5,0.2); 
  RooRealVar c3("c3","c3",-25.92,-40,5);
  RooRealVar f3("frac3_{Bkg1/Bkg2}","f3",0.5,0,1); 

  RooRealVar s_fudge2("sigma_{fudge2}","s_fudge2",1.35,1.2,1.50);
  RooFormulaVar s2("s2","0.01532*sigma_{fudge2}",RooArgList(s_fudge2));
  RooRealVar a2("#alpha2","a2",0.81,0,1.5);
  RooRealVar n2("n2","n2",25,20,50);


//Signal - P
  RooCBShape Sig_p3("Sig_p3", "Cystal Ball Function", dzero2, m, s2, a2, n2); 
//Background - P  
  RooChebychev bkg1_p3("bkg1_p3", "bkg1_p3", dzero2, RooArgList(p3));
  RooExponential bkg2_p3("bkg2_p3", "bkg2_p3", dzero2, c3);
  RooAddPdf Bkg_p3("Bkg_p3","Bkg_p3",RooArgList(bkg1_p3,bkg2_p3),RooArgList(f3));
//Full Model - P
  RooAddPdf Model_p3("Model_p3","Model_p3",RooArgList(Sig_p3,Bkg_p3),RooArgList(N_P,N_PB));

//Signal - N
  RooCBShape Sig_n3("Sig_n3", "Cystal Ball Function", dzero2, m, s2, a2, n2); 
//Background - N  
  RooChebychev bkg1_n3("bkg1_n3", "bkg1_n3", dzero2, RooArgList(p3));
  RooExponential bkg2_n3("bkg2_n3", "bkg2_n3", dzero2, c3);
  RooAddPdf Bkg_n3("Bkg_n3","Bkg_n3",RooArgList(bkg1_n3,bkg2_n3),RooArgList(f3));
//Full Model - N
  RooAddPdf Model_n3("Model_n3","Model_n3",RooArgList(Sig_n3,Bkg_n3),RooArgList(N_N,N_NB));


//***********************************************************************************************************************************************









  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");
  sample.defineType("pos3");
  sample.defineType("neg3");


  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero,dzero2),Index(sample),Import("pos",*data_p),Import("neg",*data_n),Import("pos2",*data_p2),Import("neg2",*data_n2),Import("pos3",*data_P),Import("neg3",*data_N));
  


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  simPdf.addPdf(Model_p2,"pos2");
  simPdf.addPdf(Model_n2,"neg2");
  simPdf.addPdf(Model_p3,"pos3");
  simPdf.addPdf(Model_n3,"neg3");


  //_____________________________________________________

cout<<"Things are done"<<endl;


//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());






  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(76),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));




  RooPlot *xframe_2 = dzero.frame(Bins(76),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));



  RooPlot *xframe_3 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_3,Cut("sample==sample::pos2"));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Sig_p2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Bkg_p2"),ProjWData(sample,*combData),LineStyle(kDashed));



  RooPlot *xframe_4 = dzero.frame(Bins(38),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_4,Cut("sample==sample::neg2"));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Sig_n2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Bkg_n2"),ProjWData(sample,*combData),LineStyle(kDashed));



  RooPlot *xframe_5 =dzero2.frame(Bins(60),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_5,Cut("sample==sample::pos3"));
  simPdf.plotOn(xframe_5,Slice(sample,"pos3"),ProjWData(sample,*combData));
    RooPlot* frame3 = dzero2.frame(Title("Pull Distribution")) ;
  simPdf.plotOn(xframe_5,Slice(sample,"pos2"),Components("Sig_p2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_5,Slice(sample,"pos2"),Components("Bkg_p2"),ProjWData(sample,*combData),LineStyle(kDashed));



  RooPlot *xframe_6 = dzero2.frame(Bins(60),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_6,Cut("sample==sample::neg3"));
  simPdf.plotOn(xframe_6,Slice(sample,"neg3"),ProjWData(sample,*combData));
    RooPlot* frame4 = dzero2.frame(Title("Pull Distribution")) ;
  simPdf.plotOn(xframe_6,Slice(sample,"neg2"),Components("Sig_n2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_6,Slice(sample,"neg2"),Components("Bkg_n2"),ProjWData(sample,*combData),LineStyle(kDashed));


  TCanvas* can = new TCanvas("c","c") ;
  TCanvas* can2 = new TCanvas("c2","c2") ;
  TCanvas* can3 = new TCanvas("c3","c3") ;
  can->Divide(2,1) ;
  can2->Divide(2,1) ;
  can3->Divide(2,1) ;

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;

  can2->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  can2->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_4->GetYaxis()->SetTitleOffset(1.4) ; xframe_4->Draw() ;

  can3->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_5->GetYaxis()->SetTitleOffset(1.4) ; xframe_5->Draw() ;
  can3->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_6->GetYaxis()->SetTitleOffset(1.4) ; xframe_6->Draw() ;


Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;




}












