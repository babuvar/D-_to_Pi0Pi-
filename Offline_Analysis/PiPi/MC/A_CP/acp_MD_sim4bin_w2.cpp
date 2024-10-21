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
  RooRealVar dzero("dzero","dzero",1.70,2.00,"GeV");
  RooRealVar w("w","w",0.0,15.0,"no unit") ;
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero,w));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero,w));
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(dzero,w));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(dzero,w));



void acp_MD_sim4bin_w2(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");



//  chain->Add("pipi_MC6.root");
  chain->Add("pipi_MC7.root");


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
//  h1->SetBranchAddress("Dstf",&f_Dstf);
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
  h1->SetBranchAddress("Dcharge",&f_Dcharge);



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

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.70  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
//if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD <= 2.95){
//if(f_PD > 2.1 && f_PD <= 2.5){
//if(f_PD > 2.1 && f_PD <= 2.95){
//if(f_PD > 2.1){

//if(f_Categ == 5){
//if(f_Categ != 5){


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


if(f_Df==1){w.setVal(0.39);}
else{w.setVal(1.48);}
//if(f_Df==1){w.setVal(10.0);}
//w.setVal(1.0);

//Optimized P*D* bin
if(f_PD > 2.95){

if(f_Dcharge == 1){data_p->add(RooArgSet(dzero,w));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero,w));}

}


//2nd P*D* bin
if(f_PD > 2.5 && f_PD <= 2.95){

if(f_Dcharge == 1){data_p2->add(RooArgSet(dzero,w));}
if(f_Dcharge == -1){data_n2->add(RooArgSet(dzero,w));}

}



//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}}//Photon cuts

}}}}}}//}//}
    }
 




  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",100,0,100000000);
  RooRealVar M("M","M",1,0,50);
  RooFormulaVar N_p("N_p","((0.5)*(1+A_{Raw})*N_{Sig}*M)/(1+M)",RooArgList(Araw,N_t,M));
  RooFormulaVar N_p2("N_p2","((0.5)*(1+A_{Raw})*N_{Sig})/(1+M)",RooArgList(Araw,N_t,M));
  RooFormulaVar N_n("N_n","((0.5)*(1-A_{Raw})*N_{Sig}*M)/(1+M)",RooArgList(Araw,N_t,M));
  RooFormulaVar N_n2("N_n2","((0.5)*(1-A_{Raw})*N_{Sig})/(1+M)",RooArgList(Araw,N_t,M));



  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,1000000);
  RooRealVar Mb("Mb","Mb",1,0,50);
  RooFormulaVar N_pb("N_pb","((0.5)*(1+A_{Bkg})*N_{Bkg}*Mb)/(1+Mb)",RooArgList(Abkg,N_tb,Mb));
  RooFormulaVar N_pb2("N_pb2","((0.5)*(1+A_{Bkg})*N_{Bkg})/(1+Mb)",RooArgList(Abkg,N_tb,Mb));
  RooFormulaVar N_nb("N_nb","((0.5)*(1-A_{Bkg})*N_{Bkg}*Mb)/(1+Mb)",RooArgList(Abkg,N_tb,Mb));
  RooFormulaVar N_nb2("N_nb2","((0.5)*(1-A_{Bkg})*N_{Bkg})/(1+Mb)",RooArgList(Abkg,N_tb,Mb));



  //_____________________________________________________


//Common Pars
//Signal


  RooRealVar m_fudge("mu_{fudge}","m_fudge",0.0005,-0.01,0.01);
  RooRealVar s_fudge("sigma_{fudge}","s_fudge",1.0005,0.0,3.0);
  RooFormulaVar m("m","1.86806+mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01544*sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.775);//,0,2);
  RooRealVar n("n","n",10.2);//,0,150);
//Background
//Bin-1
  RooRealVar p("p","p",-0.5,-2,0); 
  RooRealVar mean_bkg("mean_bkg","mean_bkg",1.6799,1.65,1.70);
  RooRealVar sig_bkg("sig_bkg","sig_bkg",0.03381,0.03,0.045);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
//Bin-2
  RooRealVar p2("p2","p2",-0.5,-2,0); 
  RooRealVar sig_bkg2("sig_bkg2","sig_bkg2",0.03381,0.03,0.045);
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
  RooCBShape Sig_p2("Sig_p2", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - P2 
  RooChebychev bkg1_p2("bkg1_p2", "bkg1_p2", dzero, RooArgList(p2));
  RooGaussian bkg2_p2("bkg2_p2", " Gaussian 12",dzero,mean_bkg,sig_bkg2);
//  RooExponential bkg2_p2("bkg2_p2", "bkg2_p2", dzero, c2);
  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(bkg1_p2,bkg2_p2),RooArgList(f2));
//Full Model - P2
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N2
  RooCBShape Sig_n2("Sig_n2", "Cystal Ball Function", dzero, m, s, a, n); 
//Background - N2
  RooChebychev bkg1_n2("bkg1_n2", "bkg1_n2", dzero, RooArgList(p2));
  RooGaussian bkg2_n2("bkg2_n2", " Gaussian 22",dzero,mean_bkg,sig_bkg2);
//  RooExponential bkg2_n2("bkg2_n2", "bkg2_n2", dzero, c2);
  RooAddPdf Bkg_n2("Bkg_n2","Bkg_n2",RooArgList(bkg1_n2,bkg2_n2),RooArgList(f2));
//Full Model - N2
  RooAddPdf Model_n2("Model_n2","Model_n2",RooArgList(Sig_n2,Bkg_n2),RooArgList(N_n2,N_nb2));
//--------------------------------------------------------------------------------------






  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");

  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero,w),Index(sample),Import("pos",*data_p),Import("neg",*data_n),Import("pos2",*data_p2),Import("neg2",*data_n2));
  
  RooDataSet* combData_w=new RooDataSet(combData->GetName(),combData->GetTitle(),combData,*combData->get(),0,w.GetName()) ;


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  simPdf.addPdf(Model_p2,"pos2");
  simPdf.addPdf(Model_n2,"neg2");

  //_____________________________________________________

cout<<"Things are done"<<endl;


//Fit
//   RooFitResult* fitRes = simPdf.fitTo(*combData_w,Save(),SumW2Error(kTRUE));
   RooFitResult* fitRes = simPdf.fitTo(*combData_w,Save());





cout<<"Things are done-2, A_raw is "<<Araw.getVal()<<" +- "<<Araw.getError()<<endl;


  TCanvas* can = new TCanvas("c","c") ;
  TCanvas* can2 = new TCanvas("c2","c2") ;
  can->Divide(2,1) ;
  can2->Divide(2,1) ;

  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(32),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData_w->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData_w),LineStyle(kDashed));
   Model_p.paramOn(xframe_1,data_p);

  RooPlot *xframe_2 = dzero.frame(Bins(32),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData_w->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData_w),LineStyle(kDashed));

  RooPlot *xframe_3 =dzero.frame(Bins(32),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData_w->plotOn(xframe_3,Cut("sample==sample::pos2"));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Bkg_p2"),ProjWData(sample,*combData_w),LineStyle(kDashed));
   Model_p2.paramOn(xframe_3,data_p2);

  RooPlot *xframe_4 = dzero.frame(Bins(32),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData_w->plotOn(xframe_4,Cut("sample==sample::neg2"));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),ProjWData(sample,*combData_w));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Bkg_n2"),ProjWData(sample,*combData_w),LineStyle(kDashed));




  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;
  can2->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  can2->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_4->GetYaxis()->SetTitleOffset(1.4) ; xframe_4->Draw() ;







Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;





}












