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
  RooRealVar dzero("dzero","dzero",1.68,2.00,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero));
void fit_MD2_w(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");


  chain->Add("pipi_MC0_Y4s.root");
  chain->Add("pipi_MC0_Y5s.root");

//  chain->Add("pipi_data_Y4s.root");
//  chain->Add("pipi_data_Y5s.root");
//  chain->Add("pipi_data_continuum.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Egamma1",&f_Egamma1);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Pizmom",&f_Pizmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h1->SetBranchAddress("Gam2hthe",&f_Gam2thet);

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
if(f_dzero >1.68  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
//if(f_PD > 2.95){


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

if(f_PD > 2.95){
data->add(RooArgSet(dzero));
}
if(f_PD > 2.5 && f_PD <= 2.95){
data2->add(RooArgSet(dzero));
}

num_pass++;
}}//Photon cuts

}}}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;



  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,1000000);
  RooRealVar N_sig2("N_{Sig2}","N_sig2",100,0,1000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,1000000);
  RooRealVar N_bkg2("N_{Bkg2}","N_bkg2",100,0,1000000);

  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0005,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0005,0.0,3.0);
  RooFormulaVar m("m","1.86806+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01544*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.775);//,0,2);
  RooRealVar n("n","n",10.2);//,0,150);

  RooRealVar p("p","p",-0.5,-2,0);  
  RooRealVar mean_bkg("#mu_{bkg}","mean_bkg",1.6799,1.65,1.70);
  RooRealVar sig_bkg("#sigma_{bkg}","sig_bkg",0.03381,0.03,0.04);
  RooRealVar f("frac_{Bkg1/Bkg2}","f",0.5,0,1); 
  RooRealVar p2("p2","p2",-0.5,-2,0);  
  RooRealVar mean_bkg2("#mu_{bkg2}","mean_bkg2",1.6799,1.65,1.70);
  RooRealVar sig_bkg2("#sigma_{bkg2}","sig_bkg2",0.03381,0.03,0.04);
  RooRealVar f2("frac2_{Bkg1/Bkg2}","f2",0.5,0,1); 

//Signal
  RooCBShape Sig("Sig", "Cystal Ball Function", dzero, m, s, a, n); 

//Signal-2
  RooCBShape Sig2("Sig2", "Cystal Ball Function", dzero, m, s, a, n); 


//Background
  RooChebychev bkg1("bkg1", "bkg1", dzero, RooArgList(p));
  RooGaussian bkg2("bkg2", " Gaussian 1",dzero,mean_bkg,sig_bkg);
  RooAddPdf Bkg("Bkg","Bkg",RooArgList(bkg1,bkg2),RooArgList(f));

//Background-2
  RooChebychev bkg1_2("bkg1_2", "bkg1_2", dzero, RooArgList(p2));
  RooGaussian bkg2_2("bkg2_2", " Gaussian 1_2",dzero,mean_bkg2,sig_bkg2);
  RooAddPdf Bkg2("Bkg2","Bkg2",RooArgList(bkg1_2,bkg2_2),RooArgList(f2));




//Full Model
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));

//Full Model-2
  RooAddPdf Model2("Model2","Model2",RooArgList(Sig2,Bkg2),RooArgList(N_sig2,N_bkg2));


  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("bin1");
  sample.defineType("bin2");

  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("bin",*data),Import("bin2",*data2));

  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model,"bin");
  simPdf.addPdf(Model2,"bin2");

//Fit
//   RooFitResult* fitRes = Model.fitTo(*data,Save());
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());


  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(32),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::bin"));
  simPdf.plotOn(xframe_1,Slice(sample,"bin"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_1,Slice(sample,"bin"),Components("Bkg"),ProjWData(sample,*combData),LineStyle(kDashed));
   Model.paramOn(xframe_1,data);



  RooPlot *xframe_2 = dzero.frame(Bins(32),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::bin2"));
  simPdf.plotOn(xframe_2,Slice(sample,"bin2"),ProjWData(sample,*combData));
  simPdf.plotOn(xframe_2,Slice(sample,"bin2"),Components("Bkg2"),ProjWData(sample,*combData),LineStyle(kDashed));
   Model2.paramOn(xframe_2,data2);



  TCanvas* can = new TCanvas("c","c",700,800) ;
  TCanvas* can2 = new TCanvas("c2","c2",700,800) ;

  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can2->cd() ; gPad->SetLeftMargin(0.15) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;










}

