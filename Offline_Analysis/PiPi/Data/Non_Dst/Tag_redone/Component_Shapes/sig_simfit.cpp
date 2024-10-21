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
  RooRealVar dzero("dzero","dzero",1.68,2.06,"GeV");
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero));

void sig_simfit(void)
{

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
 
  chain->Add("../Root_Files/pipi_4smc_0F.root");
  chain->Add("../Root_Files/pipi_4smc_1F.root");
  chain->Add("../Root_Files/pipi_4smc_2F.root");
  chain->Add("../Root_Files/pipi_4smc_3F.root");
  chain->Add("../Root_Files/pipi_4smc_4F.root");
  chain->Add("../Root_Files/pipi_4smc_5F.root");

  chain->Add("../Root_Files/pipi_5s_gmc0.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc1.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc2.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc3.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc4.root"); 
  chain->Add("../Root_Files/pipi_5s_gmc5.root"); 


  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
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
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
//if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD < 5.0){

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
if(f_Df==1 || f_Df == 10 ){data->add(RooArgSet(dzero));}
}


if(f_PD > 2.5 && f_PD < 2.95){
if(f_Df==1 || f_Df == 10 ){data2->add(RooArgSet(dzero));}
}


num_pass++;
}}//Photon cuts

}}}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;



  //DEFINE PDF

//Common Pars 
  RooRealVar m_dif1("m_{shift(CB2,CB1)}","m_dif1",0.0,-0.005,0.005);
  RooRealVar m_dif2("m_{shift(Gaus1,CB1)}","m_dif2",0.0,-0.06,0.06);
  RooRealVar m_dif3("m_{shift(Gaus2,CB1)}","m_dif3",0.0,-0.06,0.06);
  RooRealVar s_rat1("s_{ratio(CB2,CB1)}","s_rat1",1.0,0.5,3.0);
  RooRealVar s_rat2("s_{ratio(Gaus1,CB1)}","s_rat2",1.0,0.5,5.0);
  RooRealVar s_rat3("s_{ratio(Gaus2,CB1)}","s_rat3",1.0,0.5,5.0);



//Bin - I
  RooRealVar m("m","m", 1.86861,1.86,1.88);
  RooRealVar s("s","s", 0.01532,0.005,0.03);
  RooRealVar a("#alpha","a",0.7937,0.65,2.15);
  RooRealVar n("n","n",16.44,0,1000);
  RooFormulaVar m_g("m_g","m+m_{shift(Gaus1,CB1)}",RooArgList(m,m_dif2));
  RooFormulaVar s_g("s_g","s*s_{ratio(Gaus1,CB1)}",RooArgList(s,s_rat2));
  RooRealVar f("f","f",0.99,0.8,1.0);

//Bin - II
  RooFormulaVar m2("m2","m+m_{shift(CB2,CB1)}",RooArgList(m,m_dif1));
  RooFormulaVar s2("s2","s*s_{ratio(CB2,CB1)}",RooArgList(s,s_rat1));
  RooRealVar a2("#alpha2","a2",0.7937,0.65,2.15);
  RooRealVar n2("n2","n2",16.44,0,1000);
  RooFormulaVar m2_g("m2_g","m+m_{shift(Gaus2,CB1)}",RooArgList(m,m_dif3));
  RooFormulaVar s2_g("s2_g","s*s_{ratio(Gaus2,CB1)}",RooArgList(s,s_rat3));
  RooRealVar f2("f2","f2",0.99,0.8,1.0);


//Signal  - I
  RooCBShape sig1_a("sig1_a", "Cystal Ball Function", dzero, m, s, a, n); 
  RooGaussian sig1_b("sig1_b", "sig2_b",dzero,m_g,s_g);
  RooAddPdf Sig1("Sig1","Sig1",RooArgList(sig1_a,sig1_b),RooArgList(f));


//Signal  - II
  RooCBShape sig2_a("sig2_a", "Cystal Ball Function", dzero, m2, s2, a2, n2); 
  RooGaussian sig2_b("sig2_b", "sig2_b",dzero,m2_g,s2_g);
  RooAddPdf Sig2("Sig2","Sig2",RooArgList(sig2_a,sig2_b),RooArgList(f2));


  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("bin1");
  sample.defineType("bin2");


  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("bin1",*data),Import("bin2",*data2));
  


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Sig1,"bin1");
  simPdf.addPdf(Sig2,"bin2");



//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());




  RooPlot *xframe_1 =dzero.frame(Bins(25));
  xframe_1->SetTitle("Bin - 1");
  data->plotOn(xframe_1);
    Sig1.plotOn(xframe_1, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
    Sig1.paramOn(xframe_1);
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame() ;
    frame1->addPlotable(hpull1,"P") ;
   frame1->SetMaximum(6);
    frame1->SetMinimum(-6);
  cout<<" signalchi-1 xframe_1= "<<xframe_1->chiSquare()<<endl;


 
  RooPlot *xframe_2 =dzero.frame(Bins(25));
  xframe_2->SetTitle("Bin - 2");
  data2->plotOn(xframe_2);
    Sig2.plotOn(xframe_2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
    Sig2.paramOn(xframe_2);
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame() ;
    frame2->addPlotable(hpull2,"P") ;
   frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;



  TCanvas* can = new TCanvas("c","c",700,700) ;
TPad *pad1_m = new TPad("pad1_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad1_p = new TPad("pad1_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad1_m->Draw();
    pad1_p->Draw();

  pad1_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;
  pad1_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;


  TCanvas* can2 = new TCanvas("c2","c2",700,700) ;
TPad *pad2_m = new TPad("pad2_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad2_p = new TPad("pad2_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad2_m->Draw();
    pad2_p->Draw();

  pad2_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad2_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;


}



















