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
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(dzero));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(dzero));


void acp_MD_sim4bin_v2()
{


//  TFile *file = new TFile("KsPi_data_4s.root");
//  TFile *file = new TFile("KsPi_data_5s.root");
//  TFile *file = new TFile("KsPi_data_4s_continuum.root");
//  TTree *t1 = (TTree*)file->Get("h1");
//  Int_t nevt = (Int_t)t1->GetEntries();

//LOAD DATA FILE
  TChain* chain=new TChain("h1");
  chain->Add("KsPi_data_4s.root");
  chain->Add("KsPi_data_5s.root");
  chain->Add("KsPi_data_4s_continuum.root");

  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;


 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Ksmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df;

cout<<"done 1"<<endl;

  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);
  h1->SetBranchAddress("Egam2s",&f_Egam2s);
  h1->SetBranchAddress("Categ",&f_Categ);
  h1->SetBranchAddress("Ksmom",&f_Ksmom);
  h1->SetBranchAddress("Pimom",&f_Pimom);
  h1->SetBranchAddress("Dcharge",&f_Dcharge);





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


      h1->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
//if(f_Pizmass > 0.119  && f_Pizmass < 0.151){	//***********
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
//if(f_Pimom > 0.84 &&  f_Pimom < 1.74){



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



if(f_PD > 2.95){
if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));}
}

if(f_PD > 2.5 && f_PD <= 2.95){
if(f_Dcharge == 1){data_p2->add(RooArgSet(dzero));}
if(f_Dcharge == -1){data_n2->add(RooArgSet(dzero));}
}

//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}//Photon cuts

}}}}}//}//}//}
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
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0001,-0.01,0.01);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0005,0.0,3.0);
//  RooRealVar s_fudge2("#sigma_{fudge2}","s_fudge2",1.0005,0.0,3.0);
//  RooRealVar s_fudge3("#sigma_{fudge3}","s_fudge3",1.0005,0.0,3.0);
//  RooRealVar s_fudge4("#sigma_{fudge4}","s_fudge4",1.0005,0.0,3.0);

  RooFormulaVar mean("mean","1.869949+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar sig("sig","0.01307*#sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar sig1("sig1","0.004818*#sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar sig2("sig2","0.005236*#sigma_{fudge}",RooArgList(s_fudge));
//  RooFormulaVar sig1("sig1","0.004818*#sigma_{fudge2}",RooArgList(s_fudge2));
//  RooFormulaVar sig2("sig2","0.005236*#sigma_{fudge2}",RooArgList(s_fudge2));


  RooFormulaVar sig_2("sig_2","0.01307*0.8914*#sigma_{fudge}",RooArgList(s_fudge));
//  RooFormulaVar sig_2("sig_2","0.01307*0.8914*#sigma_{fudge3}",RooArgList(s_fudge3));
  RooFormulaVar sig12("sig12","0.004818*0.8914*#sigma_{fudge}",RooArgList(s_fudge));
  RooFormulaVar sig22("sig22","0.005236*0.8914*#sigma_{fudge}",RooArgList(s_fudge));
//  RooFormulaVar sig12("sig12","0.004818*0.8914*#sigma_{fudge4}",RooArgList(s_fudge4));
//  RooFormulaVar sig22("sig22","0.005236*0.8914*#sigma_{fudge4}",RooArgList(s_fudge4));


//  RooRealVar f_Sig("f_Sig","f_Sig",0.247);
  RooRealVar f_Sig("f_Sig","f_Sig",0.247,0,1);

//Background
//Bin-1
  RooRealVar c1("c1","c1",-9,-30,10);
//Bin-2
  RooRealVar c2("c2","c2",-9,-30,10);



// BIN-1 --------------------------------------------------------------------------------------
//Signal - P
  RooGaussian sig_p1("sig_p1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_p2("sig_p2", "signal Gaussian 2",dzero,mean,sig1,sig2);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(sig_p1,sig_p2),f_Sig);
//Background - P  
  RooExponential Bkg_p("Bkg_p", "bkg2_p", dzero, c1);
//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooGaussian sig_n1("sig_n1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_n2("sig_n2", "signal Gaussian 2",dzero,mean,sig1,sig2);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(sig_n1,sig_n2),f_Sig);
//Background - N 
  RooExponential Bkg_n("Bkg_n", "bkg2_n", dzero, c1);
//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


// BIN-2 --------------------------------------------------------------------------------------
//Signal - P
  RooGaussian sig_p12("sig_p12", "signal Gaussian 12",dzero,mean,sig);
  RooBifurGauss sig_p22("sig_p22", "signal Gaussian 22",dzero,mean,sig1,sig2);
//  RooGaussian sig_p12("sig_p12", "signal Gaussian 12",dzero,mean,sig_2);
//  RooBifurGauss sig_p22("sig_p22", "signal Gaussian 22",dzero,mean,sig12,sig22);
  RooAddPdf Sig_p2("Sig_p2","Sig_p2",RooArgList(sig_p12,sig_p22),f_Sig);
//Background - P  
  RooExponential Bkg_p2("Bkg_p2", "bkg_p2", dzero, c2);
//Full Model - P
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N
  RooGaussian sig_n12("sig_n12", "signal Gaussian 12",dzero,mean,sig);
  RooBifurGauss sig_n22("sig_n22", "signal Gaussian 22",dzero,mean,sig1,sig2);
//  RooGaussian sig_n12("sig_n12", "signal Gaussian 12",dzero,mean,sig_2);
//  RooBifurGauss sig_n22("sig_n22", "signal Gaussian 22",dzero,mean,sig12,sig22);
  RooAddPdf Sig_n2("Sig_n2","Sig_n2",RooArgList(sig_n12,sig_n22),f_Sig);
//Background - N 
  RooExponential Bkg_n2("Bkg_n2", "bkg_n2", dzero, c2);
//Full Model - N
  RooAddPdf Model_n2("Model_n2","Model_n2",RooArgList(Sig_n2,Bkg_n2),RooArgList(N_n2,N_nb2));





  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");

  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_p),Import("neg",*data_n),Import("pos2",*data_p2),Import("neg2",*data_n2));
  


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  simPdf.addPdf(Model_p2,"pos2");
  simPdf.addPdf(Model_n2,"neg2");

  //_____________________________________________________

cout<<"Things are done"<<endl;


//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());





cout<<"Things are done-2, A_raw is "<<Araw.getVal()<<" +- "<<Araw.getError()<<endl;




  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(56),Title("D^{+} #rightarrow K^{0}_{S} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame(Title("Pull Distribution")) ;
    frame1->addPlotable(hpull1,"P") ;
    frame1->SetMaximum(6);
    frame1->SetMinimum(-6);
  cout<<" signalchi-1 = "<<xframe_1->chiSquare()<<endl;
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
//   Model_p.paramOn(xframe_1,data_p);
  xframe_1->SetMaximum(4200) ;


  RooPlot *xframe_2 = dzero.frame(Bins(56),Title("D^{-} #rightarrow K^{0}_{S} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
  cout<<" signalchi-2 = "<<xframe_2->chiSquare()<<endl;
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
  xframe_2->SetMaximum(4200) ;


  RooPlot *xframe_3 =dzero.frame(Bins(28),Title("D^{+} #rightarrow K^{0}_{S} #pi^{+}"));
  combData->plotOn(xframe_3,Cut("sample==sample::pos2"));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),ProjWData(sample,*combData));
    RooHist* hpull3 = xframe_3->pullHist() ;
    RooPlot* frame3 = dzero.frame(Title("Pull Distribution")) ;
    frame3->addPlotable(hpull3,"P") ;
    frame3->SetMaximum(5);
    frame3->SetMinimum(-5);
  cout<<" signalchi-3 = "<<xframe_3->chiSquare()<<endl;
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Sig_p2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Bkg_p2"),ProjWData(sample,*combData),LineStyle(kDashed));
//   Model_p2.paramOn(xframe_3,data_p2);
  xframe_3->SetMaximum(2300) ;


  RooPlot *xframe_4 = dzero.frame(Bins(28),Title("D^{-} #rightarrow K^{0}_{S} #pi^{-}"));
  combData->plotOn(xframe_4,Cut("sample==sample::neg2"));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),ProjWData(sample,*combData));
    RooHist* hpull4 = xframe_4->pullHist() ;
    RooPlot* frame4 = dzero.frame(Title("Pull Distribution")) ;
    frame4->addPlotable(hpull4,"P") ;
    frame4->SetMaximum(5);
    frame4->SetMinimum(-5);
  cout<<" signalchi-4 = "<<xframe_4->chiSquare()<<endl;
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Sig_n2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Bkg_n2"),ProjWData(sample,*combData),LineStyle(kDashed));
  xframe_4->SetMaximum(2300) ;


  TCanvas* can = new TCanvas("c","c") ;
  TPad *pad11 = new TPad("pad11", "The pad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *pad12 = new TPad("pad12", "The pad 20% of the height",0.5,0.0,1.0,1.0,0);
    pad11->Draw();
    pad12->Draw();
  TCanvas* can2 = new TCanvas("c2","c2") ;
  TPad *pad21 = new TPad("pad21", "The pad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *pad22 = new TPad("pad22", "The pad 20% of the height",0.5,0.0,1.0,1.0,0);
    pad21->Draw();
    pad22->Draw();
  TCanvas* can3 = new TCanvas("c3","c3",700,930) ;
  TPad *pad31 = new TPad("pad31", "The pad 80% of the height",0.0,0.5,0.5,1.0,0);
  TPad *pad32 = new TPad("pad32", "The pad 20% of the height",0.5,0.5,1.0,1.0,0);
  TPad *pad33 = new TPad("pad33", "The pad 80% of the height",0.0,0.0,0.5,0.5,0);
  TPad *pad34 = new TPad("pad34", "The pad 20% of the height",0.5,0.0,1.0,0.5,0);
    pad31->Draw();
    pad32->Draw();
    pad33->Draw();
    pad34->Draw();
pad31->cd();
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad31_m->Draw();
    pad31_p->Draw();
pad32->cd();
TPad *pad32_m = new TPad("pad32_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad32_p = new TPad("pad32_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad32_m->Draw();
    pad32_p->Draw();
pad33->cd();
TPad *pad33_m = new TPad("pad33_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad33_p = new TPad("pad33_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad33_m->Draw();
    pad33_p->Draw();
pad34->cd();
TPad *pad34_m = new TPad("pad34_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad34_p = new TPad("pad34_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad34_m->Draw();
    pad34_p->Draw();

  pad11->cd() ; gPad->SetLeftMargin(0.25) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  pad12->cd() ; gPad->SetLeftMargin(0.25) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;
  pad21->cd() ; gPad->SetLeftMargin(0.25) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  pad22->cd() ; gPad->SetLeftMargin(0.25) ; TGaxis::SetMaxDigits(2); xframe_4->GetYaxis()->SetTitleOffset(1.4) ; xframe_4->Draw() ;

  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;

  pad32_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad32_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;

  pad33_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_3->Draw() ;
  pad33_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame3->Draw() ;

  pad34_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_4->Draw() ;
  pad34_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame4->Draw() ;


Float_t Asy, e_Asy;
Asy=Araw.getVal();
e_Asy=Araw.getError();

cout<<"A_raw = "<<Asy<<" +/- "<<e_Asy<<endl;





}












