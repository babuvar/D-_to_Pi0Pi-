/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #801
// 
// A Toy Monte Carlo study that perform cycles of
// event generation and fittting
//
// 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"

using namespace RooFit ;


void Toy_bkgwidth()
{


//-----------------------------------------------------------------------------------------
   RooRealVar x("x","x",1.68,2.06,"GeV");

//DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t1("N_{Sig1}","N_t1",5400,0,150000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));


  RooRealVar N_t2("N_{Sig2}","N_t2",1750,0,60000);
  RooFormulaVar N_n2("N_n2","(0.5)*(1-A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));
  RooFormulaVar N_p2("N_p2","(0.5)*(1+A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));



//  RooRealVar Abkg1("A_{Bkg1}","ABkg1",0,-1,1);
  RooRealVar Abkg1("A_{Bkg1}","ABkg1",-0.0479,-1,1);
  RooRealVar N_tb1("N_{Bkg1}","N_tb1",21600,0,10000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));

//  RooRealVar Abkg2("A_{Bkg2}","ABkg2",0,-1,1);
  RooRealVar Abkg2("A_{Bkg2}","ABkg2",-0.0146,-1,1);
  RooRealVar N_tb2("N_{Bkg2}","N_tb2",34240,0,10000000);
  RooFormulaVar N_nb2("N_nb2","(0.5)*(1-A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));
  RooFormulaVar N_pb2("N_pb2","(0.5)*(1+A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));



  //_____________________________________________________


//Common Pars
//Signal
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0,-0.005,0.005);
  RooRealVar m_fudgeBkg("#mu_{fudgeBkg}","m_fudgeBkg",-0.0,-0.04,0.04);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.000,0.5,2.5);
  RooRealVar s_fudgeBkg("#sigma_{fudgeBkg}","s_fudgeBkg",1.000,0.5,2.5);


//Bin-1
  RooFormulaVar m("m_sig","1.86805+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01530*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.805);//,0.65,1.15);
  RooRealVar n("n","n",4.92);//,0,1000);
  RooFormulaVar m_g("m_g","1.86805+#mu_{fudge}+(0.01885*#sigma_{fudge})",RooArgList(m_fudge,s_fudge));
  RooFormulaVar s_g("s_g","0.0580*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar f_sig("f_sig","f_sig",0.9799);//,0.8,1.0);

//Bin-2
  RooFormulaVar s2("s2","0.01493*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a2("#alpha2","a2",0.781);//,0.65,1.15);
  RooRealVar n2("n2","n2",3.77);//,0,1000);
  RooFormulaVar m_g2("m_g2","1.86879+#mu_{fudge}+(0.02121*#sigma_{fudge})",RooArgList(m_fudge,s_fudge));
  RooFormulaVar s_g2("s_g2","0.0558*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar f_sig2("f_sig2","f_sig2",0.9746);//,0.8,1.0);



//Combinatorial
//Bin-1
//  RooRealVar p1("p1","p1",-0.4235,-0.8,0.01);
  RooRealVar p1("p1","p1",-0.3055,-0.8,0.2);

//Bin-2
//  RooRealVar p2("p2","p2",-0.28888,-0.8,0.01);
  RooRealVar p2("p2","p2",-0.2691,-0.8,0.2);



//Peaking Background
//Bin-1
  RooFormulaVar m_bkg("#mu_{bkg}","1.6802+#mu_{fudgeBkg}",RooArgList(m_fudgeBkg));
  RooFormulaVar s_bkg("#sigma_{bkg}","0.03656*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));
  RooRealVar a_bkg("#alpha_bkg","a_bkg",-1.9014);//,-2.2,-1.80);
  RooRealVar n_bkg("n_bkg","n_bkg",1.604);//,0.8,1.3);

//Bin-2
  RooFormulaVar m_bkg2("#mu_{bkg2}","1.6802-0.00135+#mu_{fudgeBkg}",RooArgList(m_fudgeBkg));
  RooFormulaVar s_bkg2("#sigma_{bkg2}","0.999*0.03656*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));
  RooRealVar a_bkg2("#alpha_bkg2","a_bkg2",-1.9129);//,-2.2,-1.80);
  RooRealVar n_bkg2("n_bkg2","n_bkg2",1.328);//,0.8,1.3);


//Background fraction
//Bin-1
RooRealVar f_bkg_p("f_bkg_p","f_bkg_p",0.54,0.0,1.0);
RooRealVar f_bkg_n("f_bkg_n","f_bkg_n",0.476,0.0,1.0);

//Bin-2
RooRealVar f_bkg_p2("f_bkg_p2","f_bkg_p2",0.8481,0.0,1.0);
RooRealVar f_bkg_n2("f_bkg_n2","f_bkg_n2",0.8365,0.0,1.0);


//_____________________________________________________________________________________
//			BIN - 1
//_____________________________________________________________________________________

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", x, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",x,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", x, RooArgList(p1));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", x, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg_p));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig1_n("Sig1_n", "Cystal Ball Function_n", x, m, s, a, n); 
  RooGaussian Sig2_n("Sig2_n", "Sig2_n",x,m_g,s_g);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(Sig1_n,Sig2_n),RooArgList(f_sig));


//Background - N  
//Data 
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", x, RooArgList(p1));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", x, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg_n));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


//_____________________________________________________________________________________
//			BIN - 2
//_____________________________________________________________________________________


//Signal - P
  RooCBShape Sig1_p2("Sig1_p2", "Cystal Ball Function_p2", x, m, s2, a2, n2); 
  RooGaussian Sig2_p2("Sig2_p2", "Sig2_p2",x,m_g2,s_g2);
  RooAddPdf Sig_p2("Sig_p2","Sig_p2",RooArgList(Sig1_p2,Sig2_p2),RooArgList(f_sig2));


//Background - P
  RooChebychev Bkg1_p2("Bkg1_p2", "Bkg1_p2", x, RooArgList(p2));
  RooCBShape Bkg2_p2("Bkg2_p2", "Bkg2_p2", x, m_bkg2, s_bkg2, a_bkg2, n_bkg2); 
  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(Bkg1_p2,Bkg2_p2),RooArgList(f_bkg_p2));

//Full Model - P
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N
  RooCBShape Sig1_n2("Sig1_n2", "Cystal Ball Function_n2", x, m, s2, a2, n2); 
  RooGaussian Sig2_n2("Sig2_n2", "Sig2_n2",x,m_g2,s_g2);
  RooAddPdf Sig_n2("Sig_n2","Sig_n2",RooArgList(Sig1_n2,Sig2_n2),RooArgList(f_sig2));


//Background - N  
//Data 
  RooChebychev Bkg1_n2("Bkg1_n2", "Bkg1_n2", x, RooArgList(p2));
  RooCBShape Bkg2_n2("Bkg2_n2", "Bkg2_n2", x, m_bkg2, s_bkg2, a_bkg2, n_bkg2); 
  RooAddPdf Bkg_n2("Bkg_n2","Bkg_n2",RooArgList(Bkg1_n2,Bkg2_n2),RooArgList(f_bkg_n2));

//Full Model - N
  RooAddPdf Model_n2("Model_n2","Model_n2",RooArgList(Sig_n2,Bkg_n2),RooArgList(N_n2,N_nb2));






  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");



  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  simPdf.addPdf(Model_p2,"pos2");
  simPdf.addPdf(Model_n2,"neg2");

//-----------------------------------------------------------------------------------------







  // C r e a t e   m a n a g e r
  // ---------------------------

  // Instantiate RooMCStudy manager on model with x as observable and given choice of fit options
  //
  // The Silence() option kills all messages below the PROGRESS level, leaving only a single message
  // per sample executed, and any error message that occur during fitting
  //
  // The Extended() option has two effects: 
  //    1) The extended ML term is included in the likelihood and 
  //    2) A poisson fluctuation is introduced on the number of generated events 
  //
  // The FitOptions() given here are passed to the fitting stage of each toy experiment.
  // If Save() is specified, the fit result of each experiment is saved by the manager  
  //
  // A Binned() option is added in this example to bin the data between generation and fitting
  // to speed up the study at the expemse of some precision

 
  RooMCStudy* mcstudy = new RooMCStudy(simPdf,RooArgSet(x,sample),Binned(kTRUE),Silence(kTRUE),Extended(kTRUE), FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;


  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(200) ;
//  mcstudy->generateAndFit(10000) ;

  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------


  // Make plots of the distributions of mean, the error on mean and the pull of mean


  RooPlot* frame3 = mcstudy->plotPull(s_fudgeBkg,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame1 = mcstudy->plotParam(s_fudgeBkg,Bins(40)) ;
//  RooPlot* frame2 = mcstudy->plotError(s_fudgeBkg,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(s_fudgeBkg,0.0,0.06,40) ;
/*
  RooPlot* frame3 = mcstudy->plotPull(Araw,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame1 = mcstudy->plotParam(Araw,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(Araw,Bins(40)) ;
//  RooPlot* frame2 = mcstudy->plotError(Araw,0.0,0.06,40) ;
*/

  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* can = new TCanvas("MC Study","MC Study",900,900) ;
  can->Divide(2,2) ;


  can->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  can->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;



}












