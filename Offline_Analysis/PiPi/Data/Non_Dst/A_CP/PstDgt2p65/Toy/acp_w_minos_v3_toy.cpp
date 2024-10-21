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


void acp_w_minos_v3_toy()
{


//-----------------------------------------------------------------------------------------
   RooRealVar x("x","x",1.68,2.0,"GeV");

  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",76743,45000,145000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",6212635,4500000,22000000);
//  RooRealVar N_tb("N_{Bkg}","N_tb",0.0);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________

//Common Pars
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.0,-0.005,0.005);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.000);//,0.8,1.3);
  RooFormulaVar m("m_sig","1.86799+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01554*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.898);//,0.65,1.15);
  RooRealVar n("n","n",2.411);//,0,1000);
  RooFormulaVar m_g("m_g","1.86799+#mu_{fudge}+(0.02601*#sigma_{fudge})",RooArgList(m_fudge,s_fudge));
  RooFormulaVar s_g("s_g","0.0438*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar f_sig("f_sig","f_sig",0.9562);//,0.8,1.0);

//  RooRealVar p1_p("p1_p","p1_p",0.0,-1,1);
//  RooRealVar p2_p("p2_p","p2_p",0.0,-1,1);
//  RooRealVar p1_n("p1_n","p1_n",0.0,-1,1);
//  RooRealVar p2_n("p2_n","p2_n",0.0,-1,1);
  RooRealVar p1("p1","p1",-0.44936,-1,1);
  RooRealVar p2("p2","p2",0.0733,-1,1);


   RooRealVar m_b("#mu_bkg","m_b",1.6634,1.60,1.70);
  RooRealVar s_b("#sigma_bkg","s_b",0.041845,0.030,0.055);
//  RooRealVar s_b2("#sigma_bkg2","s_b2",0.041845,0.030,0.046);
  RooRealVar a_b("#alpha_bkg","a_bkg",-1.9633);//,-2.2,-1.80);
  RooRealVar n_b("n_bkg","n_bkg",1.149);//,0.8,1.3);
  RooRealVar f_bkg2("f_bkg2","f_bkg2",0.954155,0.91,0.99);
//  RooRealVar f_bkg2_p("f_bkg2_p","f_bkg2_p",0.954155,0.91,0.99);
//  RooRealVar f_bkg2_n("f_bkg2_n","f_bkg2_n",0.954155,0.91,0.99);

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", x, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",x,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", x, RooArgList(p1,p2));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", x, m_b, s_b, a_b, n_b); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg2));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig1_n("Sig1_n", "Cystal Ball Function_n", x, m, s, a, n); 
  RooGaussian Sig2_n("Sig2_n", "Sig2_n",x,m_g,s_g);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(Sig1_n,Sig2_n),RooArgList(f_sig));


//Background - N  
//Data 
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", x, RooArgList(p1,p2));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", x, m_b, s_b, a_b, n_b); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg2));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));




  //CREATE SIMULTANEOUS MODEL
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");

  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  //_____________________________________________________



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
//  mcstudy->generateAndFit(500) ;
  mcstudy->generateAndFit(10000) ;

  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------


  // Make plots of the distributions of mean, the error on mean and the pull of mean


  RooPlot* frame3 = mcstudy->plotPull(Araw,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame1 = mcstudy->plotParam(Araw,Bins(40)) ;
//  RooPlot* frame2 = mcstudy->plotError(Araw,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(Araw,0.0165,0.021,40) ;

//  RooPlot* frame6 = mcstudy->plotPull(N_t,Bins(40),FitGauss(kTRUE)) ;
//  RooPlot* frame4 = mcstudy->plotParam(N_t,Bins(40)) ;
//  RooPlot* frame5 = mcstudy->plotError(N_t,1850,1950,40) ;
//  RooPlot* frame5 = mcstudy->plotError(N_t,Bins(40)) ;



  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* can = new TCanvas("MC Study","MC Study",900,900) ;
  can->Divide(2,2) ;
//  can->Divide(3,2) ;

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  can->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;


//  can->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
//  can->cd(5) ; gPad->SetLeftMargin(0.15) ; frame5->GetYaxis()->SetTitleOffset(1.4) ; frame5->Draw() ;
//  can->cd(6) ; gPad->SetLeftMargin(0.15) ; frame6->GetYaxis()->SetTitleOffset(1.4) ; frame6->Draw() ;



}












