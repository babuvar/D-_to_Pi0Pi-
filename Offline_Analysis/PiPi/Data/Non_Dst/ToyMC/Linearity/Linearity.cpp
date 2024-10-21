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
#include <TH1F.h>
#include <TF1.h>
#include <fstream>

double v_x[100], e_x[100], v_y[100], e_y[100];
int count=-1;

using namespace RooFit ;

ofstream fout;
//int n_iters=50;
int n_iters=500;
double

double araw_cent;

void Do_Test(double araw_cent);

void Linearity()
{
fout.open("Results.txt");

for (int i=-10;i<=10;i++)
//for (int i=-4;i<=4;i++)
{
araw_cent=i*0.005;
Do_Test(araw_cent);
}
fout.close();

   TGraphErrors *graph = new TGraphErrors(21,v_x,v_y,e_x,e_y);
//   graph->SetTitle("Measured A_{CP}^{uncorr}");
   graph->SetMarkerColor(9);
   graph->SetMarkerStyle(21);
   graph->Draw("AP");
//   graph->Fit("pol1");

TF1 *func = new TF1("func","[0]+[1]*x",-8,8);
//fa->SetParameter(0,value_first_parameter);
//fa->SetParameter(1,value_second_parameter);
func->SetParName(0,"Intercept");
func->SetParName(1,"Slope");
   graph->Fit(func);


}



void Do_Test(double araw_cent){

//-----------------------------------------------------------------------------------------
   RooRealVar x("x","x",1.68,2.06,"GeV");

    //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",araw_cent,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",92600,0,145000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",7260000,0,22000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________

//Common Pars
//Fudge Factors
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",0.000,-0.005,0.005);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.0);
  RooRealVar m_diff("#mu_{diff}","m_diff",0.0313);
  RooRealVar s_ratio("#sigma_{ratio}","s_ratio",3.747);



//Signal
  RooFormulaVar m("m_sig","1.86786+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.015835*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.905);  
  RooRealVar n("n","n",2.469); 
  RooFormulaVar m_g("#mu_{G}","1.86786+#mu_{fudge}+#mu_{diff}",RooArgList(m_fudge,m_diff));
  RooFormulaVar s_g("#sigma_{G}","0.015835*#sigma_{fudge}*#sigma_{ratio}",RooArgList(s_fudge,s_ratio));
  RooRealVar f_sig("f_sig","f_sig",0.9619); 

//Misreconstructed D's
  RooRealVar m_b("#mu_{bkg}","m_b",1.6625,1.62,1.71);
  RooRealVar s_b("#sigma_{bkg}","s_b",0.04198,0.02,0.095);
  RooRealVar a_b("#alpha_bkg","a_bkg",-1.9750);
  RooRealVar n_b("n_bkg","n_bkg",1.167);


//  RooFormulaVar m_b2("#mu_bkg2","1.66034+#mu_{fudgeBkg}",RooArgList(m_fudgeBkg));
//  RooFormulaVar s_b2("#sigma_bkg2","0.04173*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));
//  RooRealVar a_b2("#alpha_bkg2","a_bkg2",-2.0131);
//  RooRealVar n_b2("n_bkg2","n_bkg2",1.190);


  RooRealVar p1("p1","p1",-0.421755,-0.7,0.01);
  RooRealVar p2("p2","p2",0.06765,0.0,0.15);

  RooRealVar f_bkg2_p("f_bkg2_p","f_bkg2_p",0.9535,0.92,0.99);
  RooRealVar f_bkg2_n("f_bkg2_n","f_bkg2_n",0.9535,0.92,0.99);

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", x, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",x,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", x, RooArgList(p1,p2));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", x, m_b, s_b, a_b, n_b); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg2_p));

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
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg2_n));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));



  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");

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
  mcstudy->generateAndFit(n_iters) ;

  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------


  // Make plots of the distributions of mean, the error on mean and the pull of mean



  RooPlot* frame1 = mcstudy->plotParam(Araw,Bins(40),FitGauss(kTRUE)) ;





  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* can = new TCanvas("MC Study","MC Study",900,900) ;



  can->cd() ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;

  TCanvas* can3 = new TCanvas("MC Study3","MC Study3",900,900) ;
  can3->cd();

    TGraphAsymmErrors  *gr =  frame1->findObject("h_fitParData_simPdf");
    gr->Draw();
    gr->Fit("gaus");

fout<<"araw_cent = "<<araw_cent<<"\t Mean = "<< gr->GetFunction("gaus")->GetParameter(1)<<" ± "<<gr->GetFunction("gaus")->GetParError(1)<<endl;



count++;

v_x[count]=araw_cent*100;
e_x[count]=0.0;
v_y[count]=gr->GetFunction("gaus")->GetParameter(1)*100;
e_y[count]=gr->GetFunction("gaus")->GetParError(1)*100;



}












