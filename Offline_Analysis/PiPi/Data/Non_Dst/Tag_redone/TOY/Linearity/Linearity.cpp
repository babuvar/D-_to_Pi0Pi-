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
int n_iters=1000;
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
  RooRealVar Araw("Araw","Araw",araw_cent,-1,1);
  RooRealVar N_t1("N_{Sig1}","N_t1",5075,0,15000);
  RooFormulaVar N_n("N_n","(0.5)*(1-Araw)*N_{Sig1}",RooArgList(Araw,N_t1));
  RooFormulaVar N_p("N_p","(0.5)*(1+Araw)*N_{Sig1}",RooArgList(Araw,N_t1));


  RooRealVar N_t2("N_{Sig2}","N_t2",1685,0,15000);
  RooFormulaVar N_n2("N_n2","(0.5)*(1-Araw)*N_{Sig2}",RooArgList(Araw,N_t2));
  RooFormulaVar N_p2("N_p2","(0.5)*(1+Araw)*N_{Sig2}",RooArgList(Araw,N_t2));



  RooRealVar Abkg1("A_{Bkg1}","ABkg1",0,-1,1);
  RooRealVar N_tb1("N_{Bkg1}","N_tb1",21740,0,100000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));

  RooRealVar Abkg2("A_{Bkg2}","ABkg2",0,-1,1);
  RooRealVar N_tb2("N_{Bkg2}","N_tb2",34340,0,100000);
  RooFormulaVar N_nb2("N_nb2","(0.5)*(1-A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));
  RooFormulaVar N_pb2("N_pb2","(0.5)*(1+A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));



  //_____________________________________________________


//Common Pars
//Signal
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.0,-0.005,0.005);
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
  RooRealVar p1("p1","p1",-0.4235,-0.8,0.01);

//Bin-2
  RooRealVar p2("p2","p2",-0.28888,-0.8,0.01);



//Peaking Background
//Bin-1
   RooRealVar m_bkg("#mu_{bkg}","m_bkg",1.6802,1.62,1.74);
  RooFormulaVar s_bkg("#sigma_{bkg}","0.03656*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));


  RooRealVar a_bkg("#alpha_bkg","a_bkg",-1.8526);//,-2.2,-1.80);
  RooRealVar n_bkg("n_bkg","n_bkg",1.575);//,0.8,1.3);

//Bin-2
   RooRealVar m_bkg2("#mu_{bkg2}","m_bkg2",1.67885,1.62,1.74);
  RooFormulaVar s_bkg2("#sigma_{bkg2}","0.03652*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));

  RooRealVar a_bkg2("#alpha_bkg2","a_bkg2",-1.8661);//,-2.2,-1.80);
  RooRealVar n_bkg2("n_bkg2","n_bkg2",1.291);//,0.8,1.3);


//Background fraction
//Bin-1
RooRealVar f_bkg("f_bkg","f_bkg",0.4525,0.0,1.0);


//Bin-2
RooRealVar f_bkg2("f_bkg2","f_bkg2",0.4836,0.0,1.0);



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
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg));

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
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg));

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
  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(Bkg1_p2,Bkg2_p2),RooArgList(f_bkg2));

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
  RooAddPdf Bkg_n2("Bkg_n2","Bkg_n2",RooArgList(Bkg1_n2,Bkg2_n2),RooArgList(f_bkg2));

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












