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


TH1D* data_p=new TH1D("pdat", "pdat", 100, 1.68, 2.06);
TH1D* data_n=new TH1D("ndat", "ndat", 100, 1.68, 2.06);

double psig=0, msig=0;

void acp_w_minos_v3(void)
{
  //LOAD DATA FILE
//  TChain* chain=new TChain("h2");
  TChain* chain=new TChain("h2");



//chain->Add("pipi_4smc_0F.root");     
//chain->Add("pipi_4smc_1F.root");    
//chain->Add("pipi_4smc_2F.root");    
//chain->Add("pipi_4smc_3F.root");    
//chain->Add("pipi_4smc_4F.root");    
chain->Add("pipi_4smc_5F.root");




  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge;


//  h2->SetBranchAddress("Deltam",&Deltam);
  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
//  h2->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
//  h2->SetBranchAddress("Dstf",&f_Dstf);
  h2->SetBranchAddress("Df",&f_Df);
//  h2->SetBranchAddress("Egam1s",&f_Egam1s);
//  h2->SetBranchAddress("Egam2s",&f_Egam2s);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
//  h2->SetBranchAddress("Categ",&f_Categ);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
//  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
//  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


  for(int i=0;i<nevt;i++)  
//  for(int i=0;i<10000;i++)
    {
photon1cutflag=0, photon2cutflag=0;

      chain->GetEntry(i);
		  dzero.setVal(f_dzero);


//if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_dzero >1.68 && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.65){


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

float w;
if(f_Df==1 || f_Df==10){w=0.70;}
else{w=1.83;}
//else{w=0.0;}


num_pass++;
//if(f_Dcharge == 1){data_p->Fill(f_dzero,w);}
//if(f_Dcharge == -1){data_n->Fill(f_dzero,w);}

if(f_Dcharge == 1){data_p->Fill(f_dzero,w); if(f_Df==1 || f_Df==10){psig++;}   }
if(f_Dcharge == -1){data_n->Fill(f_dzero,w); if(f_Df==1 || f_Df==10){msig++;}   }

//}//for MC


}//Photon cuts


}}}}//}
    }
 
cout<<"num_pass"<<num_pass<<endl;


RooDataHist("data_P","data_P",dzero,data_p) ; 
RooDataHist("data_N","data_N",dzero,data_n) ; 

 
  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t("N_{Sig}","N_t",78000,45000,145000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));

  //____________________________________________________
  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",7000000,4500000,22000000);
//  RooRealVar N_tb("N_{Bkg}","N_tb",0.0);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));

  //_____________________________________________________

//Common Pars
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.000635,-0.005,0.005);
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
  RooRealVar p1("p1","p1",0.0,-1,1);
  RooRealVar p2("p2","p2",0.0,-1,1);


   RooRealVar m_b("#mu_bkg","m_b",1.6672,1.60,1.70);
   RooRealVar m_b2("#mu_bkg2","m_b2",1.6672,1.60,1.70);
  RooRealVar s_b("#sigma_bkg","s_b",0.04050,0.030,0.046);
  RooRealVar s_b2("#sigma_bkg2","s_b2",0.04050,0.030,0.046);
//  RooRealVar s_b("#sigma_bkg","s_b",0.04050,0.025,0.051);
  RooRealVar a_b("#alpha_bkg","a_bkg",-1.9633);//,-2.2,-1.80);
//  RooRealVar a_b2("#alpha_bkg2","a_bkg2",-2.0134);//,-2.2,-1.80);
  RooRealVar n_b("n_bkg","n_bkg",1.149);//,0.8,1.3);
//  RooRealVar n_b2("n_bkg2","n_bkg2",1.189);//,0.8,1.3);
  RooRealVar f_bkg2_p("f_bkg2_p","f_bkg2_p",0.95681,0.91,0.99);
  RooRealVar f_bkg2_n("f_bkg2_n","f_bkg2_n",0.94617,0.91,0.99);

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", dzero, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",dzero,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", dzero, RooArgList(p1,p2));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", dzero, m_b, s_b, a_b, n_b); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg2_p));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig1_n("Sig1_n", "Cystal Ball Function_n", dzero, m, s, a, n); 
  RooGaussian Sig2_n("Sig2_n", "Sig2_n",dzero,m_g,s_g);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(Sig1_n,Sig2_n),RooArgList(f_sig));


//Background - N  
//Data 
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", dzero, RooArgList(p1,p2));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", dzero, m_b, s_b2, a_b, n_b); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg2_n));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));



  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
    RooDataHist* combData =new RooDataHist("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_P),Import("neg",*data_N));


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  //_____________________________________________________


 // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = simPdf.createNLL(*combData) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  RooMinuit minu(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  minu.setVerbose(kTRUE) ;

  // Call MIGRAD to minimize the likelihood
  minu.migrad() ;

  // Print values of all paramaters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  simPdf.getParameters(dzero)->Print("s") ;

  // Disable verbose logging
  minu.setVerbose(kFALSE) ;

  // Run HESSE to calculate errors from d2L/dp2
  minu.hesse() ;

  // Print value (and error) of sigma_g2 parameter, that reflects
  // value and error back propagated from MINUIT
  Araw.Print() ;

  // Run MINOS on sigma_g2 parameter only
  minu.minos(Araw) ;

//Fit
//   RooFitResult* fitRes = simPdf.fitTo(*combData,Save(),Minos(kTRUE));



  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(100),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame(Title("Pull Distribution")) ;
    frame1->addPlotable(hpull1,"P") ;
    frame1->SetMaximum(6);
    frame1->SetMinimum(-6);
  cout<<" signalchi-2 xframe_1= "<<xframe_1->chiSquare()<<endl;
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg2_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kGreen));
//  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg1_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kViolet));
   Model_p.paramOn(xframe_1,data_P);



  RooPlot *xframe_2 = dzero.frame(Bins(100),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(6);
    frame2->SetMinimum(-6);
  cout<<" signalchi-2 xframe_2= "<<xframe_2->chiSquare()<<endl;
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg2_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kGreen));
//  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg1_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kViolet));
     Model_n.paramOn(xframe_2,data_N); 



  TCanvas* can3 = new TCanvas("c3","c3") ;
  TPad *pad31 = new TPad("pad31", "The pad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *pad32 = new TPad("pad32", "The pad 20% of the height",0.5,0.0,1.0,1.0,0);
    pad31->Draw();
    pad32->Draw();





pad31->cd();
TPad *pad31_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad31_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad31_m->Draw();
    pad31_p->Draw();
pad32->cd();
TPad *pad32_m = new TPad("pad32_m", "The pad 80% of the height",0.0,0.2,1.0,1.0,0);
TPad *pad32_p = new TPad("pad32_p", "The pad 20% of the height",0.0,0.0,1.0,0.2,0);
    pad32_m->Draw();
    pad32_p->Draw();




  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;

  pad32_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad32_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;


Double_t Asy, e_Asy;
Double_t Nt,e_Nt;

Double_t Ntot, e_Ntot;
Int_t Ntot_i, e_Ntot_i;


Asy=Araw.getVal();
e_Asy=Araw.getError();

Nt=N_t.getVal();
e_Nt=N_t.getError();

Ntot_i=Nt;


cout<<"N_t(estimated)) = "<<Ntot_i<<" +/- "<<e_Nt<<endl;
Ntot=((psig+msig)*0.7); e_Ntot=sqrt((psig+msig)*0.7);
Ntot_i=Ntot; e_Ntot_i=e_Ntot;
cout<<"N_t(truth) = "<<Ntot_i<<" +/- "<<e_Ntot_i<<endl;




cout<<"A_raw (estimated) = "<<Asy<<" +/- "<<e_Asy<<endl;
cout<<"A_raw(truth) = "<<(psig-msig)/(psig+msig)<<" +/- "<<1/sqrt((psig+msig)*0.7)<<endl;


}









