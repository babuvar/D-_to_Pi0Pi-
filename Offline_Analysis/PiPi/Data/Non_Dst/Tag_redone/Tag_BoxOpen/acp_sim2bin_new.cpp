//Tagged 4-Sim fitting for tagged pipi with new model
// + weights for signal and background
//Signal is Crystal-Ball alone, without Gaussian component

 
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
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));




float nsig1=0,nsig2=0,psig1=0,psig2=0;

void acp_sim2bin_new(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");

  chain->Add("../Root_Files/pipi_dcont.root"); 
  chain->Add("../Root_Files/pipi_d4s.root"); 
  chain->Add("../Root_Files/pipi_d5s.root"); 



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
//  for(int i=0;i<5000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
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



//Optimized P*D* bin
if(f_PD > 2.95){
//if(f_PD > 2.5 && f_PD <= 2.95){

if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));   }
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero)); }


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
  RooRealVar N_t("N_{Sig}","N_t",100,0,25000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig}",RooArgList(Araw,N_t));




  RooRealVar Abkg("A_{Bkg}","ABkg",0,-1,1);
  RooRealVar N_tb("N_{Bkg}","N_tb",100,0,1000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg})*N_{Bkg}",RooArgList(Abkg,N_tb));



  //_____________________________________________________


//Common Pars
//Signal
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.000635,-0.005,0.005);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.000,0.5,2.5);
  RooRealVar m_fudgeBkg("#mu_{fudgeBkg}","m_fudgeBkg",0.0,-0.04,0.04);
  RooRealVar s_fudgeBkg("#sigma_{fudgeBkg}","s_fudgeBkg",1.000,0.5,2.5);



//Bin-1
  RooFormulaVar m("m_sig","1.86736+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01677*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a("#alpha","a",0.911); 
  RooRealVar n("n","n",4.39);//,0,1000);





//Combinatorial
//Bin-1
  RooRealVar p1("p1","p1",-0.4,-0.8,0.01);




//Peaking Background
//Bin-1
  RooFormulaVar m_bkg("#mu_{bkg}","1.6802+#mu_{fudgeBkg}",RooArgList(m_fudgeBkg));
  RooFormulaVar s_bkg("#sigma_{bkg}","0.03656*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));
  RooRealVar a_bkg("#alpha_bkg","a_bkg",-1.9015);
  RooRealVar n_bkg("n_bkg","n_bkg",1.604);





//Background fraction
//Bin-1
RooRealVar f_bkg("f_bkg","f_bkg",0.5,0.0,1.0);



//_____________________________________________________________________________________
//			BIN - 1
//_____________________________________________________________________________________

//Signal - P
  RooCBShape Sig_p("Sig_p", "Cystal Ball Function_p", dzero, m, s, a, n); 



//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", dzero, RooArgList(p1));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", dzero, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig_n("Sig_n", "Cystal Ball Function_n", dzero, m, s, a, n); 


//Background - N  
//Data 
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", dzero, RooArgList(p1));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", dzero, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));







  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");

  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_p),Import("neg",*data_n));
  


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");


  //_____________________________________________________

cout<<"Things are done"<<endl;


//Fit

   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());





  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(38),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
  xframe_1->SetMaximum(1000) ;
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame(Title("Pull Distribution")) ;
    frame1->addPlotable(hpull1,"P") ;
    frame1->SetMaximum(4);
    frame1->SetMinimum(-4);
  cout<<" signalchi-1 = "<<xframe_1->chiSquare()<<endl;
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
    RooHist* hpull3 = xframe_1->residHist() ;
    RooPlot* frame3 = dzero.frame(Title("Pull Distribution")) ;
    frame3->addPlotable(hpull3,"P") ;
//  simPdf.plotOn(frame3,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
    Sig_p.plotOn(frame3,LineColor(kRed));
  Model_p.paramOn(xframe_1);



  RooPlot *xframe_2 = dzero.frame(Bins(38),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
  xframe_2->SetMaximum(1000) ;
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(4);
    frame2->SetMinimum(-4);
  cout<<" signalchi-2 = "<<xframe_2->chiSquare()<<endl;
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
    RooHist* hpull4 = xframe_2->residHist() ;
    RooPlot* frame4 = dzero.frame(Title("Pull Distribution")) ;
    frame4->addPlotable(hpull4,"P") ;
//  simPdf.plotOn(frame4,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,hpull4),LineStyle(kDashed),LineColor(kRed));
    Sig_n.plotOn(frame4,LineColor(kRed));
  TLine *line = new TLine(1.68,0.0,2.06,0.0);
//  line->SetLineColor(kRed);
  line->Draw();
  Model_n.paramOn(xframe_2);





  TCanvas* can = new TCanvas("c","c") ;
  TPad *pad11 = new TPad("pad11", "The pad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *pad12 = new TPad("pad12", "The pad 20% of the height",0.5,0.0,1.0,1.0,0);
    pad11->Draw();
    pad12->Draw();

//-------------------------------------------------------------------------------------
pad11->cd();
TPad *pad11_m = new TPad("pad31_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad11_p = new TPad("pad31_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad11_m->Draw();
    pad11_p->Draw();
pad12->cd();
TPad *pad12_m = new TPad("pad32_m", "The pad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *pad12_p = new TPad("pad32_p", "The pad 20% of the height",0.0,0.0,1.0,0.25,0);
    pad12_m->Draw();
    pad12_p->Draw();
 


/*
  pad11->cd() ; gPad->SetLeftMargin(0.25) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  pad12->cd() ; gPad->SetLeftMargin(0.25) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;
  pad21->cd() ; gPad->SetLeftMargin(0.25) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  pad22->cd() ; gPad->SetLeftMargin(0.25) ; TGaxis::SetMaxDigits(2); xframe_4->GetYaxis()->SetTitleOffset(1.4) ; xframe_4->Draw() ;
*/
  pad11_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;  line->Draw();
  pad11_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;  line->Draw();

  pad12_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;  line->Draw();
  pad12_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;  line->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   TCanvas* can2 = new TCanvas("c2","c2") ;
  TPad *mad11 = new TPad("mad11", "The mad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *mad12 = new TPad("mad12", "The mad 20% of the height",0.5,0.0,1.0,1.0,0);
    mad11->Draw();
    mad12->Draw();

//-------------------------------------------------------------------------------------
mad11->cd();
TPad *mad11_m = new TPad("mad31_m", "The mad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *mad11_p = new TPad("mad31_p", "The mad 20% of the height",0.0,0.0,1.0,0.25,0);
    mad11_m->Draw();
    mad11_p->Draw();
mad12->cd();
TPad *mad12_m = new TPad("mad32_m", "The mad 80% of the height",0.0,0.25,1.0,1.0,0);
TPad *mad12_p = new TPad("mad32_p", "The mad 20% of the height",0.0,0.0,1.0,0.25,0);
    mad12_m->Draw();
    mad12_p->Draw();


  mad11_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame3->Draw() ;  line->Draw();
  mad11_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;  line->Draw();

  mad12_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame4->Draw() ;  line->Draw();
  mad12_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;  line->Draw();


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t Asy, e_Asy, Nt, e_Nt;
Asy=Araw.getVal();
e_Asy=Araw.getError();

Nt=N_t.getVal();
e_Nt=N_t.getError();



cout<<"Araw = "<<Asy<<" +/- "<<e_Asy<<endl;

cout<<"N_t = "<<Nt<<" +/- "<<e_Nt<<endl;

 


}












