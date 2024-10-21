//Tagged 4-Sim fitting for tagged pipi with new model
// + weights for signal and background
// + binned chi2 fit
// + untagged sample
 
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
  
TH1D* data_p=new TH1D("pdat", "pdat", 100, 1.68, 2.06);
TH1D* data_n=new TH1D("ndat", "ndat", 100, 1.68, 2.06);

TH1D* data_p2=new TH1D("pdat2", "pdat2", 100, 1.68, 2.06);
TH1D* data_n2=new TH1D("ndat2", "ndat2", 100, 1.68, 2.06);

TH1D* data_p3=new TH1D("pdat3", "pdat3", 100, 1.68, 2.06);
TH1D* data_n3=new TH1D("ndat3", "ndat3", 100, 1.68, 2.06);



float nsig1=0,nsig2=0,psig1=0,psig2=0;

void acp_MD_sim6bin_alt_v5_MC(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  TChain* chain2=new TChain("h2");


//   chain->Add("pipi_4smc_0F.root");
  chain->Add("pipi_4smc_2F.root"); 
//  chain->Add("pipi_4smc_4F.root");
//  chain->Add("pipi_4smc_1F.root"); 
//  chain->Add("pipi_4smc_3F.root"); 
//  chain->Add("pipi_4smc_5F.root");

//  chain->Add("pipi_5s_gmc0.root"); 
  chain->Add("pipi_5s_gmc2.root"); 
//  chain->Add("pipi_5s_gmc4.root");
//  chain->Add("pipi_5s_gmc1.root"); 
//  chain->Add("pipi_5s_gmc3.root"); 
//  chain->Add("pipi_5s_gmc5.root");

//   chain2->Add("pipi_4smc_0F.root");
  chain2->Add("pipi_4smc_2F.root"); 
//  chain2->Add("pipi_4smc_4F.root");
//  chain2->Add("pipi_4smc_1F.root"); 
//  chain2->Add("pipi_4smc_3F.root"); 
//  chain2->Add("pipi_4smc_5F.root");

//  chain2->Add("pipi_5s_gmc0.root"); 
  chain2->Add("pipi_5s_gmc2.root"); 
//  chain2->Add("pipi_5s_gmc4.root");
//  chain2->Add("pipi_5s_gmc1.root"); 
//  chain2->Add("pipi_5s_gmc3.root"); 
//  chain2->Add("pipi_5s_gmc5.root");



  Int_t nevt=(int)chain->GetEntries();
  Int_t nevt2=(int)chain2->GetEntries();
//  cout<<"nevt\t"<<nevt <<endl;

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


double w;

//Optimized P*D* bin
if(f_PD > 2.95){
if(f_Df==1 || f_Df==10){w=0.4132;}
else{w=1.448;}

if(f_Dcharge == 1){data_p->Fill(f_dzero,w); if(f_Df == 1 || f_Df == 10){psig1=psig1+0.4132;} }
if(f_Dcharge == -1){data_n->Fill(f_dzero,w); if(f_Df == 1 || f_Df == 10){nsig1=nsig1+0.4132;}}

}


//2nd P*D* bin
if(f_PD > 2.5 && f_PD <= 2.95){
if(f_Df==1 || f_Df==10){w=0.3050;}
else{w=1.292;}

if(f_Dcharge == 1){data_p2->Fill(f_dzero,w); if(f_Df == 1 || f_Df == 10){psig2=psig2+0.3050;} }
if(f_Dcharge == -1){data_n2->Fill(f_dzero,w); if(f_Df == 1 || f_Df == 10){nsig2=nsig2+0.3050;}}

}



//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}}//Photon cuts

}}}}}}//}//}
    }
 
RooDataHist("data_P","data_P",dzero,data_p) ; 
RooDataHist("data_N","data_N",dzero,data_n) ; 
RooDataHist("data_P2","data_P2",dzero,data_p2) ; 
RooDataHist("data_N2","data_N2",dzero,data_n2) ;

// Untagged sample
  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Df",&f_Df);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);

 photon1cutflag=0; photon2cutflag=0; sphoton1cutflag=0; sphoton2cutflag=0;




  for(int i=0;i<nevt2;i++)  
//  for(int i=0;i<10000;i++)
    {
photon1cutflag=0, photon2cutflag=0;

      chain2->GetEntry(i);
		  dzero.setVal(f_dzero);


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

double w;
if(f_Df==1 || f_Df==10){w=0.70;}
else{w=1.83;}



num_pass++;
if(f_Dcharge == 1){data_p3->Fill(f_dzero,w);}
if(f_Dcharge == -1){data_n3->Fill(f_dzero,w);}

//if(f_Dcharge == 1){data_p->Fill(f_dzero,w); if(f_Df==1 || f_Df==10){psig++;}   }
//if(f_Dcharge == -1){data_n->Fill(f_dzero,w); if(f_Df==1 || f_Df==10){msig++;}   }

//}//for MC


}//Photon cuts


}}}}//}
    }
 


RooDataHist("data_P3","data_P3",dzero,data_p3) ; 
RooDataHist("data_N3","data_N3",dzero,data_n3) ; 














  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t1("N_{Sig1}","N_t1",100,0,15000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));


  RooRealVar N_t2("N_{Sig2}","N_t2",100,0,6000);
  RooFormulaVar N_n2("N_n2","(0.5)*(1-A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));
  RooFormulaVar N_p2("N_p2","(0.5)*(1+A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));


  RooRealVar N_t3("N_{Sig3}","N_t3",78000,45000,145000);
  RooFormulaVar N_n3("N_n3","(0.5)*(1-A_{Raw})*N_{Sig3}",RooArgList(Araw,N_t3));
  RooFormulaVar N_p3("N_p3","(0.5)*(1+A_{Raw})*N_{Sig3}",RooArgList(Araw,N_t3));




  RooRealVar Abkg1("A_{Bkg1}","ABkg1",0,-1,1);
  RooRealVar N_tb1("N_{Bkg1}","N_tb1",100,0,1000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));

  RooRealVar Abkg2("A_{Bkg2}","ABkg2",0,-1,1);
  RooRealVar N_tb2("N_{Bkg2}","N_tb2",100,0,1000000);
  RooFormulaVar N_nb2("N_nb2","(0.5)*(1-A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));
  RooFormulaVar N_pb2("N_pb2","(0.5)*(1+A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));

  RooRealVar Abkg3("A_{Bkg3}","ABkg3",0,-1,1);
  RooRealVar N_tb3("N_{Bkg3}","N_tb3",7000000,4500000,22000000);
  RooFormulaVar N_nb3("N_nb3","(0.5)*(1-A_{Bkg3})*N_{Bkg3}",RooArgList(Abkg3,N_tb3));
  RooFormulaVar N_pb3("N_pb3","(0.5)*(1+A_{Bkg3})*N_{Bkg3}",RooArgList(Abkg3,N_tb3));

  //_____________________________________________________


//Common Pars
//Signal
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.000635,-0.005,0.005);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.000,0.95,1.45);

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
  RooRealVar n2("n2","n2",3.77);//,0,1000);


//Bin-3
  RooFormulaVar s3("s3","0.01554*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar a3("#alpha3","a3",0.898);//,0.65,1.15);
  RooRealVar n3("n3","n3",2.411);//,0,1000);
  RooFormulaVar m_g3("m_g3","1.86805+#mu_{fudge}+(0.02595*#sigma_{fudge})",RooArgList(m_fudge,s_fudge));
  RooFormulaVar s_g3("s_g3","0.0438*#sigma_{fudge}",RooArgList(s_fudge));
  RooRealVar f_sig3("f_sig3","f_sig3",0.9562);//,0.8,1.0);

//Combinatorial
//Bin-1
  RooRealVar p1_p("p1_p","p1_p",-0.4,-0.6,-0.1);
  RooRealVar p1_n("p1_n","p1_n",-0.4,-0.6,-0.1);
  RooRealVar p2_p("p2_p","p2_p",0.06,0.01,0.2);
  RooRealVar p2_n("p2_n","p2_n",0.06,0.01,0.2);
//Bin-2
  RooRealVar p1_p2("p1_p2","p1_p2",-0.3,-0.6,-0.1);
  RooRealVar p1_n2("p1_n2","p1_n2",-0.3,-0.6,-0.1);
  RooRealVar p2_p2("p2_p2","p2_p2",0.03,0.01,0.06);
  RooRealVar p2_n2("p2_n2","p2_n2",0.03,0.01,0.06);

//Bin-3
  RooRealVar p1_p3("p1_p3","p1_p3",-0.3,-1.0,1.0);
  RooRealVar p1_n3("p1_n3","p1_n3",-0.3,-1.0,1.0);
  RooRealVar p2_p3("p2_p3","p2_p3",0.03,-1.0,1.0);
  RooRealVar p2_n3("p2_n3","p2_n3",0.03,-1.0,1.0);


//Peaking Background
//Bin-1
   RooRealVar m_bkg_p("#mu_{bkg_p}","m_bkg_p",1.6672,1.63,1.74);
   RooRealVar m_bkg_n("#mu_{bkg_n}","m_bkg_n",1.6672,1.63,1.74);
  RooRealVar s_bkg_p("#sigma_{bkg_p}","s_bkg_p",0.04050,0.030,0.046);
  RooRealVar s_bkg_n("#sigma_{bkg_n}","s_bkg_n",0.04050,0.030,0.046);
  RooRealVar a_bkg("#alpha_bkg","a_bkg",-1.9014);//,-2.2,-1.80);
  RooRealVar n_bkg("n_bkg","n_bkg",1.604);//,0.8,1.3);

//Bin-2
  RooRealVar a_bkg2("#alpha_bkg2","a_bkg2",-1.9129);//,-2.2,-1.80);
  RooRealVar n_bkg2("n_bkg2","n_bkg2",1.328);//,0.8,1.3);

//Bin-3
   RooRealVar m_bkg_p3("#mu_{bkg_p3}","m_bkg_p3",1.6672,1.58,1.74);
   RooRealVar m_bkg_n3("#mu_{bkg_n3}","m_bkg_n3",1.6672,1.58,1.74);
  RooRealVar s_bkg_p3("#sigma_{bkg_p3}","s_bkg_p3",0.04050,0.030,0.046);
  RooRealVar s_bkg_n3("#sigma_{bkg_n3}","s_bkg_n3",0.04050,0.030,0.046);
  RooRealVar n_bkg3("n_bkg3","n_bkg3",1.104);//,0.8,1.3);




//Background fraction
//Bin-1
RooRealVar f_bkg_p("f_bkg_p","f_bkg_p",0.5,0.0,1.0);
RooRealVar f_bkg_n("f_bkg_n","f_bkg_n",0.5,0.0,1.0);

//Bin-2
RooRealVar f_bkg_p2("f_bkg_p2","f_bkg_p2",0.5,0.0,1.0);
RooRealVar f_bkg_n2("f_bkg_n2","f_bkg_n2",0.5,0.0,1.0);

//Bin-3
RooRealVar f_bkg_p3("f_bkg_p3","f_bkg_p3",0.96,0.91,0.99);
RooRealVar f_bkg_n3("f_bkg_n3","f_bkg_n3",0.95,0.91,0.99);


//_____________________________________________________________________________________
//			BIN - 1
//_____________________________________________________________________________________

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", dzero, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",dzero,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
//  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", dzero, RooArgList(p1_p,p2_p));
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", dzero, RooArgList(p1_p));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", dzero, m_bkg_p, s_bkg_p, a_bkg, n_bkg); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg_p));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig1_n("Sig1_n", "Cystal Ball Function_n", dzero, m, s, a, n); 
  RooGaussian Sig2_n("Sig2_n", "Sig2_n",dzero,m_g,s_g);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(Sig1_n,Sig2_n),RooArgList(f_sig));


//Background - N  
//Data 
//  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", dzero, RooArgList(p1_n,p2_n));
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", dzero, RooArgList(p1_n));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", dzero, m_bkg_n, s_bkg_n, a_bkg, n_bkg); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg_n));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


//_____________________________________________________________________________________
//			BIN - 2
//_____________________________________________________________________________________


//Signal - P
  RooCBShape Sig1_p2("Sig1_p2", "Cystal Ball Function_p2", dzero, m, s2, a, n2); 
  RooGaussian Sig2_p2("Sig2_p2", "Sig2_p2",dzero,m_g,s_g);
  RooAddPdf Sig_p2("Sig_p2","Sig_p2",RooArgList(Sig1_p2,Sig2_p2),RooArgList(f_sig));


//Background - P
//  RooChebychev Bkg1_p2("Bkg1_p2", "Bkg1_p2", dzero, RooArgList(p1_p2,p2_p2));
  RooChebychev Bkg1_p2("Bkg1_p2", "Bkg1_p2", dzero, RooArgList(p1_p2));
  RooCBShape Bkg2_p2("Bkg2_p2", "Bkg2_p2", dzero, m_bkg_p, s_bkg_p, a_bkg2, n_bkg2); 
  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(Bkg1_p2,Bkg2_p2),RooArgList(f_bkg_p2));

//Full Model - P
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N
  RooCBShape Sig1_n2("Sig1_n2", "Cystal Ball Function_n2", dzero, m, s2, a, n2); 
  RooGaussian Sig2_n2("Sig2_n2", "Sig2_n2",dzero,m_g,s_g);
  RooAddPdf Sig_n2("Sig_n2","Sig_n2",RooArgList(Sig1_n2,Sig2_n2),RooArgList(f_sig));


//Background - N  
//Data 
//  RooChebychev Bkg1_n2("Bkg1_n2", "Bkg1_n2", dzero, RooArgList(p1_n2,p2_n2));
  RooChebychev Bkg1_n2("Bkg1_n2", "Bkg1_n2", dzero, RooArgList(p1_n2));
  RooCBShape Bkg2_n2("Bkg2_n2", "Bkg2_n2", dzero, m_bkg_n, s_bkg_n, a_bkg2, n_bkg2); 
  RooAddPdf Bkg_n2("Bkg_n2","Bkg_n2",RooArgList(Bkg1_n2,Bkg2_n2),RooArgList(f_bkg_n2));

//Full Model - N
  RooAddPdf Model_n2("Model_n2","Model_n2",RooArgList(Sig_n2,Bkg_n2),RooArgList(N_n2,N_nb2));


//_____________________________________________________________________________________
//			BIN - 3
//_____________________________________________________________________________________

//Signal - P
  RooCBShape Sig1_p3("Sig1_p3", "Cystal Ball Function_p", dzero, m, s3, a3, n3); 
  RooGaussian Sig2_p3("Sig2_p3", "Sig2_p3",dzero,m_g3,s_g3);
  RooAddPdf Sig_p3("Sig_p3","Sig_p3",RooArgList(Sig1_p3,Sig2_p3),RooArgList(f_sig3));


//Background - P
  RooChebychev Bkg1_p3("Bkg1_p3", "Bkg1_p3", dzero, RooArgList(p1_p3,p2_p3));
  RooCBShape Bkg2_p3("Bkg2_p3", "Bkg2_p3", dzero, m_bkg_p3, s_bkg_p3, a_bkg, n_bkg3); 
  RooAddPdf Bkg_p3("Bkg_p3","Bkg_p3",RooArgList(Bkg1_p3,Bkg2_p3),RooArgList(f_bkg_p3));

//Full Model - P
  RooAddPdf Model_p3("Model_p3","Model_p3",RooArgList(Sig_p3,Bkg_p3),RooArgList(N_p3,N_pb3));

//Signal - N
  RooCBShape Sig1_n3("Sig1_n3", "Cystal Ball Function_n", dzero, m, s3, a3, n3); 
  RooGaussian Sig2_n3("Sig2_n3", "Sig2_n3",dzero,m_g3,s_g3);
  RooAddPdf Sig_n3("Sig_n3","Sig_n3",RooArgList(Sig1_n3,Sig2_n3),RooArgList(f_sig3));


//Background - N  
//Data 
  RooChebychev Bkg1_n3("Bkg1_n3", "Bkg1_n3", dzero, RooArgList(p1_n3,p2_n3));
  RooCBShape Bkg2_n3("Bkg2_n3", "Bkg2_n3", dzero, m_bkg_n3, s_bkg_n3, a_bkg, n_bkg3);  
  RooAddPdf Bkg_n3("Bkg_n3","Bkg_n3",RooArgList(Bkg1_n3,Bkg2_n3),RooArgList(f_bkg_n3));

//Full Model - N
  RooAddPdf Model_n3("Model_n3","Model_n3",RooArgList(Sig_n3,Bkg_n3),RooArgList(N_n3,N_nb3));



//----------------------------------------------------------------------------------------------------












  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("pos");
  sample.defineType("neg");
  sample.defineType("pos2");
  sample.defineType("neg2");
  sample.defineType("pos3");
  sample.defineType("neg3");


//  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero,w),Index(sample),Import("pos",*data_p),Import("neg",*data_n),Import("pos2",*data_p2),Import("neg2",*data_n2));

      RooDataHist* combData =new RooDataHist("combData","combData",RooArgSet(dzero),Index(sample),Import("pos",*data_P),Import("neg",*data_N),Import("pos2",*data_P2),Import("neg2",*data_N2),Import("pos3",*data_P3),Import("neg3",*data_N3));




  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model_p,"pos");
  simPdf.addPdf(Model_n,"neg");
  simPdf.addPdf(Model_p2,"pos2");
  simPdf.addPdf(Model_n2,"neg2");
  simPdf.addPdf(Model_p3,"pos3");
  simPdf.addPdf(Model_n3,"neg3");


  //_____________________________________________________

cout<<"Things are done"<<endl;


//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());






  //DeltaM PLOTING
  RooPlot *xframe_1 =dzero.frame(Bins(100),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_1,Cut("sample==sample::pos"));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),ProjWData(sample,*combData));
//  xframe_1->SetMaximum(520) ;
    RooHist* hpull1 = xframe_1->pullHist() ;
    RooPlot* frame1 = dzero.frame(Title("Pull Distribution")) ;
    frame1->addPlotable(hpull1,"P") ;
    frame1->SetMaximum(4);
    frame1->SetMinimum(-4);
  cout<<" signalchi-1 = "<<xframe_1->chiSquare()<<endl;
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Sig_p"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_1,Slice(sample,"pos"),Components("Bkg_p"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_p.paramOn(xframe_1);



  RooPlot *xframe_2 = dzero.frame(Bins(100),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_2,Cut("sample==sample::neg"));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),ProjWData(sample,*combData));
//  xframe_2->SetMaximum(520) ;
    RooHist* hpull2 = xframe_2->pullHist() ;
    RooPlot* frame2 = dzero.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(hpull2,"P") ;
    frame2->SetMaximum(4);
    frame2->SetMinimum(-4);
  cout<<" signalchi-2 = "<<xframe_2->chiSquare()<<endl;
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Sig_n"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_2,Slice(sample,"neg"),Components("Bkg_n"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_n.paramOn(xframe_2);


  RooPlot *xframe_3 =dzero.frame(Bins(100),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_3,Cut("sample==sample::pos2"));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),ProjWData(sample,*combData));
//  xframe_3->SetMaximum(820) ;
    RooHist* hpull3 = xframe_3->pullHist() ;
    RooPlot* frame3 = dzero.frame(Title("Pull Distribution")) ;
    frame3->addPlotable(hpull3,"P") ;
    frame3->SetMaximum(4);
    frame3->SetMinimum(-4);
  cout<<" signalchi-3 = "<<xframe_3->chiSquare()<<endl;
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Sig_p2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_3,Slice(sample,"pos2"),Components("Bkg_p2"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_p2.paramOn(xframe_3);


  RooPlot *xframe_4 = dzero.frame(Bins(100),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_4,Cut("sample==sample::neg2"));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),ProjWData(sample,*combData));
//  xframe_4->SetMaximum(820);
    RooHist* hpull4 = xframe_4->pullHist() ;
    RooPlot* frame4 = dzero.frame(Title("Pull Distribution")) ;
    frame4->addPlotable(hpull4,"P") ;
    frame4->SetMaximum(4);
    frame4->SetMinimum(-4);
  cout<<" signalchi-4 = "<<xframe_4->chiSquare()<<endl;
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Sig_n2"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_4,Slice(sample,"neg2"),Components("Bkg_n2"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_n2.paramOn(xframe_4);


  RooPlot *xframe_5 =dzero.frame(Bins(100),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  combData->plotOn(xframe_5,Cut("sample==sample::pos3"));
  simPdf.plotOn(xframe_5,Slice(sample,"pos3"),ProjWData(sample,*combData));
    RooHist* hpull5 = xframe_5->pullHist() ;
    RooPlot* frame5 = dzero.frame(Title("Pull Distribution")) ;
    frame5->addPlotable(hpull5,"P") ;
    frame5->SetMaximum(6);
    frame5->SetMinimum(-6);
  cout<<" signalchi-5 = "<<xframe_5->chiSquare()<<endl;
  simPdf.plotOn(xframe_5,Slice(sample,"pos3"),Components("Sig_p3"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_5,Slice(sample,"pos3"),Components("Bkg_p3"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_p3.paramOn(xframe_5);


  RooPlot *xframe_6 = dzero.frame(Bins(100),Title("D^{-} #rightarrow #pi^{0} #pi^{-}"));
  combData->plotOn(xframe_6,Cut("sample==sample::neg3"));
  simPdf.plotOn(xframe_6,Slice(sample,"neg3"),ProjWData(sample,*combData));
    RooHist* hpull6 = xframe_6->pullHist() ;
    RooPlot* frame6 = dzero.frame(Title("Pull Distribution")) ;
    frame6->addPlotable(hpull6,"P") ;
    frame6->SetMaximum(6);
    frame6->SetMinimum(-6);
  cout<<" signalchi-6 = "<<xframe_6->chiSquare()<<endl;
  simPdf.plotOn(xframe_6,Slice(sample,"neg3"),Components("Sig_n3"),ProjWData(sample,*combData),LineStyle(kDashed),LineColor(kRed));
  simPdf.plotOn(xframe_6,Slice(sample,"neg3"),Components("Bkg_n3"),ProjWData(sample,*combData),LineStyle(kDashed));
  Model_n3.paramOn(xframe_6);








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
  TCanvas* can4 = new TCanvas("c4","c4") ;
  TPad *pad41 = new TPad("pad41", "The pad 80% of the height",0.0,0.0,0.5,1.0,0);
  TPad *pad42 = new TPad("pad42", "The pad 20% of the height",0.5,0.0,1.0,1.0,0);
    pad41->Draw();
    pad42->Draw();

  pad11->cd() ; gPad->SetLeftMargin(0.25) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  pad12->cd() ; gPad->SetLeftMargin(0.25) ; xframe_2->GetYaxis()->SetTitleOffset(1.4) ; xframe_2->Draw() ;
  pad21->cd() ; gPad->SetLeftMargin(0.25) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  pad22->cd() ; gPad->SetLeftMargin(0.25) ; xframe_4->GetYaxis()->SetTitleOffset(1.4) ; xframe_4->Draw() ;
  pad41->cd() ; gPad->SetLeftMargin(0.25) ; xframe_5->GetYaxis()->SetTitleOffset(1.4) ; xframe_5->Draw() ;
  pad42->cd() ; gPad->SetLeftMargin(0.25) ; xframe_6->GetYaxis()->SetTitleOffset(1.4) ; xframe_6->Draw() ;


  pad31_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_1->Draw() ;
  pad31_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame1->Draw() ;

  pad32_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_2->Draw() ;
  pad32_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame2->Draw() ;

  pad33_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_3->Draw() ;
  pad33_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame3->Draw() ;

  pad34_m->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_4->Draw() ;
  pad34_p->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); frame4->Draw() ;





////////////////////////////////////////////////////////////////////////////////////////


Float_t Asy, e_Asy, Nt1, Nt2, e_Nt1, e_Nt2;
Asy=Araw.getVal();
e_Asy=Araw.getError();

/*
Nt1=N_t1.getVal();
e_Nt1=N_t1.getError();

Nt2=N_t2.getVal();
e_Nt2=N_t2.getError();
*/

cout<<"Araw = "<<Asy<<" +/- "<<e_Asy<<endl;
cout<<"Araw(truth)= "<<((psig1+psig2)-(nsig1+nsig2))/((psig1+psig2)+(nsig1+nsig2))<<" +/- "<<(1/sqrt(nsig1+psig1+nsig2+psig2))<<endl<<endl;
/*
cout<<"N_t1 = "<<Nt1<<" +/- "<<e_Nt1<<endl;
cout<<"N_t1(truth) = "<<nsig1+psig1<<" +/- "<<sqrt(nsig1+psig1)<<endl<<endl;
cout<<"N_t2 = "<<Nt2<<" +/- "<<e_Nt2<<endl;
cout<<"N_t2(truth) = "<<nsig2+psig2<<" +/- "<<sqrt(nsig2+psig2)<<endl;
*/
}












