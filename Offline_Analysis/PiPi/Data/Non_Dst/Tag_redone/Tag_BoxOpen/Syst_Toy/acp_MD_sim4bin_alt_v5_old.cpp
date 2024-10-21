//Tagged 4-Sim fitting for tagged pipi with new model
// + weights for signal and background

 
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
#include <fstream>

using namespace RooFit;
void DoFit(int ii);
double result[4][10000];
double low_lim=1.68;

double mudif1=0.00073;//0.00028
double mudif2= 0.0188;// 
double mudif3= 0.022;// 
double sigrat1= 0.975;//
double sigrat2= 3.79;//
double sigrat3= 3.65;//
double f1sig= 0.9799;//
double f2sig= 0.9746;//
double al1sig= 0.805;//
double al2sig= 0.781;//
double n1sig= 4.92;//
double n2sig= 3.77;//
double al1bkg= -1.9014;//
double al2bkg= -1.9129;//
double n1bkg= 1.604;//
double n2bkg= 1.328;//
double mudifbkg= -0.00135;//
double sigratbkg= 0.999;//


  RooRealVar dzero("dzero","M_{D}  ",low_lim,2.06,"GeV");
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dzero));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dzero));
  RooDataSet* data_p2=new RooDataSet("data_p2","data_p2",RooArgSet(dzero));
  RooDataSet* data_n2=new RooDataSet("data_n2","data_n2",RooArgSet(dzero));

  TH1D* Araw_res=new TH1D("Araw_res", "Araw_res", 50, 0.3, 0.8);

float nsig1=0,nsig2=0,psig1=0,psig2=0;

void acp_MD_sim4bin_alt_v5(char * filename)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");

   chain->Add("../../Root_Files/pipi_dcont.root"); 
  chain->Add("../../Root_Files/pipi_d4s.root"); 
  chain->Add("../../Root_Files/pipi_d5s.root"); 


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
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >low_lim  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){


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


if(f_Dcharge == 1){data_p->add(RooArgSet(dzero));   }
if(f_Dcharge == -1){data_n->add(RooArgSet(dzero));  }


}


//2nd P*D* bin
if(f_PD > 2.5 && f_PD <= 2.95){

if(f_Dcharge == 1){data_p2->add(RooArgSet(dzero));  }
if(f_Dcharge == -1){data_n2->add(RooArgSet(dzero)); }


}



//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}}//Photon cuts

}}}}}}//}//}
    }
TRandom *seed = new TRandom(0);
seed->SetSeed(0);

TRandom *ran = new TRandom(seed->Rndm());
TRandom *ran1 = new TRandom(seed->Rndm());
TRandom *ran2 = new TRandom(seed->Rndm());
TRandom *ran3 = new TRandom(seed->Rndm());
TRandom *ran4 = new TRandom(seed->Rndm());
TRandom *ran5 = new TRandom(seed->Rndm());
TRandom *ran6 = new TRandom(seed->Rndm());
TRandom *ran7 = new TRandom(seed->Rndm());
TRandom *ran8 = new TRandom(seed->Rndm());
TRandom *ran9 = new TRandom(seed->Rndm());
TRandom *ran10 = new TRandom(seed->Rndm());
TRandom *ran11 = new TRandom(seed->Rndm());
TRandom *ran12 = new TRandom(seed->Rndm());
TRandom *ran13 = new TRandom(seed->Rndm());
TRandom *ran14 = new TRandom(seed->Rndm());
TRandom *ran15 = new TRandom(seed->Rndm());
TRandom *ran16 = new TRandom(seed->Rndm());
TRandom *ran17 = new TRandom(seed->Rndm());





//double num_iters=1000;
double num_iters=5;

for(int ii=0;ii<num_iters;ii++){

mudif1=ran->Gaus(0.00073,0.00028);
mudif2= ran1->Gaus(0.0188,0.0086);
mudif3=   ran2->Gaus(0.022,0.011);
sigrat1=  ran3->Gaus(0.975,0.015);
sigrat2=  ran4->Gaus(3.79,0.23);
sigrat3=  ran5->Gaus(3.65,0.3);
f1sig=  ran6->Gaus(0.9799,0.003);
f2sig=  ran7->Gaus(0.9746,0.0051);
al1sig=  ran8->Gaus(0.805,0.014);
al2sig=  ran9->Gaus(0.781,0.026);
n1sig=  ran10->Gaus(4.92,0.20);
n2sig=  ran11->Gaus(3.77,0.23);
al1bkg=  ran12->Gaus(-1.9014,0.031);
al2bkg=  ran13->Gaus(-1.9129,0.038);
n1bkg=  ran14->Gaus(1.604,0.064);
n2bkg=  ran15->Gaus(1.328,0.067);
mudifbkg=  ran16->Gaus(-0.00135,0.0022);
sigratbkg=  ran17->Gaus(0.999,0.032);

DoFit(ii);
//Araw_res->Fill(result[0][ii]);
}


/*

//    TF1 *func = new TF1("myGauss","([0]*exp(-0.5*((x-[2])/[3])**2)) +( [1]*exp(-0.5*((x-[2])/[4])**2))",3.4,4.2);
    TF1 *func = new TF1("myGauss","([0]*exp(-0.5*((x-[2])/[3])**2)) +( [1]*exp(-0.5*((x-[2])/[4])**2))",0.25,1.75);
    // parameter names
    func->SetParNames("Factor1","Factor2","Mean","Sigma1","Sigma2");
    //Set Parameters:
    func->SetParameter(0,90000);
    func->SetParLimits(0,0,100000);
    func->SetParameter(1,40000);
    func->SetParLimits(1,0,50000);
    func->SetParameter(2,0.52);//mean
    func->SetParLimits(2,0.48,0.56);//mean
    func->SetParameter(3,0.1);//Sigma1
    func->SetParLimits(3,0.01,0.5);//Sigma1
    func->SetParameter(4,0.1);
    func->SetParLimits(4,0.01,2.0);


TCanvas* can = new TCanvas("can","can") ;
Araw_res->SetMarkerStyle(20);
//Araw_res->SetMarkerColor(kBlue);

Araw_res->Fit(func,"R");// use option "R" to restrict to a certain region for fitting
//Araw_res->Fit("gaus");
can->cd();
Araw_res->Draw("EP");
*/
ofstream fout;
//fout.open("out1.txt");
fout.open(filename);
//Print values
for(int ii=0;ii<num_iters;ii++){
fout<<result[0][ii]<<endl;
}
fout.close();

exit();

}


void DoFit(int ii)
{
  //DEFINE PDF
//Common
  //____________________________________________________
  RooRealVar Araw("A_{Raw}","Araw",0,-1,1);
  RooRealVar N_t1("N_{Sig1}","N_t1",100,0,150000);
  RooFormulaVar N_n("N_n","(0.5)*(1-A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));
  RooFormulaVar N_p("N_p","(0.5)*(1+A_{Raw})*N_{Sig1}",RooArgList(Araw,N_t1));


  RooRealVar N_t2("N_{Sig2}","N_t2",100,0,60000);
  RooFormulaVar N_n2("N_n2","(0.5)*(1-A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));
  RooFormulaVar N_p2("N_p2","(0.5)*(1+A_{Raw})*N_{Sig2}",RooArgList(Araw,N_t2));



  RooRealVar Abkg1("A_{Bkg1}","ABkg1",0,-1,1);
  RooRealVar N_tb1("N_{Bkg1}","N_tb1",100,0,10000000);
  RooFormulaVar N_nb("N_nb","(0.5)*(1-A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));
  RooFormulaVar N_pb("N_pb","(0.5)*(1+A_{Bkg1})*N_{Bkg1}",RooArgList(Abkg1,N_tb1));

  RooRealVar Abkg2("A_{Bkg2}","ABkg2",0,-1,1);
  RooRealVar N_tb2("N_{Bkg2}","N_tb2",100,0,10000000);
  RooFormulaVar N_nb2("N_nb2","(0.5)*(1-A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));
  RooFormulaVar N_pb2("N_pb2","(0.5)*(1+A_{Bkg2})*N_{Bkg2}",RooArgList(Abkg2,N_tb2));



  //_____________________________________________________


//Common Pars
 
  RooRealVar m_dif1("m_{shift(CB2,CB1)}","m_dif1",mudif1);
  RooRealVar m_dif2("m_{shift(Gaus1,CB1)}","m_dif2",mudif2);
  RooRealVar m_dif3("m_{shift(Gaus2,CB1)}","m_dif3",mudif3);
  RooRealVar s_rat1("s_{ratio(CB2,CB1)}","s_rat1",sigrat1);
  RooRealVar s_rat2("s_{ratio(Gaus1,CB1)}","s_rat2",sigrat2);
  RooRealVar s_rat3("s_{ratio(Gaus2,CB1)}","s_rat3",sigrat3);
  RooRealVar m_difBkg("m_{shiftBkg}","m_difBkg",mudifbkg);
  RooRealVar s_ratBkg("s_{ratioBkg}","s_ratBkg",sigratbkg);

  RooRealVar a("#alpha","a",al1sig);
  RooRealVar n("n","n",n1sig);
  RooRealVar f_sig("f_sig","f_sig",f1sig);
  RooRealVar a2("#alpha2","a2",al2sig);
  RooRealVar n2("n2","n2",n2sig);
  RooRealVar f_sig2("f_sig2","f_sig2",f2sig);

  RooRealVar a_bkg("#alpha_bkg","a_bkg",al1bkg);
  RooRealVar n_bkg("n_bkg","n_bkg",n1bkg);
  RooRealVar a_bkg2("#alpha_bkg2","a_bkg2",al2bkg);
  RooRealVar n_bkg2("n_bkg2","n_bkg2",n2bkg);





//Signal
  RooRealVar m_fudge("#mu_{fudge}","m_fudge",-0.000635,-0.005,0.005);
  RooRealVar m_fudgeBkg("#mu_{fudgeBkg}","m_fudgeBkg",-0.0,-0.04,0.04);
  RooRealVar s_fudge("#sigma_{fudge}","s_fudge",1.000,0.5,2.5);
  RooRealVar s_fudgeBkg("#sigma_{fudgeBkg}","s_fudgeBkg",1.000,0.5,2.5);


//Bin-1
  RooFormulaVar m("m_sig","1.86805+#mu_{fudge}",RooArgList(m_fudge));
  RooFormulaVar s("s","0.01530*#sigma_{fudge}",RooArgList(s_fudge));


  RooFormulaVar m_g("m_g","1.86805+#mu_{fudge}+(m_{shift(Gaus1,CB1)}*#sigma_{fudge})",RooArgList(m_fudge,s_fudge,m_dif2));
  RooFormulaVar s_g("s_g","0.01530*#sigma_{fudge}*s_{ratio(Gaus1,CB1)}",RooArgList(s_fudge,s_rat2));


//Bin-2
  RooFormulaVar m2("m_sig2","1.86805+#mu_{fudge}+m_{shift(CB2,CB1)}",RooArgList(m_fudge,m_dif1));
  RooFormulaVar s2("s2","0.01530*#sigma_{fudge}*s_{ratio(CB2,CB1)}",RooArgList(s_fudge,s_rat1));


  RooFormulaVar m_g2("m_g2","1.86805+#mu_{fudge}+(m_{shift(Gaus2,CB1)}*#sigma_{fudge})",RooArgList(m_fudge,s_fudge,m_dif3));
  RooFormulaVar s_g2("s_g2","0.0558*#sigma_{fudge}*s_{ratio(Gaus2,CB1)}",RooArgList(s_fudge,s_rat3));


//Combinatorial
//Bin-1
  RooRealVar p1("p1","p1",-0.4,-0.8,0.01);
  RooRealVar p1_2("p1_2","p1_2",-0.4,-0.8,0.01);

//Bin-2
  RooRealVar p2("p2","p2",-0.3,-0.8,0.01);
  RooRealVar p2_2("p2_2","p2_2",-0.3,-0.8,0.01);



//Peaking Background
//Bin-1
  RooFormulaVar m_bkg("#mu_{bkg}","1.6802+#mu_{fudgeBkg}",RooArgList(m_fudgeBkg));
  RooFormulaVar s_bkg("#sigma_{bkg}","0.03656*#sigma_{fudgeBkg}",RooArgList(s_fudgeBkg));



//Bin-2
  RooFormulaVar m_bkg2("#mu_{bkg2}","1.6802+#mu_{fudgeBkg}+m_{shiftBkg}",RooArgList(m_fudgeBkg,m_difBkg));
  RooFormulaVar s_bkg2("#sigma_{bkg2}","0.03656*#sigma_{fudgeBkg}*s_{ratioBkg}",RooArgList(s_fudgeBkg,s_ratBkg));



//Background fraction
//Bin-1
RooRealVar f_bkg("f_bkg","f_bkg",0.5,0.0,1.0);


//Bin-2
RooRealVar f_bkg2("f_bkg2","f_bkg2",0.5,0.0,1.0);



//_____________________________________________________________________________________
//			BIN - 1
//_____________________________________________________________________________________

//Signal - P
  RooCBShape Sig1_p("Sig1_p", "Cystal Ball Function_p", dzero, m, s, a, n); 
  RooGaussian Sig2_p("Sig2_p", "Sig2_p",dzero,m_g,s_g);
  RooAddPdf Sig_p("Sig_p","Sig_p",RooArgList(Sig1_p,Sig2_p),RooArgList(f_sig));


//Background - P
  RooChebychev Bkg1_p("Bkg1_p", "Bkg1_p", dzero, RooArgList(p1));
  RooCBShape Bkg2_p("Bkg2_p", "Bkg2_p", dzero, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_p("Bkg_p","Bkg_p",RooArgList(Bkg1_p,Bkg2_p),RooArgList(f_bkg));

//Full Model - P
  RooAddPdf Model_p("Model_p","Model_p",RooArgList(Sig_p,Bkg_p),RooArgList(N_p,N_pb));

//Signal - N
  RooCBShape Sig1_n("Sig1_n", "Cystal Ball Function_n", dzero, m, s, a, n); 
  RooGaussian Sig2_n("Sig2_n", "Sig2_n",dzero,m_g,s_g);
  RooAddPdf Sig_n("Sig_n","Sig_n",RooArgList(Sig1_n,Sig2_n),RooArgList(f_sig));


//Background - N  
//Data 
  RooChebychev Bkg1_n("Bkg1_n", "Bkg1_n", dzero, RooArgList(p1));
  RooCBShape Bkg2_n("Bkg2_n", "Bkg2_n", dzero, m_bkg, s_bkg, a_bkg, n_bkg); 
  RooAddPdf Bkg_n("Bkg_n","Bkg_n",RooArgList(Bkg1_n,Bkg2_n),RooArgList(f_bkg));

//Full Model - N
  RooAddPdf Model_n("Model_n","Model_n",RooArgList(Sig_n,Bkg_n),RooArgList(N_n,N_nb));


//_____________________________________________________________________________________
//			BIN - 2
//_____________________________________________________________________________________


//Signal - P
  RooCBShape Sig1_p2("Sig1_p2", "Cystal Ball Function_p2", dzero, m2, s2, a2, n2); 
  RooGaussian Sig2_p2("Sig2_p2", "Sig2_p2",dzero,m_g2,s_g2);
  RooAddPdf Sig_p2("Sig_p2","Sig_p2",RooArgList(Sig1_p2,Sig2_p2),RooArgList(f_sig2));


//Background - P
  RooChebychev Bkg1_p2("Bkg1_p2", "Bkg1_p2", dzero, RooArgList(p2));
  RooCBShape Bkg2_p2("Bkg2_p2", "Bkg2_p2", dzero, m_bkg2, s_bkg2, a_bkg2, n_bkg2); 
  RooAddPdf Bkg_p2("Bkg_p2","Bkg_p2",RooArgList(Bkg1_p2,Bkg2_p2),RooArgList(f_bkg2));

//Full Model - P
  RooAddPdf Model_p2("Model_p2","Model_p2",RooArgList(Sig_p2,Bkg_p2),RooArgList(N_p2,N_pb2));

//Signal - N
  RooCBShape Sig1_n2("Sig1_n2", "Cystal Ball Function_n2", dzero, m2, s2, a2, n2); 
  RooGaussian Sig2_n2("Sig2_n2", "Sig2_n2",dzero,m_g2,s_g2);
  RooAddPdf Sig_n2("Sig_n2","Sig_n2",RooArgList(Sig1_n2,Sig2_n2),RooArgList(f_sig2));


//Background - N  
//Data 
  RooChebychev Bkg1_n2("Bkg1_n2", "Bkg1_n2", dzero, RooArgList(p2));
  RooCBShape Bkg2_n2("Bkg2_n2", "Bkg2_n2", dzero, m_bkg2, s_bkg2, a_bkg2, n_bkg2); 
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




Float_t Asy, e_Asy, Nt1, Nt2, e_Nt1, e_Nt2;
//Asy=Araw.getVal();
//e_Asy=Araw.getError();

//Nt1=N_t1.getVal();
//e_Nt1=N_t1.getError();

//Nt2=N_t2.getVal();
//e_Nt2=N_t2.getError();


//cout<<"Araw = "<<Asy<<" +/- "<<e_Asy<<endl;

//cout<<"N_t1 = "<<Nt1<<" +/- "<<e_Nt1<<endl;

//cout<<"N_t2 = "<<Nt2<<" +/- "<<e_Nt2<<endl;

 result[0][ii]=Araw.getVal()*100;
 result[1][ii]=Araw.getError()*100;
// result[2][ii]=N_t.getVal();
// result[3][ii]=N_t.getError();

}












