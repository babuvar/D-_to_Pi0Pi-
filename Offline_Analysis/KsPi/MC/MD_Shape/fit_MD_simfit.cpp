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
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(dzero));
  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero));
int numsign=0;

void fit_MD_simfit()
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");


  chain->Add("KsPi_GMC0_4s.root");
  chain->Add("KsPi_GMC0_5s.root");
  chain->Add("KsPi_GMC1_4s.root");
  chain->Add("KsPi_GMC1_5s.root");


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
  h1->SetBranchAddress("Df",&f_Df);

int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;



  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<10000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);
		  dzero.setVal(f_dzero);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){

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
data->add(RooArgSet(dzero));
//if(f_Df==1){data->add(RooArgSet(dzero));}
}


if(f_PD > 2.5 && f_PD < 2.95){
data2->add(RooArgSet(dzero));
//if(f_Df==1){data2->add(RooArgSet(dzero));}
}
numsign++;

//cout<<"fabs(f_Pizmass - 0.135) = "<<fabs(f_Pizmass - 0.135)<<endl;
//cout<<"fabs(f_Pizsmass - 0.135) = "<<fabs(f_Pizsmass - 0.135)<<endl;
//cout<<"---------------------------"<<endl;


}//Photon cuts

}}}}}
    }
 cout<<"numsign = "<<numsign<<endl;


  //DEFINE PDF
//Common
  RooRealVar N_sig("N_{Sig}","N_sig",100,0,1000000);
  RooRealVar N_bkg("N_{Bkg}","N_bkg",100,0,1000000);
  RooRealVar N_sig2("N_{Sig2}","N_sig2",100,0,1000000);
  RooRealVar N_bkg2("N_{Bkg2}","N_bkg2",100,0,1000000);

  RooRealVar mean("mean","mean",1.87,1.86,1.88);
  RooRealVar sig("sig","sig",0.009,0.008,0.016);
  RooRealVar sig1("sig1","sig1",0.005,0.002,0.006);
  RooRealVar sig2("sig2","sig2",0.005,0.002,0.006);
  RooRealVar f_Sig("f_Sig","f_Sig",0.3,0.15,0.5);


//Signal
  RooGaussian sig_p1("sig_p1", "signal Gaussian 1",dzero,mean,sig);
  RooBifurGauss sig_p2("sig_p2", "signal Gaussian 2",dzero,mean,sig1,sig2);
  RooAddPdf Sig("Sig_p","Sig",RooArgList(sig_p1,sig_p2),f_Sig);

//Signal-2
  RooRealVar s_fact("#sigma_{factor}","s_fact",1.0005,0.0,2.0);
  RooFormulaVar sig_2("#sigma2","sig * #sigma_{factor}",RooArgList(sig,s_fact));
  RooFormulaVar sig12("#sigma12","sig1 * #sigma_{factor}",RooArgList(sig1,s_fact));
  RooFormulaVar sig22("#sigma22","sig2 * #sigma_{factor}",RooArgList(sig2,s_fact));


  RooGaussian sig_p12("sig_p12", "signal Gaussian 12",dzero,mean,sig_2);
  RooBifurGauss sig_p22("sig_p22", "signal Gaussian 22",dzero,mean,sig12,sig22);
  RooAddPdf Sig2("Sig2","Sig2",RooArgList(sig_p12,sig_p22),f_Sig);

//Background
   RooRealVar p("p","p",0.5,-1,1);   
  RooChebychev Bkg("Bkg", "Bkg", dzero, RooArgList(p));

//Background-2
   RooRealVar p2("p2","p2",0.5,-1,1);   
  RooChebychev Bkg2("Bkg2", "Bkg2", dzero, RooArgList(p2));

//Full Model
  RooAddPdf Model("Model","Model",RooArgList(Sig,Bkg),RooArgList(N_sig,N_bkg));
  RooAddPdf Model2("Model2","Model2",RooArgList(Sig2,Bkg2),RooArgList(N_sig2,N_bkg2));

  //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("bin1");
  sample.defineType("bin2");


  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dzero),Index(sample),Import("bin1",*data),Import("bin2",*data2));
  


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(Model,"bin1");
  simPdf.addPdf(Model2,"bin2");
//  simPdf.addPdf(Sig,"bin1");
//  simPdf.addPdf(Sig2,"bin2");


//Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());


  TCanvas* can = new TCanvas("c","c",700,800) ;
           can->Divide(2,1);

           can->cd(1);
  RooPlot* dzero_frame = dzero.frame(Bins(30),Title("Mass of D candidate"));
  dzero_frame->SetTitle("Mass of D candidate");
  data->plotOn(dzero_frame);
  Model.plotOn(dzero_frame, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
  Model.plotOn(dzero_frame,Components("Bkg"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
  Model.paramOn(dzero_frame);
  dzero_frame->Draw();

           can->cd(2);
  RooPlot* dzero_frame2 = dzero.frame(Bins(30),Title("Mass of D candidate"));
  dzero_frame2->SetTitle("Mass of D candidate");
  data2->plotOn(dzero_frame2);
  Model2.plotOn(dzero_frame2, LineColor(kBlue), LineStyle(kSolid),LineWidth(2));
  Model2.plotOn(dzero_frame2,Components("Bkg2"), LineColor(kBlue), LineStyle(kDashed),LineWidth(2));
  Model2.paramOn(dzero_frame2);
  dzero_frame2->Draw();


}

