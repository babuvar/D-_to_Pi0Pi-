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
  RooRealVar w("w","w",0.0,1.0,"no unit") ;
  RooRealVar CosTheta("CosTheta","CosTheta",-1.0,1.0,"no unit");
  RooDataSet* data_pipi=new RooDataSet("data_pipi","data_pipi",RooArgSet(CosTheta,w));
  RooDataSet* data_kspi=new RooDataSet("data_kspi","data_kspi",RooArgSet(CosTheta,w));
  RooDataSet* data_pipi2=new RooDataSet("data_pipi2","data_pipi2",RooArgSet(CosTheta,w));
  RooDataSet* data_kspi2=new RooDataSet("data_kspi2","data_kspi2",RooArgSet(CosTheta,w));


 

void plots_DstTheta_forPstBins(void)
{
  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
//  TChain* chain=new TChain("h2");



   chain->Add("mer_pipi.root");

int nsig_pipi=0, nsig_kspi=0, nsig_pipi2=0, nsig_kspi2=0;

  Int_t nevt=(int)chain->GetEntries();
//  cout<<"nevt\t"<<nevt <<endl;

  chain->Add("mer_kspi.root");
Int_t nevt2=(int)chain->GetEntries();
 

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df, f_Ksmom, f_Picthet, f_Phimom, f_Cthetdst;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dmass",&f_dzero);
  h1->SetBranchAddress("Pizsmass",&f_Pizsmass);
  h1->SetBranchAddress("Pizmass",&f_Pizmass);
  h1->SetBranchAddress("Dstf",&f_Dstf);
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
  h1->SetBranchAddress("Picthet",&f_Picthet);
  h1->SetBranchAddress("Cthetdst",&f_Cthetdst);


  h1->SetBranchAddress("Ksmom",&f_Ksmom);
  h1->SetBranchAddress("Phimom",&f_Phimom);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


//for(int i=0;i<numbins;i++){cout<<"cut["<<i<<"] = "<<cut[i]<<endl;}




  for(int i=0;i<nevt;i++) 
//  for(int i=0;i<100000;i++)
//  for(int i=0;i<10;i++)
    {


      chain->GetEntry(i);
		  CosTheta.setVal(f_Cthetdst);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
//if(f_Pizsmass > 0.11  && f_Pizsmass < 0.16){
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 && f_Pimom < 5.84){
//if(f_PD > 2.95){
if(f_Df == 1){
//if(f_Dstf == 1){



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

if(f_PD > 2.95){
w.setVal(0.00002018);
data_pipi->add(RooArgSet(CosTheta,w));
nsig_pipi++;}

if(f_PD > 2.5 && f_PD < 2.95){
w.setVal(0.000060376);
data_pipi2->add(RooArgSet(CosTheta,w));
nsig_pipi2++;}


}}//Photon cuts

}}}}}}}//}
    }
 cout<<"nsig_pipi = "<<nsig_pipi<<endl;
 cout<<"nsig_pipi2 = "<<nsig_pipi2<<endl;


  for(int i=nevt;i<nevt2;i++) 
    {


      chain->GetEntry(i);
		  CosTheta.setVal(f_Cthetdst);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.80  && f_dzero < 1.94){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Ksmom > 1.06 ){				//***********
if(f_Pimom > 0.84 ){
//if(f_PD > 2.95){
if(f_Df == 1){

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
w.setVal(0.000019589);
data_kspi->add(RooArgSet(CosTheta,w));
nsig_kspi++;}

if(f_PD > 2.5 && f_PD < 2.95){
 w.setVal(0.000061652);
data_kspi2->add(RooArgSet(CosTheta,w));
nsig_kspi2++;}

}//Photon cuts

}}}}}}//}

    }

 cout<<"nsig_kspi = "<<nsig_kspi<<endl;
 cout<<"nsig_kspi2 = "<<nsig_kspi2<<endl;

 








  RooDataSet* data_pipi_w=new RooDataSet(data_pipi->GetName(),data_pipi->GetTitle(),data_pipi,*data_pipi->get(),0,w.GetName()) ;
  RooDataSet* data_kspi_w=new RooDataSet(data_kspi->GetName(),data_kspi->GetTitle(),data_kspi,*data_kspi->get(),0,w.GetName()) ;
  RooDataSet* data_pipi2_w=new RooDataSet(data_pipi2->GetName(),data_pipi2->GetTitle(),data_pipi2,*data_pipi2->get(),0,w.GetName()) ;
  RooDataSet* data_kspi2_w=new RooDataSet(data_kspi2->GetName(),data_kspi2->GetTitle(),data_kspi2,*data_kspi2->get(),0,w.GetName()) ;

 

  TCanvas* can = new TCanvas("c","c") ;
  TCanvas* can2 = new TCanvas("c2","c2") ;
  TCanvas* can3 = new TCanvas("c3","c3") ;
  TCanvas* can4 = new TCanvas("c4","c4") ;
  TCanvas* can5 = new TCanvas("c5","c5") ;



  //DeltaM PLOTING
  RooPlot *xframe =CosTheta.frame(Bins(50),Title(" "));
  data_pipi_w->plotOn(xframe, MarkerColor(kRed));
  data_kspi_w->plotOn(xframe, MarkerColor(kGreen));
  data_pipi2_w->plotOn(xframe, MarkerColor(7));
  data_kspi2_w->plotOn(xframe, MarkerColor(kBlue));
 

  RooPlot *xframe2 =CosTheta.frame(Bins(50),Title(" "));
  data_pipi2_w->plotOn(xframe2, MarkerColor(7));
  data_kspi2_w->plotOn(xframe2, MarkerColor(kBlue));

  RooPlot *xframe3 =CosTheta.frame(Bins(50),Title(" "));
   data_pipi_w->plotOn(xframe3, MarkerColor(kRed));
  data_kspi_w->plotOn(xframe3, MarkerColor(kGreen));


  RooPlot *xframe4 =CosTheta.frame(Bins(50),Title(" "));
   data_pipi_w->plotOn(xframe4, MarkerColor(kRed));
  data_pipi2_w->plotOn(xframe4, MarkerColor(7));


  RooPlot *xframe5 =CosTheta.frame(Bins(50),Title(" "));
  data_kspi_w->plotOn(xframe5, MarkerColor(kGreen));
   data_kspi2_w->plotOn(xframe5, MarkerColor(kBlue));



can ->cd();  xframe->Draw();
can2 ->cd();  xframe2->Draw();
can3 ->cd();   xframe3->Draw();
can4 ->cd();  xframe4->Draw();
can5 ->cd();  xframe5->Draw();




}















