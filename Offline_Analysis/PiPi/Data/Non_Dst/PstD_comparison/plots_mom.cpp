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
  RooRealVar DMom("DMom","DMom",2.5,5.0,"GeV");
  RooDataSet* data_tag=new RooDataSet("data_tag","data_tag",RooArgSet(DMom,w));
  RooDataSet* data_untag=new RooDataSet("data_untag","data_untag",RooArgSet(DMom,w));


void plots_mom(void)
{
int nsig_kspi=0, nsig_tag=0, nsig_untag=0;

  //LOAD DATA FILE
  TChain* chain=new TChain("h1");
  TChain* chain2=new TChain("h2");


  chain->Add("pipi_4smc_0F.root");
  chain->Add("pipi_SMC_2M_4s.root");

  chain2->Add("pipi_4smc_0F.root");

  Int_t nevt=(int)chain->GetEntries();
  Int_t nevt2=(int)chain2->GetEntries();


//  cout<<"nevt\t"<<nevt <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Dcharge, f_Pizmass, f_Df, f_Ksmom, f_Phimom,f_Pstdst;


  h1->SetBranchAddress("Deltam",&Deltam);
  h1->SetBranchAddress("Pstdst",&f_Pstdst);
  h1->SetBranchAddress("Pstd",&f_PD);
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
 


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  float bin;int bin1; 

int num_pass=0;


  for(int i=0;i<nevt;i++) 
    {


      chain->GetEntry(i);
		  DMom.setVal(f_PD);

if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84){
if(f_Pstdst > 2.95){
if(f_Df == 1){




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

w.setVal(0.00001672884);
data_tag->add(RooArgSet(DMom,w));


nsig_tag++;


}}//Photon cuts

}}}}}}}}//}
    }
 cout<<"nsig_tag = "<<nsig_tag<<endl;



  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Df",&f_Df);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Pstd",&f_PD);


  for(int i=0;i<nevt2;i++) 
    {


      chain2->GetEntry(i);
		  DMom.setVal(f_PD);


if(f_dzero >1.68  && f_dzero < 2.06){       //~3sigma range to estimate F.O.M.
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84){
if(f_PD > 2.65){
if(f_Df == 1){




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


w.setVal(0.00000913025);
data_untag->add(RooArgSet(DMom,w));


nsig_untag++;


}//Photon cuts

}}}}}}//}
    }
 cout<<"nsig_untag = "<<nsig_untag<<endl;



  RooDataSet* data_tag_w=new RooDataSet(data_tag->GetName(),data_tag->GetTitle(),data_tag,*data_tag->get(),0,w.GetName()) ;
  RooDataSet* data_untag_w=new RooDataSet(data_untag->GetName(),data_untag->GetTitle(),data_untag,*data_untag->get(),0,w.GetName()) ;



  TCanvas* can = new TCanvas("c","c") ;
  //DeltaM PLOTING
  RooPlot *xframe =DMom.frame(Bins(50),Title("D^{+} #rightarrow #pi^{0} #pi^{+}"));
  data_tag_w->plotOn(xframe, MarkerColor(kRed));
  data_untag_w->plotOn(xframe, MarkerColor(kGreen));
  xframe->Draw();



}
















