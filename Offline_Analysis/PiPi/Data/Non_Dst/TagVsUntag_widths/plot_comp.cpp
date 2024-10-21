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
  RooRealVar dzero("dzero","M_{D}   ",1.68,2.0,"GeV");
  RooRealVar pimom("pimom","P_{#pi^{#pm}}  ",0.0,5.0,"GeV");
  RooRealVar pizmom("pizmom","P_{#pi^{0}}  ",0.0,5.0,"GeV");


  RooRealVar w("w","w",0.0,0.1,"no unit") ;

  RooDataSet* data1=new RooDataSet("data1","data1",RooArgSet(dzero,pimom,pizmom,w));
  RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(dzero,pimom,pizmom,w));

int nsig1=0, nsig2=0;

void plot_comp(void)
{
  //LOAD DATA FILE
  TChain* chain1=new TChain("h1");
  TChain* chain2=new TChain("h2");



  chain1->Add("pipi_4smc_0F.root");
  chain1->Add("pipi_4smc_1F.root");

  chain2->Add("pipi_5s_gmc0.root");
  chain2->Add("pipi_5s_gmc1.root");



  Int_t nevt1=(int)chain1->GetEntries();
  cout<<"nevt1\t"<<nevt1 <<endl;

  Int_t nevt2=(int)chain2->GetEntries();
  cout<<"nevt2\t"<<nevt2 <<endl;

 Float_t Deltam, f_PD, f_dzero, f_Pizsmass, f_Dstf, f_Egam1s, f_Egam2s, f_Egamma1, f_Egamma2, f_Categ, f_Pizmom, f_Pimom, f_Gam1thet, f_Gam2thet, f_Pizmass, f_Df; 
    Float_t   f_Dcharge;



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

int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;

  for(int i=0;i<nevt1;i++) 
    {


      chain1->GetEntry(i);
		  dzero.setVal(f_dzero);
		  pimom.setVal(f_Pimom);
		  pizmom.setVal(f_Pizmom);


if(Deltam > 0.139 && Deltam < 0.142){       //optimized for M_D fitting
if(f_dzero >1.68  && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_Pizsmass > 0.125  && f_Pizsmass < 0.143){
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){


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

w.setVal(0.000072459);
if(f_Df==1 || f_Df == 10){
data1->add(RooArgSet(dzero,pimom,pizmom,w));
nsig1++;
}



}}//Photon cuts

}}}}}}}
    }
cout<<"nsig1 = "<<nsig1<<endl;



  h2->SetBranchAddress("Pstd",&f_PD);
  h2->SetBranchAddress("Dmass",&f_dzero);
  h2->SetBranchAddress("Pizmass",&f_Pizmass);
  h2->SetBranchAddress("Df",&f_Df);
  h2->SetBranchAddress("Egamma1",&f_Egamma1);
  h2->SetBranchAddress("Egamma2",&f_Egamma2);
  h2->SetBranchAddress("Pizmom",&f_Pizmom);
  h2->SetBranchAddress("Pimom",&f_Pimom);
  h2->SetBranchAddress("Gam1hthe",&f_Gam1thet);
  h2->SetBranchAddress("Gam2hthe",&f_Gam2thet);
  h2->SetBranchAddress("Dcharge",&f_Dcharge);


int photon1cutflag=0, photon2cutflag=0, sphoton1cutflag=0, sphoton2cutflag=0;




  for(int i=0;i<nevt2;i++) 
    {


      chain2->GetEntry(i);
		  dzero.setVal(f_dzero);
 		  pimom.setVal(f_Pimom);
		  pizmom.setVal(f_Pizmom);


if(f_dzero >1.68 && f_dzero < 2.0){       //~3sigma range to estimate F.O.M.
if(f_Pizmass > 0.119  && f_Pizmass < 0.151){
if(f_Pizmom > 1.06 ){
if(f_Pimom > 0.84 ){
if(f_PD > 2.5){


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

w.setVal(0.000042653);

//Photon cuts
if(photon1cutflag == 1 && photon2cutflag == 1){
if(f_Df==1 || f_Df==10){
data2->add(RooArgSet(dzero,pimom,pizmom,w));
nsig2++;
}
}

}}}}}
    }

 cout<<"nsig2 = "<<nsig2<<endl;


  RooDataSet* data1_w=new RooDataSet(data1->GetName(),data1->GetTitle(),data1,*data1->get(),0,w.GetName()) ;
  RooDataSet* data2_w=new RooDataSet(data2->GetName(),data2->GetTitle(),data2,*data2->get(),0,w.GetName()) ;


  TCanvas* can = new TCanvas("c","c",700,700) ;
can->Divide(2,2);

  //DeltaM PLOTING
  RooPlot *xframe =dzero.frame(Bins(38));
  data1_w->plotOn(xframe,RooFit::Name("data1"), MarkerColor(kGreen));
  data2_w->plotOn(xframe,RooFit::Name("data2"), MarkerColor(kBlue));
  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;

   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(xframe->findObject("data1"),"Tagged D #rightarrow #pi #pi sample","p");
   legend->AddEntry(xframe->findObject("data2"),"Untagged D #rightarrow #pi #pi sample","p");
   legend->Draw();



  //PiMom PLOTING
  RooPlot *xframe2 =pimom.frame(Bins(50));
  data1_w->plotOn(xframe2,RooFit::Name("data1_2"), MarkerColor(kGreen));
  data2_w->plotOn(xframe2,RooFit::Name("data2_2"), MarkerColor(kBlue));
  can->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;

   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(xframe2->findObject("data1_2"),"Tagged D #rightarrow #pi #pi sample","p");
   legend->AddEntry(xframe2->findObject("data2_2"),"Untagged D #rightarrow #pi #pi sample","p");
   legend->Draw();


  //PizMom PLOTING
  RooPlot *xframe3 =pizmom.frame(Bins(50));
  data1_w->plotOn(xframe3,RooFit::Name("data1_3"), MarkerColor(kGreen));
  data2_w->plotOn(xframe3,RooFit::Name("data2_3"), MarkerColor(kBlue));
  can->cd(4) ; gPad->SetLeftMargin(0.15) ; xframe3->GetYaxis()->SetTitleOffset(1.4) ; xframe3->Draw() ;

   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(xframe3->findObject("data1_3"),"Tagged D #rightarrow #pi #pi sample","p");
   legend->AddEntry(xframe3->findObject("data2_3"),"Untagged D #rightarrow #pi #pi sample","p");
   legend->Draw();


}










