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
float x[100],yb[100],ex[100],ys[100],eys[100],eyb[100];
float  yb_c[100] ,ys_c[100],eys_c[100]=0,eyb_c[100]=0;

void ComPlot_DatMC(void)
{


float numbins=50;
float binwidth=0.050/numbins;

for(int i=0;  i<=99; i++){x[i]=0.110+((i+0.5)*binwidth); ex[i]=0;}

 
//Reading
//Text file with results
ifstream fin;
fin.open("Results_Signal.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>ys[i]>>eys[i];
}
fin.close();
fin.open("Results_Background.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>yb[i]>>eyb[i];
}
fin.close();

//Text file with results (counted)
//fin.open("Results_Signal_count.txt");
fin.open("Results_Signal_data.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>ys_c[i]>>eys_c[i];
}
fin.close();
//fin.open("Results_Background_count.txt");
fin.open("Results_Background_data.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>yb_c[i]>>eyb_c[i];
}
fin.close();




// PLOTS OF ESTIMATED SIGNAL AND BACKGROUND
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;
      	TCanvas* cnv2 = new TCanvas("cnv2","cnv2") ;


   TGraphErrors *gr1 = new TGraphErrors(numbins,x,ys,ex,eys);
   TGraphErrors *gr2 = new TGraphErrors(numbins,x,yb,ex,eyb);
   TGraphErrors *gr3 = new TGraphErrors(numbins,x,ys_c,ex,eys_c);
   TGraphErrors *gr4 = new TGraphErrors(numbins,x,yb_c,ex,eyb_c);

   	gr1->SetFillColor(38);   	gr3->SetFillColor(38);
   	gr2->SetFillColor(kRed-3);   gr2->SetMinimum(0.0);     	gr4->SetFillColor(kRed-3);   gr4->SetMinimum(0.0);


//	gr1->Fit("gaus");
//   	cnv1->cd(); gr1->Draw("AB");
//   	cnv2->cd(); gr2->Draw("AB");

//      	TCanvas* cnv3 = new TCanvas("cnv3","cnv3") ;   	cnv3->cd();
   gr1->SetMarkerColor(38);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);
   gr2->SetMarkerColor(38);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);

   gr3->SetMarkerColor(kRed-3);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr4->SetMarkerColor(kRed-3);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);


cnv1->cd();
     TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Mass of slow #pi^{0}");
     mg->Add(gr1,"p");
     mg->Add(gr3,"p");
     mg->Draw("A");
  leg = new TLegend(0.6,0.7,0.89,0.89);
//  leg->AddEntry(gr1,"Signal (Binned Fits)","lep");
//  leg->AddEntry(gr3,"Signal Truth Matched","lep");
  leg->AddEntry(gr1,"Signal MC (Binned Fits)","lep");
  leg->AddEntry(gr3,"Signal data (Binned Fits)","lep");

  leg->Draw();

cnv2->cd();
     TMultiGraph *mg2 = new TMultiGraph();
     mg2->SetTitle("Mass of slow #pi^{0}");
     mg2->Add(gr2,"p");
     mg2->Add(gr4,"p");
     mg2->Draw("A");
  leg2 = new TLegend(0.6,0.7,0.89,0.89);
//  leg2->AddEntry(gr2,"Background (Binned Fits)","lep");
//  leg2->AddEntry(gr4,"Background Truth Matched","lep");
  leg2->AddEntry(gr2,"Background MC (Binned Fits)","lep");
  leg2->AddEntry(gr4,"Background data (Binned Fits)","lep");
  leg2->Draw();



   cnv1->Update();//  cnv1->SaveAs("Signal.png");
   cnv2->Update();//  cnv2->SaveAs("Background.png");
//   cnv3->Update();  cnv3->SaveAs("Overlay.png");





}






