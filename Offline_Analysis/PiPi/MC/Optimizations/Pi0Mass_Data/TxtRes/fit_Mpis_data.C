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

void fit_Mpis_data()
{


float numbins=50;
float binwidth=0.050/numbins;

for(int i=0;  i<=99; i++){x[i]=0.110+((i+0.5)*binwidth); ex[i]=0;}

 

//Text file with results (counted)
ifstream fin;
//fin.open("Results_Signal_count.txt");
fin.open("Results_Signal_data.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>ys_c[i]>>eys_c[i];
}
fin.close();





// PLOTS OF ESTIMATED SIGNAL AND BACKGROUND
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;


   TGraphErrors *gr3 = new TGraphErrors(numbins,x,ys_c,ex,eys_c);




	gr3->Fit("gaus");
//   	cnv1->cd(); gr1->Draw("AB");
//   	cnv2->cd(); gr2->Draw("AB");

//      	TCanvas* cnv3 = new TCanvas("cnv3","cnv3") ;   	cnv3->cd();


   gr3->SetMarkerColor(kRed-3);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
     gr3->Draw("AP");

cnv1->cd();




   cnv1->Update();// 

}






