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

void fit_Mpis_data2()
{


float numbins=50;
float binwidth=0.050/numbins;

for(int i=0;  i<=99; i++){x[i]=0.110+((i+0.5)*binwidth); ex[i]=0;}

 

//Text file with results (counted)
ifstream fin;
fin.open("Results_Signal.txt");
//fin.open("Results_Signal_data.txt");
for(int i =0; i < numbins; i++){
fin>>x[i]>>ys_c[i]>>eys_c[i];
}
fin.close();





// PLOTS OF ESTIMATED SIGNAL AND BACKGROUND
      	TCanvas* cnv1 = new TCanvas("cnv1","cnv1") ;


   TGraphErrors *gr3 = new TGraphErrors(numbins,x,ys_c,ex,eys_c);
Double_t par[9];

/*
    // Fitting using user-defined functions

    TF1 *func = new TF1("myGauss","[0]*exp(-0.5*((x-[1])/[2])**2)",0.116,0.156);
    func->SetParNames("Factor","Mean","Sigma");
    func->SetParameter(0,100000);
    func->SetParameter(1,0.135);
    func->SetParameter(2,0.005);
*/

//    TF1 *func = new TF1("myGauss","[0]*([1]*exp(-0.5*((x-[2])/[3])**2) + (1-[1])*exp(-0.5*((x-[2])/[4])**2))",0.11,0.16);
//    TF1 *func = new TF1("myGauss","[0]*(([1]*exp(-0.5*((x-[2])/[3])**2)) + ((1.0-[1])*exp(-0.5*((x-[2])/[4])**2)))",0.116,0.156);
    TF1 *func = new TF1("myGauss","([0]*exp(-0.5*((x-[2])/[3])**2)) +( [1]*exp(-0.5*((x-[2])/[4])**2))",0.115,0.155);



    // parameter names
    func->SetParNames("Factor1","Factor2","Mean","Sigma1","Sigma2");
    //Set Parameters:


    func->SetParameter(0,90000);
    func->SetParLimits(0,0,100000);
    func->SetParameter(1,40000);
    func->SetParLimits(1,0,50000);
    func->SetParameter(2,0.135);
    func->SetParameter(3,0.0055);
    func->SetParLimits(3,0.001,0.01);
    func->SetParameter(4,0.007);
    func->SetParLimits(4,0.005,0.03);

    gr3->Fit(func,"R");// use option "R" to restrict to a certain region for fitting





cnv1->cd();

   gr3->SetMarkerColor(kRed-3);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
     gr3->Draw("AP");

  func->Draw("same");


}






