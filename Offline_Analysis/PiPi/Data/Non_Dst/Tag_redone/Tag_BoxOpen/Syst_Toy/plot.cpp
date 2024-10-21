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


  TH1D* Araw_res=new TH1D("Araw_res", "Araw_res", 50, 0.3, 0.8);



void plot()
{

int count = -1;

double araw_value[6000];

ifstream fin;
//fin.open("full.txt");
fin.open("full_wCorr.txt");
while(!fin.eof()){
count++;
fin>>araw_value[count];
}
fin.close();


for(int i=0; i<5000 ; i++){
//cout<<araw_value[i]<<endl;
Araw_res->Fill(araw_value[i]);

}




    TF1 *func = new TF1("myGauss","([0]*exp(-0.5*((x-[2])/[3])**2)) +( [1]*exp(-0.5*((x-[2])/[4])**2))",0.3,0.8);
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



}











