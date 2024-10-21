#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TRandom.h>

void corr_rand() {


TCanvas *c = new TCanvas("myCanvas","My Canvas",700,700);
TRandom *ran = new TRandom(450);
TH2F *h1 = new TH2F("my_hist1","Gauss",70,-700,700,70,-700,700);
TH2F *h2 = new TH2F("my_hist2","Gauss2",70,-700,700,70,-700,700);
double x,y,z;

h1->SetMarkerColor(4);
h2->SetMarkerColor(6);

for(int ii=0;ii<10000;ii++){
x=ran->Gaus(-100,150);
y=ran->Gaus(-100,50);
//Correlation coefficient
//corr_coeff=0.5;
//z=(corr_coeff*x)+((1-corr_coeff)*y);
z=x+y;
//h1->Fill(x,y);
h1->Fill(x,z);

x=ran->Gaus(100,150);
y=ran->Gaus(100,50);
z=x+y;
h2->Fill(x,z);
}


h1->Draw();
h2->Draw("sames");


}
