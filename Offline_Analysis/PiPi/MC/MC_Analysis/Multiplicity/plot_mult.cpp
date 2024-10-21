#include<cmath>
void plot_mult(){





    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
//    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");



  TChain* chain2=new TChain("h1");


//  TChain* chain=new TChain("h3");
//  TChain* chain2=new TChain("h4");

//  chain->Add("pipi_1M_SMC.root");
//  chain2->Add("pipi_1M_SMC.root");
//  chain2->Add("pipi_MC0_Y5s.root");
//    chain2->Add("pipi_MC1_Y5s.root");

    chain2->Add("pipi_MC0_Y4s.root");
    chain2->Add("pipi_MC1_Y4s.root");


//c1->cd();
//h3->Draw("Tmult","Tmult > 0 && Tmult < 10");
//c2->cd();
//h4->Draw("Tmult","Tmult > 0");

//c1->cd();
//h4->Draw("Tmult","Tmult > 0 && Tmult < 10");



//To fing out BCS success rate
c1->cd();

//h1->Draw("Dstf","Multflag > 1 && Truevt ==1 && Dstf ==1");
//h1->Draw("Dstf","Multflag == 1 && Truevt ==1");



//h1->Draw("Dstf","Multflag > 1 && Truevt ==1");
h1->Draw("Dstf","Multflag > 1 && Truevt !=1");


}












