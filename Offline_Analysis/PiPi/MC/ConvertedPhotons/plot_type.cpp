#include<cmath>
void plot_type(){

    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");



  TChain* chain=new TChain("h1");


  chain->Add("DtoPiZPi_50000x2_Y4s_conv.root");



c1->cd();
h1->Draw("Type");

}












