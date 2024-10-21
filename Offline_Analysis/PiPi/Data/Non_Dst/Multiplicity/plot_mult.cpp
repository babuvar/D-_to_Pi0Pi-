#include<cmath>
void plot_mult(){



    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
    c1->cd();



//    TChain* chain2=new TChain("h4");
//    TChain* chain2=new TChain("h5");

//    TChain* chain2=new TChain("h2");
    TChain* chain2=new TChain("h1");


    chain2->Add("pipi_mc.root");
//    chain2->Add("ks_mc.root");


//h5->Draw("Tmult","Tmult > 0 && Tmult < 10"); //Un-Tagged
//h5->Draw("Tmult","Tmult > 1 && Tmult < 10"); //Un-Tagged

//h4->Draw("Tmult","Tmult > 0 && Tmult < 10"); //Tagged
//h4->Draw("Tmult","Tmult > 1 && Tmult < 10"); //Tagged



//To fing out BCC success rate
//c1->cd();

//h1->Draw("Dstf","Multflag > 1 && Truevt ==1 && Dstf ==1");
h1->Draw("Dstf","Multflag > 1 && Truevt ==1");


//h2->Draw("Df","Multflag > 1 && Truevt ==1 && Df ==1");
//h2->Draw("Df","Multflag > 1 && Truevt ==1");


}












