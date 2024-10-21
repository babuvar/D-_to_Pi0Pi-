#include<cmath>
void BkgSrc(){





    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
//    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");

TCut t1="Did == 421 || Did == -421";
TCut t2="Did == 411 || Did == -411";
//TCut t=t1||t2;
TCut t3="Dau1id != 0  && Dau1id > -1000 && Dau1id < 1000";
//TCut t3="Dau1id != 0 && Dau1id > 205 && Dau1id < 215";
TCut t4="Dcharge == 1";
TCut t4p="Dcharge == -1";


//TCut T=t1&&t3&&t4;
TCut T=t2&&t3&&t4;
//TCut Tp=t1&&t3&&t4p;
TCut Tp=t2&&t3&&t4p;

  TChain* chain=new TChain("h1");


  chain->Add("pipi_charm_BkgSrc.root");


//c1->cd();
//h1->Draw("Dmass",T);





c1->cd();
h1->SetLineColor(kBlue);
h1->Draw("Dau1id",Tp);
h1->SetLineColor(kRed);
h1->Draw("Dau1id",T,"same");




c2->cd();
h1->SetLineColor(kBlue);
h1->Draw("Dau2id",Tp);
h1->SetLineColor(kRed);
h1->Draw("Dau2id",T,"same");



/*
c3->cd();
h1->SetLineColor(kBlue);
h1->Draw("Dau3id",T);
h1->SetLineColor(kRed);
h1->Draw("Dau3id",Tp,"same");
*/


}












