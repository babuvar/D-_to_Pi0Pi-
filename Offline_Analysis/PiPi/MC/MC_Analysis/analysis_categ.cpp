#include<cmath>
void analysis_categ(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);


  TChain* chain=new TChain("h1");
  chain->Add("pipi_MC6.root");




//	h1->SetLineColor(kBlue);
	h1->SetLineColor(kRed);
	
	c1->cd();
//	h1->Draw("Categ","Df == 1 && Categ== 1");
	h1->Draw("Categ","Df != 1 && Categ== 5");


   	gPad->Update();


}

