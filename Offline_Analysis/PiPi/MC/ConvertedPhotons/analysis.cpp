#include<cmath>
void analysis(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);
    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",700,700);

  TChain* chain=new TChain("h1");


//  chain->Add("DtoPiZPi_50000x2_Y4s_conv.root");
  chain->Add("DtoPiZPi_exp55_Data_conv.root");


	h1->SetLineColor(kBlue);
	
	c1->cd();
	h1->Draw("Deltam","Type == 2");
	c2->cd();
	h1->Draw("Dmass","Type == 2");




   	gPad->Update();

//	c1->SaveAs("a.eps");
//	c2->SaveAs("b.eps");



}

