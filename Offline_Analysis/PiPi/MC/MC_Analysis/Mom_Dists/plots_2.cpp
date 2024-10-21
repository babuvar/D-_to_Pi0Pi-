#include<cmath>
void plots_2(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);
    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",700,700);
    TCanvas *c3 = new TCanvas("myCanvas3","myCanvas3",700,700);


  TChain* chain=new TChain("h1");
  chain->Add("PiPi_1M.root");



//-------------------------------------------------------

	h1->SetLineColor(kBlue);
	c1->cd();
//	h1->Draw("Pimom","Dstf==1 && Pimom > 0.75 && Pizmom > 0.75");
	h1->Draw("Pimom","Dstf==1");
//	h1->Draw("Genpimom","Dstf==1");
  TH1F * myh1 =(TH1F *) htemp->Clone();
  myh1->SetName("myh1");
  myh1->SetOption("E");
  myh1->Sumw2();
	myh1->SetLineColor(kBlue);       

//-------------------------------------------------------




	c2->cd();
//	h1->Draw("Pizmom","Dstf==1 && Pimom > 0.75 && Pizmom > 0.75");
	h1->Draw("Pizmom","Dstf==1");
//	h1->Draw("Genpizmo","Dstf==1");
  TH1F * myh2 =(TH1F *) htemp->Clone();
  myh2->SetName("myh2");
  myh2->SetOption("E");
  myh2->Sumw2();
	myh2->SetLineColor(kRed);       

//-------------------------------------------------------



//Draw
  c3->cd();
  myh1->Draw(); 
  myh2->Draw("same"); 


  leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->AddEntry(myh1,"P_{#pi^{#pm}}","lep");
  leg->AddEntry(myh2,"P_{#pi^{0}}","lep");
  leg->Draw();


   	gPad->Update();


}

