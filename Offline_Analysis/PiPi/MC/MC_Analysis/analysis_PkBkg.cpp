#include<cmath>
void analysis_PkBkg(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);
    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",700,700);
    TCanvas *c3 = new TCanvas("myCanvas3","myCanvas3",700,700);
    TCanvas *c4 = new TCanvas("myCanvas4","myCanvas4",700,700);
    TCanvas *c5 = new TCanvas("myCanvas5","myCanvas5",700,700);


  TChain* chain=new TChain("h1");
  chain->Add("pipi_MC6.root");
  chain->Add("pipi_MC7.root");


TCut t1="Deltam > 0.139 && Deltam < 0.142";
//TCut t1="Deltam > 0.148 && Deltam < 0.155";
TCut t2="Dmass > 1.70  && Dmass < 2.0";
TCut t3="fabs(Pizsmass - 0.135) < 0.012";
TCut t4="fabs(Pizmass - 0.135) < 0.012";
TCut t5="Pizmom > 0.95";
TCut t6="Pimom > 0.74";
TCut t7="Pstdst > 2.9";
//Photon cuts
TCut t8="(Gam1hthe < -60 && Egamma1 > 0.150) || (Gam1hthe > 73 && Egamma1 > 0.100) || (Gam1hthe > -60 && Gam1hthe < 73 && Egamma1 > 0.050)";
TCut t9="(Gam2hthe < -60 && Egamma2 > 0.150) || (Gam2hthe > 73 && Egamma2 > 0.100) || (Gam2hthe > -60 && Gam2hthe < 73 && Egamma2 > 0.050)";
//--------------------------
TCut t10="(Categ == 5 && Egam1s > 0.044) || (Categ == 2 && Egam1s > 0.066) || (Categ == 6 && Egam1s > 0.044)";
TCut t11="(Categ == 5 && Egam2s > 0.044) || (Categ == 2 && Egam2s > 0.036) || (Categ == 6 && Egam2s > 0.054)";


TCut T=t1&&t2&&t3&&t4&&t5&&t6&&t7&&t8&&t9&&t10&&t11;


	h1->SetLineColor(kBlue);
//	h1->SetLineColor(kRed);
	
	c1->cd();
//	h1->Draw("Df","Df != 1");



//	h1->Draw("Dmass",T);

//	h1->Draw("Dmass",T&&"Df != 1");

	h1->Draw("Dmass",T&&"Df == 0");
	c2->cd();
	h1->SetLineColor(kRed);
	h1->Draw("Dmass",T&&"Df == -1");

	c3->cd();
	h1->SetLineColor(kGreen);
	h1->Draw("Dmass",T&&"Df == 3");

	c4->cd();
	h1->SetLineColor(kYellow);
	h1->Draw("Dmass",T&&"Df == -11");

	c5->cd();
	h1->SetLineColor(kPink);
	h1->Draw("Dmass",T&&"Df != 1 && Df != 0 && Df != -1 && Df != -11 && Df != 3 && Df != 10");


//	h1->Draw("Df",T&&"Df != 1 && Df != 0 && Df != -1 && Df != -11 && Df != 3");


/*
	h1->Draw("Dmass",T&&"Df == 10");

	c2->cd();
	h1->SetLineColor(kRed);
	h1->Draw("Dmass",T&&"Df == 21");

	c3->cd();
	h1->SetLineColor(kGreen);
	h1->Draw("Dmass",T&&"Df == 5");
*/


//	h1->Draw("Dmass",T&&"Df == 1");
//	h1->Draw("Dmass",T&&"Df == 3");


   	gPad->Update();


}

