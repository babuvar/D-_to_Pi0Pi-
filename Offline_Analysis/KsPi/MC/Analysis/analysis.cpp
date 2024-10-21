#include<cmath>
void analysis(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);


  TChain* chain=new TChain("h1");
  chain->Add("KsPi_GMC0_4s.root");


TCut t1="Deltam > 0.139 && Deltam < 0.142";
TCut t2="Dmass >1.80  && Dmass < 1.94";
TCut t3="Pizsmass > 0.125  && Pizsmass < 0.143";
TCut t4="Ksmom > 1.06";
TCut t5="Pimom > 0.84";
TCut t6="Pstdst > 2.95";
TCut t7="(Categ == 5 && Egam1s > 0.046) || (Categ == 2 && Egam1s > 0.068) || (Categ == 6 && Egam1s > 0.030)";
TCut t8="(Categ == 5 && Egam2s > 0.046) || (Categ == 2 && Egam2s > 0.036) || (Categ == 6 && Egam2s > 0.044)";

//TCut t9="Df==1";
TCut t9="Df!=1";

TCut T=t1&&t2&&t3&&t4&&t5&&t6&&t7&&t8&&t9;


	h1->SetLineColor(kBlue);
//	h1->SetLineColor(kRed);
	
	c1->cd();
	h1->Draw("Dmass",T);


   	gPad->Update();


}

