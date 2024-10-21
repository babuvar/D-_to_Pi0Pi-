#include<cmath>
void compare(){




    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",800,800);
c1->Divide(2,2);
    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",800,800);
c2->Divide(2,2);

TCut t1="Deltam > 0.139 && Deltam < 0.143";
TCut t1p="Deltam > 0.15";
TCut t2="Pstdst > 2.9";


TCut t3="Dmass > 1.7 && Dmass < 2.0";
//TCut t3="Dmass > 1.82 && Dmass < 1.9";
TCut t3p="Dmass > 1.9";

TCut t4="Egam1s > 0.03 & Egam2s > 0.03";
TCut t5="(Gam1thet < -60 && Egamma1 > 0.150) || (Gam1thet > 73 && Egamma1 > 0.100) || (Gam1thet > -60 && Gam1thet < 73 && Egamma1 > 0.050)";
TCut t6="(Gam2thet < -60 && Egamma2 > 0.150) || (Gam2thet > 73 && Egamma2 > 0.100) || (Gam2thet > -60 && Gam2thet < 73 && Egamma2 > 0.050)";

TCut t7="Type == 2";



TCut t8=t1&&t2&&t3&&t4&&t5&&t6&&t7;

TCut t8p=t1p&&t2&&t3p&&t4&&t5&&t6;

TCut t9="Conmas2 < 0.002";
TCut t10=t8&&t9;

  TChain* chain=new TChain("h1");
//  chain->Add("DtoPiZPi_50000x2_Y4s_conv.root");
  chain->Add("DtoPiZPi_exp55_Data_conv.root");

	h1->SetLineColor(kBlue);


	c1->cd(1);
	h1->Draw("Dmass",t8);
	c1->cd(2);
	h1->Draw("Deltam",t8);
	c1->cd(3);
	h1->Draw("Conmas2",t10);
	c1->cd(4);
	h1->Draw("Pstdst",t8);


	c2->cd(1);
	h1->Draw("Pimom",t8);
	c2->cd(2);
	h1->Draw("Pizmom",t8);
	c2->cd(3);
	h1->Draw("Pstpizs",t8);
	c2->cd(4);
	h1->Draw("Pizmass",t8);



   	gPad->Update();





}

