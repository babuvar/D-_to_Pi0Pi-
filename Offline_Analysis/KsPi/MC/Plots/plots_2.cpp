#include<cmath>
void plots_2(){


    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",700,700);
    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",700,700);
    TCanvas *c3 = new TCanvas("myCanvas3","myCanvas3",700,700);
//    TCanvas *c4 = new TCanvas("myCanvas4","myCanvas4",700,700);

  TChain* chain=new TChain("h1");
//  chain->Add("KsPi_GMC0_4s.root");
//  chain->Add("KsPi_GMC1_4s.root");
//  chain->Add("SMC_100K_kspi.root");

  chain->Add("KsPi_GMC0_4s.root");
  chain->Add("KsPi_GMC1_4s.root");
Int_t nevt=(int)chain->GetEntries();

//  chain->Add("pipi_MC6.root");
//  chain->Add("pipi_MC7.root");
//  chain->Add("SMC_100K_pipi.root");
//  chain->Add("SMC_1M_pipi.root");

  chain->Add("pipi_MC0_Y4s.root");
  chain->Add("pipi_MC1_Y4s.root");
  chain->Add("pipi_SMC_2M_4s.root");
Int_t nevt2=(int)chain->GetEntries();



cout<<"nevt = "<<nevt<<endl;
cout<<"nevt2 = "<<nevt2<<endl;


TCut t1="Deltam > 0.139 && Deltam < 0.142";
TCut t2="Dmass >1.80  && Dmass < 1.94";
TCut t3="Pizsmass > 0.125  && Pizsmass < 0.143";
TCut t4="Ksmom > 1.06";
TCut t5="Pimom > 0.84";
TCut t6="Pstdst > 2.95";
TCut t7="(Categ == 5 && Egam1s > 0.046) || (Categ == 2 && Egam1s > 0.068) || (Categ == 6 && Egam1s > 0.030)";
TCut t8="(Categ == 5 && Egam2s > 0.046) || (Categ == 2 && Egam2s > 0.036) || (Categ == 6 && Egam2s > 0.044)";
TCut t12="Df == 1";

TCut t2p="Dmass >1.70  && Dmass < 2.0";
TCut t9="Pizmass > 0.119  && Pizmass < 0.151";
TCut t4p="Pizmom > 1.06";
TCut t10="(Gam1hthe < -60 && Egamma1 > 0.150) || (Gam1hthe > 73 && Egamma1 > 0.100) || (Gam1hthe > -60 && Gam1hthe < 73 && Egamma1 > 0.050)";
TCut t11="(Gam2hthe < -60 && Egamma2 > 0.150) || (Gam2hthe > 73 && Egamma2 > 0.100) || (Gam2hthe > -60 && Gam2hthe < 73 && Egamma2 > 0.050)";


TCut T=t1&&t2&&t3&&t4&&t5&&t6&&t7&&t8&&t12;

TCut Tp=t1&&t2p&&t3&&t4p&&t5&&t6&&t7&&t8&&t9&&t10&&t11&&t12;





//-------------------------------------------------------

	h1->SetLineColor(kBlue);
	c1->cd();
	h1->Draw("Pimom",T,"",nevt,0);
//	h1->Draw("Picthet",T,"",nevt,0);
//	Cloning
  TH1F * myh1 =(TH1F *) htemp->Clone();
  myh1->SetName("myh1");
//  myh1->Rebin(5);
Double_t norm = myh1->GetEntries();
  myh1->SetOption("E");
  myh1->Sumw2();
  myh1->Scale(1/norm);
	myh1->SetLineColor(kBlue);       

//-------------------------------------------------------




	c2->cd();
	h1->Draw("Pimom",Tp,"",(nevt2-nevt),nevt);
//	h1->Draw("Picthet",Tp,"",(nevt2-nevt),nevt);
//	Cloning
  TH1F * myh2 =(TH1F *) htemp->Clone();
  myh2->SetName("myh2");
//  myh2->Rebin(5);
Double_t norm2 = myh2->GetEntries();
  myh2->SetOption("E");
  myh2->Sumw2();
  myh2->Scale(1/norm2);
	myh2->SetLineColor(kRed);       

//-------------------------------------------------------



//Draw
  c3->cd();
  myh2->Draw(); 
  myh1->Draw("same"); 
//  myh1->Draw(); 
//  myh2->Draw("same"); 


  leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->AddEntry(myh1,"D #rightarrow K_{S} #pi","lep");
  leg->AddEntry(myh2,"D #rightarrow #pi #pi","lep");
  leg->Draw();




   	gPad->Update();


}
