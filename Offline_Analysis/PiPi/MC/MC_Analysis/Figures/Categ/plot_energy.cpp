#include<cmath>
void plot_energy(){




    TCanvas *c1 = new TCanvas("myCanvas1","myCanvas1",800,800);
//    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",800,800);
//c1->Divide(1,2);
//    TCanvas *c2 = new TCanvas("myCanvas2","myCanvas2",800,800);


TCut t1="Deltam > 0.139 && Deltam < 0.143";
TCut t1p="Deltam > 0.15";
TCut t2="Pstdst > 2.9";

//TCut t3="Dmass > 1.7 && Dmass < 2.0";
TCut t3="Dmass > 1.82 && Dmass < 1.9";
TCut t3p="Dmass > 1.9";

TCut t4="Egam1s > 0.03 & Egam2s > 0.03";
TCut t5="(Gam1thet < -60 && Egamma1 > 0.150) || (Gam1thet > 73 && Egamma1 > 0.100) || (Gam1thet > -60 && Gam1thet < 73 && Egamma1 > 0.050)";
TCut t6="(Gam2thet < -60 && Egamma2 > 0.150) || (Gam2thet > 73 && Egamma2 > 0.100) || (Gam2thet > -60 && Gam2thet < 73 && Egamma2 > 0.050)";





TCut t7=t1&&t2&&t3&&t4&&t5&&t6;

TCut t7p=t1p&&t2&&t3p&&t4&&t5&&t6;

  TChain* chain=new TChain("h1");
  chain->Add("pipi_MC7.root");
  Int_t nevt=(int)chain->GetEntries();

 Float_t  f_Dstf,f_Egamma2,f_Egam1s;

  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Egamma2",&f_Egamma2);
  h1->SetBranchAddress("Egam1s",&f_Egam1s);

    TH1F *h11 = new TH1F("my_hist11","hist11",100,0.0,2.5);  //deltam
    TH1F *h2 = new TH1F("my_hist2","hist2",100,0.0,2.5);  //deltam

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);
if(f_Dstf==1){
h11->Fill(f_Egamma2);
h2->Fill(f_Egam1s);}
    }

	c1->cd(1);
h2->Draw();
h11->Draw("same");

/*
	//c1->cd(1);
	c1->cd();
	h1->SetLineColor(kBlue);
	h1->Draw("Pizmom","Dstf==1");
	//c1->cd(1);
	c2->cd();
	h1->Draw("Pizsmom","Dstf==1","same");
*/


	//c1->cd(1);
/*
	h1->SetLineColor(kBlue);
//	h1->SetTitle("Average energy of the hard pion decay photons");

	c1->cd();
	h1->Draw("Egamma2","Dstf==1");
	h1->Draw("Egam1s","Dstf==1","same");

  TH1F * myh1 =(TH1F *) htemp->Clone();
  myh1->SetName("myh1");


	h1->Draw("Egam1s","Dstf==1","same");
  TH1F * myh2 =(TH1F *) htemp->Clone();
  myh2->SetName("myh2");

	c1->cd();
  myh1->Draw(); 
  myh2->Draw("same"); 


  leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->AddEntry(myh1,"D #rightarrow K_{S} #pi","lep");
  leg->AddEntry(myh2,"D #rightarrow #pi #pi","lep");
  leg->Draw();
*/

//	h1->Draw("Egamma1","Dst//	
//	h1->Draw("Egamma2","Dstf==1","same");

//	h1->Draw("(Egamma2+Egamma2)/2","Dstf==1");

	//c1->cd(1);
//	c2->cd();
//	h1->SetTitle("Average energy of the soft pion decay photons");
//	h1->Draw("Egam2s","Dstf==1","same");
//	h1->Draw("Egamma2","Dstf==1");//*******

//gpad->Update();

}












