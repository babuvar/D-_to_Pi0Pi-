#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;
void MCTruth(void)
{
    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");


    TH1F *h11 = new TH1F("my_hist11","hist11",50,0.135,0.155);  //deltam
    TH1F *h2 = new TH1F("my_hist2","hist2",50,0.135,0.155);  //deltam
    TH1F *h3 = new TH1F("my_hist3","hist3",50,0.135,0.155);  //deltam
    TH1F *ha = new TH1F("my_hista","hista",50,0.135,0.155);  //deltam
    TH1F *hd = new TH1F("my_histd","histd",50,0.135,0.155);  //deltam
//    TH1F *h4 = new TH1F("my_hist4","hist4",50,0.135,0.155);  //deltam

    TH1F *h5 = new TH1F("my_hist5","hist5",50,1.70,2.0);  //d0 mass
    TH1F *h6 = new TH1F("my_hist6","hist6",50,1.70,2.0);  //d0 mass
    TH1F *h7 = new TH1F("my_hist7","hist7",50,1.70,2.0);  //d0 mass
    TH1F *hb = new TH1F("my_histb","histb",50,1.70,2.0);  //d0 mass
    TH1F *he = new TH1F("my_histe","histe",50,1.70,2.0);  //d0 mass
//    TH1F *h8 = new TH1F("my_hist8","hist8",50,1.70,2.0);  //d0 mass


    TH1F *h8 = new TH1F("my_hist8","hist5",50,2.5,5.0);  //P*(D*)
    TH1F *h9 = new TH1F("my_hist9","hist6",50,2.5,5.0);  //P*(D*)
    TH1F *h10 = new TH1F("my_hist10","hist7",50,2.5,5.0);  //P*(D*)
    TH1F *hc = new TH1F("my_histc","histc",50,2.5,5.0);  //P*(D*)
    TH1F *hf = new TH1F("my_histf","histf",50,2.5,5.0);  //P*(D*)


  h11->SetFillColor(4);
  h2->SetFillColor(8);
  h3->SetFillColor(2);
//  h4->SetFillColor(6);
  ha->SetFillColor(6);
  hd->SetFillColor(5);


  h5->SetFillColor(4);
  h6->SetFillColor(8);
  h7->SetFillColor(2);
//  h8->SetFillColor(6);
  hb->SetFillColor(6);
  he->SetFillColor(5);


  h8->SetFillColor(4);
  h9->SetFillColor(8);
  h10->SetFillColor(2);
  hc->SetFillColor(6);
  hf->SetFillColor(5);



  //LOAD DATA FILE
  TChain* chain=new TChain("h1");

//  chain->Add("DtoPiZPi_50000x2_Y4s_noBCS.root");
  chain->Add("pipi_gmc55_noBCS.root");
//  chain->Add("DtoPiZPi_50000x2_Y4s_wBCS.root");
//  chain->Add("pipi_gmc55_wBCS.root");



  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;

 Float_t  f_Deltam, f_Dmass, f_Dstf, f_Df, f_Dstf, f_PD, f_Did, f_Dstid;

  h1->SetBranchAddress("Deltam",&f_Deltam);
  h1->SetBranchAddress("Dmass",&f_Dmass);
  h1->SetBranchAddress("Pstdst",&f_PD);
  h1->SetBranchAddress("Dstf",&f_Dstf);
  h1->SetBranchAddress("Df",&f_Df);
  h1->SetBranchAddress("Dstid",&f_Dstid);
  h1->SetBranchAddress("Did",&f_Did);


  THStack *hs1 = new THStack("hs1","M_{D*} - M_{D}");
  THStack *hs2 = new THStack("hs2","D Mass");
  THStack *hs3 = new THStack("hs3","P*_{D*}");
 

  Int_t nevt_ac_p(0);
  Int_t nevt_ac_n(0);
  for(int i=0;i<nevt;i++) 
    {
      chain->GetEntry(i);

if(f_Deltam > 0.135 && f_Deltam < 0.155){
if(f_Dmass > 1.70 && f_Dmass < 2.70){
/*
if(f_Dstf == 1){h11->Fill(f_Deltam); h5->Fill(f_Dmass); h8->Fill(f_PD);}
else{
if(f_Df == 1){h2->Fill(f_Deltam); h6->Fill(f_Dmass); h9->Fill(f_PD);}
else{h3->Fill(f_Deltam); h7->Fill(f_Dmass); h10->Fill(f_PD);}
}

}}
*/
if(f_Dstf == 1){h11->Fill(f_Deltam); h5->Fill(f_Dmass); h8->Fill(f_PD);}
else if(f_Df == 1){h2->Fill(f_Deltam); h6->Fill(f_Dmass); h9->Fill(f_PD);}
else if(f_Did == 411 || f_Did == -411  ){ ha->Fill(f_Deltam); hb->Fill(f_Dmass); hc->Fill(f_PD); }
else if(f_Did == 421 || f_Did == -421  ){ hd->Fill(f_Deltam); he->Fill(f_Dmass); hf->Fill(f_PD); }
else{ h3->Fill(f_Deltam); h7->Fill(f_Dmass); h10->Fill(f_PD); }


}}

    }




  hs1->Add(ha);
  hs1->Add(hd);
  hs1->Add(h2);
  hs1->Add(h3);
  hs1->Add(h11);


  hs2->Add(h7);
  hs2->Add(hb);
  hs2->Add(he);
  hs2->Add(h6);
  hs2->Add(h5);


  hs3->Add(h10);
  hs3->Add(h9);
  hs3->Add(h8);
  hs3->Add(hc);
  hs3->Add(hf);


  leg = new TLegend(0.6,0.7,0.89,0.89);
//  leg = new TLegend();
  leg->SetHeader("Reconstructed Object");
  leg->AddEntry(h11,"Signal","f");
  leg->AddEntry(h2,"Wrong slow #pi^{0}","f");
  leg->AddEntry(h3,"combinatorial","f");
  leg->AddEntry(ha,"D+ contamination","f");
  leg->AddEntry(hd,"D0 contamination","f");

c1->cd(); hs1->Draw();  leg->Draw();
c2->cd(); hs2->Draw();  leg->Draw();
c3->cd(); hs3->Draw();  leg->Draw();




}

