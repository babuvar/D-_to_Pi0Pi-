#include<cmath>
void BkgSrc2(){





//    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");
//    TCanvas *c2 = new TCanvas("myCanvas2","My Canvas2");
//    TCanvas *c3 = new TCanvas("myCanvas3","My Canvas3");

TCut T="Did == 421 || Did == -421";

Float_t f_Did, f_Dau1id, f_Dau2id, f_Dau3id;


  TChain* chain=new TChain("h1");
  chain->Add("pipi_charm_BkgSrc.root");
  Int_t nevt=(int)chain->GetEntries();
  cout<<"nevt\t"<<nevt <<endl;



  h1->SetBranchAddress("Did",&f_Did);
  h1->SetBranchAddress("Dau1id",&f_Dau1id);
  h1->SetBranchAddress("Dau2id",&f_Dau2id);
  h1->SetBranchAddress("Dau3id",&f_Dau3id);

  for(int i=0;i<nevt;i++) {
      chain->GetEntry(i);
if((f_Did==421 || f_Did==-421) && f_Dau1id != 0){ 
//if((f_Did==411 || f_Did==-411) && f_Dau1id != 0){


cout<<f_Dau1id<<"\t"<<f_Dau2id<<"\t"<<f_Dau3id<<endl;}


}


/*
c1->cd();
h1->Draw("Dau1id",T);

c2->cd();
h1->Draw("Dau2id",T);

c3->cd();
h1->Draw("Dau3id",T);
*/
}












