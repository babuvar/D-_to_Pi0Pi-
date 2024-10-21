void gerrors() {

   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

   c1->GetFrame()->SetBorderSize(12);

   const Int_t n = 5;
//   Float_t y[n]  = {1.3415,1.38903,1.52425,1.67227,1.76009,1.8673,1.82663,2.17421};
//   Float_t ey[n]  = {1.94106,1.7456,1.78908,1.97322,2.0031,2.00878,2.01167,1.98092};

   Float_t y[n]  = {1.38188,1.40174,1.53851,1.686,1.7479};
   Float_t ey[n]  = {1.95787,1.97082,1.9842,1.96921,1.97861};


   Float_t x[n] = {1.68,1.69,1.70,1.71,1.72};
   Float_t ex[n] = {0,0,0,0,0};
   TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);

//  TF1 *f = new TF1("f", "[0]");
//  gr->Fit(f);



   gr->SetTitle("Measured A_{CP}^{uncorr}");
   gr->SetMarkerColor(9);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");




   gr->SetMarkerSize(1.4);


//  leg = new TLegend(0.6,0.7,0.89,0.89);
//  leg->AddEntry(gr,"MC (1-D Fit)","lep");
//  leg->Draw();



   c1->Update();
}
