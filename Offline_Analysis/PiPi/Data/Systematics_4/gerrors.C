void gerrors() {

   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

   c1->GetFrame()->SetBorderSize(12);

   const Int_t n = 5;

   Float_t y[n]  = {0.00522197,0.00355482 , 0.00439223,0.00603567 , 0.0082023};
   Float_t ey[n]  = {0.0193663, 0.019442, 0.0194512, 0.0192487, 0.0195775};


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
