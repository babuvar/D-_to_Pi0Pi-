void plot_eff(){
int count=-1;
float a[50], b[50],c[50],d[50],e_a[50],dummy;
float a2[50], b2[50],c2[50],d2[50],e_a2[50];

ifstream fin;
//fin.open("Results_data.txt");
//fin.open("Results_GMC.txt");
//fin.open("Results_data_5s.txt");
//fin.open("Results_GMC_5s.txt");
fin.open("Results_data_continuum.txt");
while(!fin.eof()){
count++;
fin>>a[count]>>b[count]>>c[count];
}//while
fin.close();
count=-1;
//fin.open("lum_table.txt");
//fin.open("lum_table_5s.txt");
fin.open("lum_table_continuum.txt");

while(!fin.eof()){
count++;
fin>>dummy>>d[count];
}//while
fin.close();




for(int i=0;i<count;i++){
//cout<<a[i]<<"\t"<<b[i]<<"\t"<<c[i]<<"\t"<<d[i]<<endl;
b[i]=(b[i]*1000)/d[i];
c[i]=(c[i]*1000)/d[i];

e_a[i]=0.5;

//cout<<a[i]<<"\t"<<b[i]<<"\t"<<c[i]<<endl;
}


//Graph and plot
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   TGraphErrors *gr = new TGraphErrors(count,a,b,e_a,c);
//   gr->SetTitle("Yield/fb^{-1} by experiment for the mode D #rightarrow K #pi #pi (Data)");
   gr->SetTitle("Yield/fb^{-1} by experiment for the mode D #rightarrow K #pi #pi (MC)");
//   gr->SetMarkerColor(9);
   gr->SetMarkerColor(6);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
//   gr->GetYaxis()->SetRangeUser(0.0,1900.0);
   gr->GetYaxis()->SetRangeUser(0.0,2300.0);
//   gr->Draw("AP");

   c1->Update();

}






