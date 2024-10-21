// 4S and 5S together

void plot_eff_tog(){
int count=-1, count2=-1 ;
float a[50], b[50],c[50],d[50],e_a[50],dummy;
float a2[50], b2[50],c2[50],d2[50],e_a2[50];

ifstream fin;


//fin.open("Results_data.txt");
fin.open("Results_GMC.txt");
while(!fin.eof()){
count++;
fin>>a[count]>>b[count]>>c[count];
}//while
fin.close();



//fin.open("Results_data_5s.txt");
fin.open("Results_GMC_5s.txt");
while(!fin.eof()){
count2++;
fin>>a2[count2]>>b2[count2]>>c2[count2];
}//while
fin.close();
//-------------------------------------------------

count=-1;
fin.open("lum_table.txt");
while(!fin.eof()){
count++;
fin>>dummy>>d[count];
}//while
fin.close();

count2=-1;
fin.open("lum_table_5s.txt");
while(!fin.eof()){
count2++;
fin>>dummy>>d2[count2];
}//while
fin.close();





for(int i=0;i<count;i++){
//cout<<a[i]<<"\t"<<b[i]<<"\t"<<c[i]<<"\t"<<d[i]<<endl;
b[i]=(b[i]*1000)/d[i];
c[i]=(c[i]*1000)/d[i];
e_a[i]=0.5;
//cout<<a[i]<<"\t"<<b[i]<<"\t"<<c[i]<<endl;
}


for(int i=0;i<count2;i++){
b2[i]=(b2[i]*1000)/d2[i];
c2[i]=(c2[i]*1000)/d2[i];
e_a2[i]=0.5;
}



//Graph and plot
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   TGraphErrors *gr = new TGraphErrors(count,a,b,e_a,c);
   TGraphErrors *gr2 = new TGraphErrors(count2,a2,b2,e_a2,c2);

   gr->SetMarkerColor(6);
   gr->SetMarkerStyle(21);

   gr2->SetMarkerColor(kRed);
   gr2->SetMarkerStyle(21);

     TMultiGraph *mg = new TMultiGraph();
//     mg->SetTitle("Yield/fb^{-1} by experiment for the mode D #rightarrow K #pi #pi (Data)");
     mg->SetTitle("Yield/fb^{-1} by experiment for the mode D #rightarrow K #pi #pi (MC)");
     mg->Add(gr,"p");
     mg->Add(gr2,"p");
     mg->Draw("A");

  leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->AddEntry(gr,"#Upsilon(4S)","lep");
  leg->AddEntry(gr2,"#Upsilon(5S)","lep");
  leg->Draw();



//   gr->GetYaxis()->SetRangeUser(0.0,1900.0);
//   gr->GetYaxis()->SetRangeUser(0.0,2300.0);


   c1->Update();

}






