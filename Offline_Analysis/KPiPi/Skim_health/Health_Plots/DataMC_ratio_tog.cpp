// Y4S and Y5S together
#include<cmath>
void DataMC_ratio_tog(){

int count=-1, count2=-1;
float exp[50], e_exp[50], lum_MC[50], e_lum_MC[50], lum_Data[50], e_lum_Data[50], ratio[50], e_ratio[50], sqsum;
float exp2[50], e_exp2[50], lum_MC2[50], e_lum_MC2[50], lum_Data2[50], e_lum_Data2[50], ratio2[50], e_ratio2[50], sqsum2;


ifstream fin;
fin.open("Results_GMC.txt");
while(!fin.eof()){
count++;
fin>>exp[count]>>lum_MC[count]>>e_lum_MC[count];
}//while
fin.close();

fin.open("Results_GMC_5s.txt");
while(!fin.eof()){
count2++;
fin>>exp2[count2]>>lum_MC2[count2]>>e_lum_MC2[count2];
}//while
fin.close();




count=-1;
fin.open("Results_data.txt");
while(!fin.eof()){
count++;
fin>>exp[count]>>lum_Data[count]>>e_lum_Data[count];
}//while
fin.close();


count2=-1;
fin.open("Results_data_5s.txt");
while(!fin.eof()){
count2++;
fin>>exp2[count2]>>lum_Data2[count2]>>e_lum_Data2[count2];
}//while
fin.close();


//---------------------------------------------------------------------------------------------------------------

for(int i=0;i<count;i++){
ratio[i]=lum_Data[i]/lum_MC[i];
sqsum=((e_lum_MC[i]/lum_MC[i])*(e_lum_MC[i]/lum_MC[i]))+((e_lum_Data[i]/lum_Data[i])*(e_lum_Data[i]/lum_Data[i]));
e_ratio[i]=ratio[i]*sqrt(sqsum);
e_exp[i]=0.0;
}


for(int i=0;i<count2;i++){
ratio2[i]=lum_Data2[i]/lum_MC2[i];
sqsum2=((e_lum_MC2[i]/lum_MC2[i])*(e_lum_MC2[i]/lum_MC2[i]))+((e_lum_Data2[i]/lum_Data2[i])*(e_lum_Data2[i]/lum_Data2[i]));
e_ratio2[i]=ratio2[i]*sqrt(sqsum2);
e_exp2[i]=0.0;
}

//---------------------------------------------------------------------------------------------------------------

//Graph and plot
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   TGraphErrors *gr = new TGraphErrors(count,exp,ratio,e_exp,e_ratio);
   TGraphErrors *gr2 = new TGraphErrors(count2,exp2,ratio2,e_exp2,e_ratio2);

   gr->SetMarkerColor(6);
   gr->SetMarkerStyle(21);

   gr2->SetMarkerColor(kRed);
   gr2->SetMarkerStyle(21);

     TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Data : MC ratio of yields for the mode D #rightarrow K #pi #pi");
     mg->Add(gr,"p");
     mg->Add(gr2,"p");
     mg->Draw("A");


  leg = new TLegend(0.6,0.7,0.89,0.89);
  leg->AddEntry(gr,"#Upsilon(4S)","lep");
  leg->AddEntry(gr2,"#Upsilon(5S)","lep");
  leg->Draw();



   c1->Update();

}






