#include<cmath>
void DataMC_ratio(){

int count=-1;
float exp[50], e_exp[50], lum_MC[50], e_lum_MC[50], lum_Data[50], e_lum_Data[50], ratio[50], e_ratio[50], sqsum;

ifstream fin;
//fin.open("Results_GMC.txt");
fin.open("Results_GMC_5s.txt");
while(!fin.eof()){
count++;
fin>>exp[count]>>lum_MC[count]>>e_lum_MC[count];
}//while
fin.close();
count=-1;

//fin.open("Results_data.txt");
fin.open("Results_data_5s.txt");
while(!fin.eof()){
count++;
fin>>exp[count]>>lum_Data[count]>>e_lum_Data[count];
}//while
fin.close();




for(int i=0;i<count;i++){
ratio[i]=lum_Data[i]/lum_MC[i];
sqsum=((e_lum_MC[i]/lum_MC[i])*(e_lum_MC[i]/lum_MC[i]))+((e_lum_Data[i]/lum_Data[i])*(e_lum_Data[i]/lum_Data[i]));
e_ratio[i]=ratio[i]*sqrt(sqsum);
e_exp[i]=0.0;
}


//Graph and plot
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   TGraphErrors *gr = new TGraphErrors(count,exp,ratio,e_exp,e_ratio);
   gr->SetTitle("Data : MC ratio of yields for the mode D #rightarrow K #pi #pi");
   gr->SetMarkerColor(6);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");

//   gr->GetYaxis()->SetRangeUser(0.0,2300.0);


   c1->Update();

}






