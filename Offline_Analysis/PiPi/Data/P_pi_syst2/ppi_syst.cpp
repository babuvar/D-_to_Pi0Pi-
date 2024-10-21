#include<fstream>
#include<string>
void ppi_syst(void)
{
float st, sum=0;
int i, j, count=-1;
float a[15][15], b[15][15];

ifstream fin;
fin.open("Map_excess.txt");
while(!fin.eof())
{
count++;
fin>>st;
j=count%6;
i=floor(count/6);
a[i][j]=st;
}
fin.close();

//--------------------------------------------
count=-1;
fin.open("Map_Asy.txt");
while(!fin.eof())
{
count++;
fin>>st;
j=count%6;
i=floor(count/6);
b[i][j]=st;
}
fin.close();


for(i=0; i<=9; i++){
for(j=0; j<=5; j++){
//sum=sum+(b[i][j]*a[i][j]);
sum=sum+(b[i][j]*a[i][j]);
}
cout<<endl;}

cout<<"sum = "<<sum/100<<endl;

}













