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
#include <fstream>

using namespace RooFit;
const int num =100000;

double mean(double x[num]);
double stdev(double x[num]);
double cov(double x1[num],double x2[num]);
double corr(double x1[num],double x2[num]);


void corr_test(void)
{
double x1[num], x2[num], x3[num], x4[num], x5[num];
double y1[num], y2[num], y3[num], y4[num], y5[num];

int count=-1;
int i,j;
int dim=5;



TMatrixD a(dim,dim);

TRandom *ran = new TRandom(50);
TRandom *ran1 = new TRandom(150);
TRandom *ran2 = new TRandom(250);
TRandom *ran3 = new TRandom(355);
TRandom *ran4 = new TRandom(450);

for(int ii=0;ii<num;ii++){

 x1[ii]=ran->Gaus(0.0,1.0);
 x2[ii]=ran1->Gaus(0.0,1.0);
 x3[ii]=ran2->Gaus(0.0,1.0);
 x4[ii]=ran3->Gaus(0.0,1.0);
 x5[ii]=ran4->Gaus(0.0,1.0);
}


//-------------------------------------------------------------------------------------------------------------------------------------
TMatrixD c(dim,dim);
TMatrixD d(dim,dim);
TMatrixD X(dim,dim);
TMatrixD Y(dim,dim);



ifstream fin;
fin.open("Corr_untag_sig.txt");//dim=5;

while(!fin.eof())
{count++;i=count/dim;j=count%dim;
fin>>a[i][j];
if(count==(dim*dim)-1)){break;}
}
fin.close();

a.Print();

TDecompChol b(a);
b.Decompose();
c=b.GetU();//U
//c.Print();
d.Transpose(c);//U^T
d.Print();



TVectorD vec(5);
TVectorD vec2(5);

//vec.Print();
for(int ii=0;ii<num;ii++){

vec[0]= x1[ii];
vec[1]= x2[ii];
vec[2]= x3[ii];
vec[3]= x4[ii];
vec[4]= x5[ii];

vec2=vec;
vec2 *= d;

y1[ii] = vec2[0];
y2[ii] = vec2[1];
y3[ii] = vec2[2];
y4[ii] = vec2[3];
y5[ii] = vec2[4];

}

/*
cout<<"----------------------------------------------------------------------"<<endl;
cout<<corr(x1,x1)<<"\t"<<corr(x1,x2)<<"\t"<<corr(x1,x3)<<"\t"<<corr(x1,x4)<<"\t"<<corr(x1,x5)<<endl;
cout<<corr(x2,x1)<<"\t"<<corr(x2,x2)<<"\t"<<corr(x2,x3)<<"\t"<<corr(x2,x4)<<"\t"<<corr(x2,x5)<<endl;
cout<<corr(x3,x1)<<"\t"<<corr(x3,x2)<<"\t"<<corr(x3,x3)<<"\t"<<corr(x3,x4)<<"\t"<<corr(x3,x5)<<endl;
cout<<corr(x4,x1)<<"\t"<<corr(x4,x2)<<"\t"<<corr(x4,x3)<<"\t"<<corr(x4,x4)<<"\t"<<corr(x4,x5)<<endl;
cout<<corr(x5,x1)<<"\t"<<corr(x5,x2)<<"\t"<<corr(x5,x3)<<"\t"<<corr(x5,x4)<<"\t"<<corr(x5,x5)<<endl;

cout<<"----------------------------------------------------------------------"<<endl;

cout<<corr(y1,y1)<<"\t"<<corr(y1,y2)<<"\t"<<corr(y1,y3)<<"\t"<<corr(y1,y4)<<"\t"<<corr(y1,y5)<<endl;
cout<<corr(y2,y1)<<"\t"<<corr(y2,y2)<<"\t"<<corr(y2,y3)<<"\t"<<corr(y2,y4)<<"\t"<<corr(y2,y5)<<endl;
cout<<corr(y3,y1)<<"\t"<<corr(y3,y2)<<"\t"<<corr(y3,y3)<<"\t"<<corr(y3,y4)<<"\t"<<corr(y3,y5)<<endl;
cout<<corr(y4,y1)<<"\t"<<corr(y4,y2)<<"\t"<<corr(y4,y3)<<"\t"<<corr(y4,y4)<<"\t"<<corr(y4,y5)<<endl;
cout<<corr(y5,y1)<<"\t"<<corr(y5,y2)<<"\t"<<corr(y5,y3)<<"\t"<<corr(y5,y4)<<"\t"<<corr(y5,y5)<<endl;
cout<<"----------------------------------------------------------------------"<<endl;
*/

X[0][0]=corr(x1,x1);	X[0][1]=corr(x1,x2);	X[0][2]=corr(x1,x3);	X[0][3]=corr(x1,x4);	X[0][4]=corr(x1,x5);
X[1][0]=corr(x2,x1);	X[1][1]=corr(x2,x2);	X[1][2]=corr(x2,x3);	X[1][3]=corr(x2,x4);	X[1][4]=corr(x2,x5);
X[2][0]=corr(x3,x1);	X[2][1]=corr(x3,x2);	X[2][2]=corr(x3,x3);	X[2][3]=corr(x3,x4);	X[2][4]=corr(x3,x5);
X[3][0]=corr(x4,x1);	X[3][1]=corr(x4,x2);	X[3][2]=corr(x4,x3);	X[3][3]=corr(x4,x4);	X[3][4]=corr(x4,x5);
X[4][0]=corr(x5,x1);	X[4][1]=corr(x5,x2);	X[4][2]=corr(x5,x3);	X[4][3]=corr(x5,x4);	X[4][4]=corr(x5,x5);


Y[0][0]=corr(y1,y1);	Y[0][1]=corr(y1,y2);	Y[0][2]=corr(y1,y3);	Y[0][3]=corr(y1,y4);	Y[0][4]=corr(y1,y5);
Y[1][0]=corr(y2,y1);	Y[1][1]=corr(y2,y2);	Y[1][2]=corr(y2,y3);	Y[1][3]=corr(y2,y4);	Y[1][4]=corr(y2,y5);
Y[2][0]=corr(y3,y1);	Y[2][1]=corr(y3,y2);	Y[2][2]=corr(y3,y3);	Y[2][3]=corr(y3,y4);	Y[2][4]=corr(y3,y5);
Y[3][0]=corr(y4,y1);	Y[3][1]=corr(y4,y2);	Y[3][2]=corr(y4,y3);	Y[3][3]=corr(y4,y4);	Y[3][4]=corr(y4,y5);
Y[4][0]=corr(y5,y1);	Y[4][1]=corr(y5,y2);	Y[4][2]=corr(y5,y3);	Y[4][3]=corr(y5,y4);	Y[4][4]=corr(y5,y5);

X.Print();
Y.Print();



}//void corr_test(void)



double mean(double x[num])
{
double sum=0;
for(int i=0;i<num;i++)
{sum=sum+x[i];}
sum=sum/num;
return sum;
}


double stdev(double x[num])
{
double sum=0; double mn=mean(x);
for(int i=0;i<num;i++)
{sum=sum+pow((x[i]-mn),2);}
sum=sum/num;
sum=sqrt(sum);
return sum;
}

double cov(double x1[num],double x2[num])
{
double sum=0; double mn1=mean(x1); double mn2=mean(x2);
for(int i=0;i<num;i++)
{sum=sum+((x1[i]-mn1)*(x2[i]-mn2));}
sum=sum/num;
return sum;
}

double corr(double x1[num],double x2[num])
{
double cv=cov(x1,x2);
double st1=stdev(x1);
double st2=stdev(x2);
double cr = cv/(st1*st2);
return cr;
}






