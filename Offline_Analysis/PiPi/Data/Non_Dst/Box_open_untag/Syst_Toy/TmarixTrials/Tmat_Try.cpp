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
#include "TVectorT.h"
#include "TMatrixD.h"



using namespace RooFit;
void Tmat_Try()
{

TMatrixD a(3,3);
TMatrixD c(3,3);
TMatrixD d(3,3);
TMatrixD e(3,3);

a[0][0]=1.0;
a[0][1]=0.5;
a[0][2]=-0.5;
a[1][0]=0.5;
a[1][1]=1.0;
a[1][2]=0.5;
a[2][0]=-0.5;
a[2][1]=0.5;
a[2][2]=1.0;

a.Print();
//TDecompChol b(3);
TDecompChol b(a);
b.Decompose();
b.Print();
c=b.GetU();
c.Print();//U
d.Transpose(c);//U^T
e=d*c;

e.Print();


}









