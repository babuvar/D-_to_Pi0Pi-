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
void transpose()
{

TMatrixD a(2,2);
TMatrixD b(2,2);
a[0][0]=1.0;
a[0][1]=1.0;
a[1][0]=2.0;
a[1][1]=2.0;

b[0][0]=1.0;
b[0][1]=1.0;
b[1][0]=3.0;
b[1][1]=3.0;

b.Transpose(a);

a.Print();
b.Print();

}

