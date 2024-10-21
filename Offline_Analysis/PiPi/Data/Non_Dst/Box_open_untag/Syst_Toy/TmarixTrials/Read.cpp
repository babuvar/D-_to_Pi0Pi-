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
#include <fstream>


using namespace RooFit;
void Read()
{
int count=-1;
int i,j;
//int dim=5;
//int dim=2;
int dim=12;
//int dim=6;


ifstream fin;
//fin.open("Corr_untag_sig.txt"); int dim=5;
//fin.open("Corr_untag_bkg.txt"); int dim=2;
fin.open("Corr_tag_sig.txt"); int dim=12;
//fin.open("Corr_tag_bkg.txt"); int dim=6;

TMatrixD a(dim,dim);
TMatrixD c(dim,dim);
TMatrixD d(dim,dim);
TMatrixD e(dim,dim);



while(!fin.eof())
{count++;i=count/dim;j=count%dim;
fin>>a[i][j];
if(count==(dim*dim)-1)){break;}
}
fin.close();

a.Print();

TDecompChol b(a);
b.Decompose();
//b.Print();
c=b.GetU();
//c.Print();//U
d.Transpose(c);//U^T
d.Print();
e=d*c;

//e.Print();


}









