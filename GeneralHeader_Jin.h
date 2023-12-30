//check if x not nan or inf: isfininte(x);
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TF2.h>
#include <TH3.h>
#include <TF3.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TMath.h>
#include <TFile.h>
#include <TLine.h>
#include <vector>
#include <TString.h>
#include <TVector.h>
#include <TVectorT.h>
#include <TFormula.h>
#include <TClass.h>
#include <Riostream.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TDirectory.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <cmath>
#include <math.h>
#include <sys/stat.h>
#include <TAxis.h>
#include <TLine.h>
#include "TLegend.h"
#include "TBox.h"
#include "TEllipse.h"
#define sqrt2 1.41421356237
#define pi 3.14159265358979323846
//  double cont_sigma[5]={2.295749,6.180074,11.829158,19.333909,28.743702};  // delta_chi2 values for 1 - 5 sigma in 2-dimensional (2 degrees of freedom) space.
//  const int ncont = 5;
using namespace std;
using namespace TMath;
bool fexist(const char* file1){
 ifstream infile(file1);
 return infile.good();
}
int fsize(const char* file1){
 int size=0, exist;
 exist = fexist(file1);
 if(exist!=1) {
  printf("\n File %s doesn't exist\n",file1);
  return 0;
 }
 streampos begin, end;
 ifstream myfile(file1, ios::binary);
 begin = myfile.tellg();
 myfile.seekg (0, ios::end);
 end = myfile.tellg();
 myfile.close();
 size = end - begin;
  return size;
}
double lin_int(double xx,double x1,double x2,double y1,double y2){ // linear interpolation
  double m, y;
  if(xx==x1) return y1;
  if(xx==x2) return y2;
  if(x1 == x2 && y1 == y2) return y1;
  else if(x1==x2 && y1!=y2){ printf("\n Bad input to lin_int in my_h.h: x1 = x2, y1 != y2. Returning y1"); return y1;}
  m = (y2-y1)/(x2-x1);
  y = m*(xx-x1) + y1;
  return y;
}

/*
//Make spectrum prediction tgraph TGraphErrors(n,x,y,ex,ey);
TGraphErrors *make_spec_graph(char hname[], char htitle[], double *dda){ //dda will have the y-values of the spectral
 const int numb_ent = 24;
 double xerr[numb_ent],yerr[numb_ent],xx[numb_ent], melec=0.511; //electron mass
 for(int jj=0; jj<21; jj++){
        xx[jj] = 3.75 + jj*0.5 - melec; //4 to 14 MeV, then convert to kinetic energy
 }
 xx[21] = 14.5 - melec;
 xx[22] = 15.5 - melec;
 xx[23] = 18.0 - melec;
 for(int jj=0; jj<24; jj++){
        yerr[jj]=0.;
        xerr[jj] = 0.25;
        if(jj==23) xerr[jj]=2.0;  //16-20
        if(jj==22 || jj==21) xerr[jj]=0.5; //14-15 & 15-16
 }
 TGraphErrors *tge_tmp = new TGraphErrors(numb_ent,xx,dda,xerr,yerr);
 tge_tmp->SetName(hname);
 tge_tmp->SetTitle(htitle);
 return tge_tmp;
}
*/
