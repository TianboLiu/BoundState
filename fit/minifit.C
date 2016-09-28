#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

const int NPOINT = 36;
double rs[NPOINT], ds[NPOINT], err[NPOINT];
const double Mp = 0.938272046;
const double Mphi = 1.019455;


double fitfunc(const double * A){
  double qc, Q, sigma;
  double sum = 0.0;
  for (int i = 0; i < NPOINT; i++){
    qc = (rs[i]*rs[i] - Mp*Mp) / (2.0 * rs[i]);
    Q = sqrt(pow(rs[i]*rs[i] - Mphi*Mphi - Mp*Mp, 2) - 4.0 * Mphi*Mphi * Mp*Mp) / (2.0 * rs[i]);
    sigma = 389.379 * 2.0 * M_PI * A[0] * A[0] * M_PI * sqrt(Mp*Mp + qc*qc) / (sqrt(Mp*Mp + qc*qc) + qc) * M_PI * Q * sqrt(Mphi*Mphi + Q*Q) * sqrt(Mp*Mp + Q*Q) / (sqrt(Mphi*Mphi + Q*Q) + sqrt(Mp*Mp + Q*Q));
    cout << sigma << endl;
    sum += pow((sigma - ds[i])/err[i], 2);
  }
  return sum;
}

double fitfunc2(const double * A){
  double qc, Q, sigma;
  double sum = 0.0;
  for (int i = 0; i < NPOINT; i++){
    qc = (rs[i]*rs[i] - Mp*Mp) / (2.0 * rs[i]);
    Q = sqrt(pow(rs[i]*rs[i] - Mphi*Mphi - Mp*Mp, 2) - 4.0 * Mphi*Mphi * Mp*Mp) / (2.0 * rs[i]);
    sigma = 389.379 * A[0] * A[0] / (32.0 * M_PI) * Q / (qc * (sqrt(Mp*Mp + qc*qc) + qc)) / (sqrt(Mphi*Mphi + Q*Q) + sqrt(Mp*Mp + Q*Q));
    cout << sigma << endl;
    sum += pow((sigma - ds[i])/err[i], 2);
  }
  return sum;
}

int getdata(){
  ifstream infile("set0.dat");
  double tmp;
  for (int i = 0; i < NPOINT; i++){
    infile >> rs[i] >> tmp >> ds[i] >> err[i];
  }
  infile.close();
  return 0;
}

int minimizer(const char * minName = "Minuit", const char * algoName = ""){
  getdata();
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&fitfunc, 1);
  min->SetFunction(f);
  min->SetVariable(0, "A", 0.008, 1e-4);
  min->Minimize();
  cout << min->MinValue() << endl;
  return 0;
}

#ifndef __CINT__

int main(){
  return minimizer();
}

#endif
