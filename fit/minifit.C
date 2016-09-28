#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

double ds[36], err[36];

double fitfunc(const double * A){
  double sum = 0.0;
  for (int i = 0; i < 36; i++){
    sum += pow((A[0] - ds[i])/err[i], 2);
  }
  return sum;
}

int getdata(){
  ifstream infile("set0.dat");
  double tmp;
  for (int i = 0; i < 36; i++){
    infile >> tmp >> tmp >> ds[i] >> err[i];
  }
  infile.close();
  return 0;
}

int minimizer(const char * minName = "Minuit", const char * algoName = ""){
  getdata();
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&fitfunc, 1);
  min->SetFunction(f);
  min->SetVariable(0, "A", 0.05, 1e-3);
  min->Minimize();
  return 0;
}

#ifndef __CINT__

int main(){
  return minimizer();
}

#endif
