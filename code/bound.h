#ifndef _BOUND_H_
#define _BOUND_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>

#include "TSystem.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/GaussIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"


/************* Constants *****************/
const double GeVfm = 1.0 / 0.1973269718;//GeV * fm
const double RadToDeg = 180.0 / M_PI;
const double DegToRad = M_PI / 180.0;
const double Mp = 0.9382720;
const double Mn = 0.9395654;
const double MN = (Mp + Mn) / 2.0;
const double Mphi = 1.019455;
/************* End of Constants **********/

/************* Parameters ****************/
double NA;//nucleon number
double Lambda = 3.0 * GeVfm;//cut off parameter
double Md;//bound state mass
/************* End of Parameters *********/

//Wave function and effective potential
const int _MAX1_ = 166;
double _x1[_MAX1_], _y1[_MAX1_];
const int _MAX2_ = 40;
double _x2[_MAX2_], _y2[_MAX2_];
const int _MAX3_ = 301;
double _x3[_MAX3_], _y3[_MAX3_];

ROOT::Math::Interpolator ip1(_MAX1_, ROOT::Math::Interpolation::kCSPLINE);
ROOT::Math::Interpolator ip2(_MAX2_, ROOT::Math::Interpolation::kCSPLINE);
ROOT::Math::Interpolator ip3(_MAX3_, ROOT::Math::Interpolation::kCSPLINE);

int readgrid(){
  double tmp;
  //wave function
  std::ifstream f1("wf.dat");
  if (!f1.is_open()){
    std::cerr << "File wf1.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX1_; i++)
    f1 >> _x1[i] >> tmp >> _y1[i];
  f1.close();
  ip1.SetData(_MAX1_, _x1, _y1);
  //potential
  std::ifstream f2("Veff.dat");
  if (!f2.is_open()){
    std::cerr << "File Veff.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX2_; i++)
    f2 >> _x2[i] >> _y2[i] >> tmp >> tmp;
  f2.close();
  ip2.SetData(_MAX2_, _x2, _y2);
  return 0;
}

//wave function u(r) with normalization int u(r)^2 dr = 1
double _Nur_ = 1.0;//wave function renormalization
double ur(double r, void * par = 0){//r in unit GeV^-1
  double rfm = r / GeVfm;//trans to fm unit
  double result = 0.0;
  if (rfm < 8.20)
    result = ip1.Eval(rfm);
  result = _Nur_ * result;//renormalize
  return result;//in unit GeV^1/2
}

double ur2(double r, void * par = 0){//r in unit GeV^-1
  double result = pow(ur(r), 2);
  return result;//in unit GeV
}

double urone(){//check ur normalization
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ig.SetFunction(&ur2);
  ig.SetRelTolerance(0.001);
  double su = ig.Integral(0.0, 18.20 * GeVfm);
  _Nur_ = 1.0 / sqrt(su);
  //std::cout << _Nur_ << std::endl;
  return ig.Integral(0.0, 8.5 * GeVfm);
}

//Effective potential
double Veff(double r, void * par = 0){//r in unit GeV^-1
  double rfm = r / GeVfm;//trans to fm unit
  double result = 0.0;
  if (rfm < 2.0)
    result = ip2.Eval(rfm);
  result = result /1000.0;//trans to GeV unit
  return result;//in unit GeV
}

//FQ integrand
double FQint(double r, void * pQ = 0){//r in unit GeV^-1
  double Q = * ((double *) pQ);//Q in unit GeV
  double result;
  if (Q == 0)
    result = r * Veff(r) * ur(r) / pow(2.0 * M_PI, 1.5);
  else
    result = sin(Q * r) / Q * Veff(r) * ur(r) * exp(-Lambda*Lambda*Q*Q) / pow(2.0 * M_PI, 1.5);
  return result;//in unit GeV^1/2
}

//Generate FQ.dat and interpolate
int makeFQ(){//Q in unit GeV;
  FILE * fp;
  fp = fopen("FQ.dat","w");
  double Q;
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ig.SetFunction(&FQint, &Q);
  ig.SetRelTolerance(0.001);
  for (Q = 0.0; Q < 0.3; Q+=0.001){
    fprintf(fp, "%.4E  %.8E\n", Q, ig.Integral(0.0, 2.5 * GeVfm));
  }
  fprintf(fp, "%.4E  %.8E\n", 0.3, 0.0);
  fclose(fp);
  std::ifstream f3("FQ.dat");
  for (int i = 0; i < _MAX3_; i++)
    f3 >> _x3[i] >> _y3[i];
  f3.close();
  ip3.SetData(_MAX3_, _x3, _y3);
  return 0;
}

//phi N -> d transition amplitude
double FQ(double Q, void * par = 0){//Q in unit GeV
  double result = 0.0;
  if (Q < 0.3)
    result = ip3.Eval(Q);
  return result;//in unit GeV^-1/2
}



#endif
