#ifndef _BOUND_H_
#define _BOUND_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>

#include "TSystem.h"
#include "Math/Functor.h"
#include "Math/ParamFunctor.h"
#include "Math/Factory.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/GaussIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/GSLMCIntegrator.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph2D.h"


/************* Constants *****************/
const double GeVfm = 1.0 / 0.1973269718;//GeV * fm
const double RadToDeg = 180.0 / M_PI;
const double DegToRad = M_PI / 180.0;
const double Mp = 0.9382720;
const double Mn = 0.9395654;
const double MN = (Mp + Mn) / 2.0;
const double Mphi = 1.019455;
const double Mkaon = 0.493677;
const double Mpion = 0.13957018;
const double MLambda = 1.115683;
/************* End of Constants **********/

/************* Parameters ****************/
double NA = 12.0;//nucleon number
double Lambda = 3.0 * GeVfm;//cut off parameter
double Md = MN + Mphi - 0.0448;//bound state mass
double Ebeam = 1.45;//photon beam energy
double b = 1.64 * GeVfm;//C12 shell model harmonic oscillator length
/************* End of Parameters *********/

//Wave function and effective potential
const int _MAX1_ = 161;
double _x1[_MAX1_], _y1[_MAX1_];
const int _MAX2_ = 60;
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
    f1 >> _x1[i] >> _y1[i];
  f1.close();
  ip1.SetData(_MAX1_, _x1, _y1);
  //potential
  std::ifstream f2("Veff.dat");
  if (!f2.is_open()){
    std::cerr << "File Veff.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX2_; i++)
    f2 >> _x2[i] >> _y2[i];
  f2.close();
  ip2.SetData(_MAX2_, _x2, _y2);
  return 0;
}

//wave function u(r) with normalization int u(r)^2 dr = 1
double _Nur_ = 1.0;//wave function renormalization
double ur(double r, void * par = 0){//r in unit GeV^-1
  double rfm = r / GeVfm;//trans to fm unit
  double result = 0.0;
  if (rfm < 8.0)
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
  double su = ig.Integral(0.0, 18.20 * GeVfm) / pow(_Nur_, 2);
  _Nur_ = 1.0 / sqrt(su);
  //std::cout << _Nur_ << std::endl;
  return ig.Integral(0.0, 8.5 * GeVfm);
}

//Effective potential
double Veff(double r, void * par = 0){//r in unit GeV^-1
  double rfm = r / GeVfm;//trans to fm unit
  double result = 0.0;
  if (rfm < 5.0)
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
  ig.SetAbsTolerance(0.0);
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

int LoadFQ(){//load FQ grid from existing file
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

//gamma N -> phi N transition amplitude
double tQ(double Q = 0, void * par = 0){
  return 0.005;//in unit GeV^-2
}

//C12 nuclear wave function
double wfC12(double p, void * par = 0){//p in unit GeV
  //Normalization: int wfC12^2 p^2 sin(theta) dp dtheta dphi = 1
  double Normalization = sqrt(8.0 * pow(b, 3.0) / 3.0 / sqrt(M_PI));
  //double Normalization = sqrt(4.0 * pow(b, 3.0) / sqrt(M_PI));
  double result = 1.0 / sqrt(4.0 * M_PI) * Normalization * b * p * exp(-b*b*p*p/2.0);
  //double result = 1.0 / sqrt(4.0 * M_PI) * Normalization * exp(-b*b*p*p/2.0);
  return result;
}

//Total transition amplitude 
double Tk(TLorentzVector pf[3]){//4-momentum of k, p, pd in unit GeV
  TLorentzVector k = pf[0];//4-momentum of intermediate phi
  TLorentzVector p = pf[1];//4-momentum of final nucleon
  TLorentzVector pd = pf[2];//4-momentum of bound state d
  TLorentzVector q(0.0, 0.0, Ebeam, Ebeam);//4-momentum of photon beam
  TLorentzVector p1 = p + k - q;//4-momentum of struck nucleon
  p1.SetE(sqrt(pow(p1.P(), 2) + MN*MN));//set energy on-shell
  TLorentzVector p2 = pd - k;//4-momentum of reaction nucleon
  p2.SetE(sqrt(pow(p2.P(), 2) + MN*MN));//set energy on-shell
  TLorentzVector kp2 = k + p2;//4-momentum of phi-N system
  double s = kp2.M2();//invariant mass square of phi-N system
  double Q2 = (pow(s-Mphi*Mphi-MN*MN, 2) - 4.0*Mphi*Mphi*MN*MN) / (4.0 * s);
  if (Q2 <= 0)
    return 0;
  double Q = sqrt(Q2);
  double Fvalue = FQ(Q);
  double tvalue = tQ(Q);
  double Edenominator = Ebeam + NA * MN - sqrt(pow(p1.P(), 2) + pow((NA-1)*MN, 2)) - k.E() - p.E();//Energy denominator of intermediate state
  double result = wfC12(p1.P()) * wfC12(p2.P()) * Fvalue * tvalue / Edenominator;
  return result;
}
  
double Tint(const double * kk, const double * par){
  //kk: k, theta, phi
  //par (every 3): p, pd
  TLorentzVector pf[3];
  pf[0].SetXYZM(kk[0] * sin(kk[1]) * cos(kk[2]), kk[0] * sin(kk[1]) * sin(kk[2]), kk[0] * cos(kk[1]), Mphi);//set 4-momentum of intermediate phi meson
  pf[1].SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), MN);//set 4-momentum of final nucleon
  pf[2].SetXYZM(par[3] * sin(par[4]) * cos(par[5]), par[3] * sin(par[4]) * sin(par[5]), par[3] * cos(par[4]), Md);//set 4-momentum of final bound state
  double result = Tk(pf) * kk[0] * kk[0] * sin(kk[1]);
  if (isnan(result)) return 0;
  return result;
}

double Tfi(TLorentzVector p, TLorentzVector pd){
  double par[6] = {p.P(), p.Theta(), p.Phi(), pd.P(), pd.Theta(), pd.Phi()};
  double xl[3] = {0, 0, -M_PI};//integration lower bound
  double xu[3] = {2.0, M_PI, M_PI};//integration upper bound
  ROOT::Math::WrappedParamFunction<> wf(&Tint, 3, 6);
  wf.SetParameters(par);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 100000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}

double T2(const double * omega, const double * par){
  TLorentzVector pd;
  pd.SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), Md);//set 4-momentum of final bound state
  double pp2 = pow(Ebeam + 2.0 * MN - pd.E(), 2) - MN*MN;
  if (pp2 <= 0.0)
    return 0;
  double pp = sqrt(pp2);//momentum of final nucleon constrained by energy conservation
  TLorentzVector p;
  p.SetXYZM(pp * sin(omega[0]) * cos(omega[1]), pp * sin(omega[0]) * sin(omega[1]), pp * cos(omega[0]), MN);//set 4-momentum of final nucleon
  double result = pow(Tfi(p, pd), 2) * sin(omega[0]) * p.P() * p.E();
  if (isnan(result)) return 0;
  return result;
}

double dsigma(const double * Pd){//ds / dpd dcostheta
  TLorentzVector pd;
  pd.SetXYZM(Pd[0] * sin(Pd[1]), 0.0, Pd[0] * cos(Pd[1]), Md);
  double par[3] = {pd.P(), pd.Theta(), pd.Phi()};
  double xl[2] = {0.0, -M_PI};//integration lower boundary
  double xu[2] = {M_PI, M_PI};//integration upper boundary
  ROOT::Math::WrappedParamFunction<> wf(&T2, 2, 3);
  wf.SetParameters(par);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.01, 500);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  if (isnan(result)) return 0;
  return result * pow(2.0 * M_PI, 5) * pow(pd.P(), 2);
}

double sigmaint(const double * Pd){
  double result = dsigma(Pd) * sin(Pd[1]);
  std::cout << result << std::endl;
  return result;
}

double T2MC(const double * omega, void * vpar){//
  double * par = (double *) vpar;
  return T2(omega, par);
}

double dsigmaMC(const double * Pd, double * err = 0){
  TLorentzVector pd;
  pd.SetXYZM(Pd[0] * sin(Pd[1]), 0.0, Pd[0] * cos(Pd[1]), Md);
  double par[3] = {pd.P(), pd.Theta(), pd.Phi()};
  double xl[2] = {0.0, -M_PI};//integration lower boundary
  double xu[2] = {M_PI, M_PI};//integration upper boundary
  ROOT::Math::WrappedParamFunction<> wf(&T2, 2, 3);
  wf.SetParameters(par);
  ROOT::Math::GSLMCIntegrator ig(ROOT::Math::IntegrationMultiDim::kVEGAS, 1.0e-20, 0.02, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  if (err != 0) err[0] = ig.Error()/result;
  return result * pow(2.0 * M_PI, 5) * pow(pd.P(), 2);
}

double sigma(){//Total cross section
  double xl[2] = {0.0, 0.0};//integration lower boundary
  double xu[2] = {1.5, M_PI};//integration upper boundary
  ROOT::Math::Functor wf(&sigmaint, 2);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
  ig.SetFunction(wf);
  ig.SetAbsTolerance(0.0);
  ig.SetRelTolerance(0.01);
  double result = ig.Integral(xl, xu);
  return result;//in unit GeV^-2
}

//generate 2D table of dsigma for interpolation
int makeDS(const char * filename = "ds.dat"){//pd in GeV, costheta, dsigma / dpd dcostheta
  FILE * fp;
  fp = fopen(filename, "w");
  double Pd[2];
  double ds;
  std::cout << "Generating dsigma grid ..." << std::endl;
  for (Pd[0] = 0.0; Pd[0] < 1.4; Pd[0] += 0.02){
    for (Pd[1] = 0.0; Pd[1] <= 0.501 * M_PI; Pd[1] += 0.01 * M_PI){
      ds = dsigma(Pd);
      std::cout << Pd[0] << "   " << Pd[1] << "  " << cos(Pd[1]) << "   " << ds << std::endl;
      fprintf(fp, "%.6E  %.6E  %.6E\n", Pd[0], cos(Pd[1]), ds);
    }
  }
  fclose(fp);
  return 0;
}

//decay events generator
int decay0(TLorentzVector pd, TLorentzVector * pfs){//NKK decay channel
  double masses[3] = {Mp, Mkaon, Mkaon};
  TGenPhaseSpace event;
  event.SetDecay(pd, 3, masses);
  event.Generate();
  TLorentzVector * p;
  p = event.GetDecay(0);
  pfs[0].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//p
  p = event.GetDecay(1);
  pfs[1].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K+
  p = event.GetDecay(2);
  pfs[2].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K-
  return 0;
}

int decay1(TLorentzVector pd, TLorentzVector * pfs){//LambdaK decay channel
  TLorentzVector * p;
  double mass1[2] = {MLambda, Mkaon};
  TGenPhaseSpace step1;
  step1.SetDecay(pd, 2, mass1);
  step1.Generate();
  p = step1.GetDecay(1);
  pfs[2].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K+
  p = step1.GetDecay(0);
  TLorentzVector pL(p->Px(), p->Py(), p->Pz(), p->E());//Lambda
  double mass2[2] = {Mp, Mpion};
  TGenPhaseSpace step2;
  step2.SetDecay(pL, 2, mass2);
  step2.Generate();
  p = step2.GetDecay(0);
  pfs[0].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//p
  p = step2.GetDecay(1);
  pfs[1].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//pi-
  return 0;
}

//grid for dsigma interpolation
TGraph2D ds2D(1);
int LoadDS(const char * datfile){
  ifstream f0(datfile);
  if (!f0.is_open()){
    std::cerr << "ds file does not exist!" << std::endl;
    return 1;
  }
  int nline = 0;
  std::string tmp;
  while (std::getline(f0, tmp))
    nline++;
  f0.close();
  double x, y, z;  
  ifstream f1(datfile);
  ds2D.Set(nline);
  for (int i = 0; i < nline; i++){
    f1 >> x >> y >> z;
    ds2D.SetPoint(i, x, y, z);
  }
  f1.close();
  return 0;
}

//ds
double ds(const double * pd){//interpolation of ds/dpd dcostheta
  //pd: pd in GeV, costheta
  return ds2D.Interpolate(pd[0], pd[1]);
}
  
double sigmatotal(){//total cross section of bound state production in GeV^-2
  double xl[2] = {0.0, -1.0};//set lower integration boundary
  double xu[2] = {1.4, 1.0};//set upper integration boundary
  ROOT::Math::Functor wf(&ds, 2);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}



#endif
