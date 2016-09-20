#ifndef _BOUND_H_
#define _BOUND_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <omp.h>

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
#include "TTree.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TRandom3.h"


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
double Md = MN + Mphi - 0.0077;//bound state mass
double Ebeam = 1.45;//photon beam energy
double b = 1.64 * GeVfm;//C12 shell model harmonic oscillator length
int L1 = 1;
int L2 = 1; 
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

int readgrid(const char * fwf = "wf.dat", const char * fVeff = "Veff.dat"){
  //double tmp;
  //wave function
  std::ifstream f1(fwf);
  if (!f1.is_open()){
    std::cerr << "File wf1.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX1_; i++)
    f1 >> _x1[i] >> _y1[i];
  f1.close();
  ip1.SetData(_MAX1_, _x1, _y1);
  //potential
  std::ifstream f2(fVeff);
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
int makeFQ(const char * fFQ = "FQ.dat"){//Q in unit GeV;
  FILE * fp;
  fp = fopen(fFQ,"w");
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
  std::ifstream f3(fFQ);
  for (int i = 0; i < _MAX3_; i++)
    f3 >> _x3[i] >> _y3[i];
  f3.close();
  ip3.SetData(_MAX3_, _x3, _y3);
  return 0;
}

int LoadFQ(const char * fFQ = "FQ.dat"){//load FQ grid from existing file
  std::ifstream f3(fFQ);
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
double wfC12s(double p, void * par = 0){//p in unit GeV
  //Normalization: int wfC12^2 p^2 sin(theta) dp dtheta dphi = 1
  double Normalization = sqrt(4.0 * pow(b, 3.0) / sqrt(M_PI));
  double result = 1.0 / sqrt(4.0 * M_PI) * Normalization * exp(-b*b*p*p/2.0);
  return result;
}

double wfC12p(double p, void * par = 0){//p in unit GeV
  //Normalization: int wfC12^2 p^2 sin(theta) dp dtheta dphi = 1
  double Normalization = sqrt(8.0 * pow(b, 3.0) / 3.0 / sqrt(M_PI));
  double result = 1.0 / sqrt(4.0 * M_PI) * Normalization * b * p * exp(-b*b*p*p/2.0);
  return result;
}

double wfC12(double p, int L){
  if (L == 0)
    return wfC12s(p);
  else if (L == 1)
    return wfC12p(p);
  else
    return 0;
}

//Total transition amplitude 
double Tk(const double * kk, const TLorentzVector p, const TLorentzVector pd){//4-momentum of k, p, pd in unit GeV
  //kk: theta_k, phi_k
  TLorentzVector q(0.0, 0.0, Ebeam, Ebeam);//4-momentum of photon beam
  double Es = Ebeam + NA * MN - p.E();//Remaining Energy for phi and A-1
  if (Es < (NA - 1) * MN + Mphi) return 0;//below phi production threshold
  TLorentzVector l = p - q;//4-momentum p1 - k
  double CA0 = cos(kk[0]) * cos(l.Theta()) + sin(kk[0]) * sin(l.Theta()) * cos(kk[1] - l.Phi());
  double bb = pow((NA - 1) * MN, 2) + pow(l.P(), 2);
  double dcm = bb*bb - 2.0*bb*(Es*Es+Mphi*Mphi) + pow(Es, 4) - 2.0*Es*Es*Mphi*Mphi + 4.0*CA0*CA0*pow(l.P(), 2)*Mphi*Mphi + pow(Mphi, 4);
  if (dcm <= 0) return 0;//no match k value
  double k1 = (l.P() * CA0 * (bb - Es*Es - Mphi*Mphi) + Es * sqrt(dcm)) / (2.0 * (Es*Es - CA0*CA0*l.P()*l.P()));
  double k2 = (l.P() * CA0 * (bb - Es*Es - Mphi*Mphi) - Es * sqrt(dcm)) / (2.0 * (Es*Es - CA0*CA0*l.P()*l.P()));
  //std::cout << k1 << " " << k2 << std::endl; 
  double result1 = 0.0;
  double result2 = 0.0;
  if (k1 > 0){
    double res = 1.0 / (- (k1 + l.P() * CA0) / sqrt(k1*k1 + 2.0*k1*l.P()*CA0 + bb) - k1 / sqrt(k1*k1 + Mphi*Mphi));
    TLorentzVector k;
    k.SetXYZM(k1 * sin(kk[0]) * cos(kk[1]), k1 * sin(kk[0]) * sin(kk[1]), k1 * cos(kk[0]), Mphi);
    TLorentzVector p1 = p + k - q;//4-momentum of struck nucleon
    TLorentzVector p2 = pd - k;//4-momentum of reaction nucleon
    p2.SetE(sqrt(p2.P()*p2.P() + MN*MN));//set p2 on-shell energy
    TLorentzVector kp2 = k + p2;//4-momentum of phi-N system
    double s = kp2.M2();//invariant mass square of phi-N system
    double Q2 = (pow(s-Mphi*Mphi-MN*MN, 2) - 4.0*Mphi*Mphi*MN*MN) / (4.0 * s);
    if (Q2 > 0){
      double Q = sqrt(Q2);
      double Fvalue = FQ(Q);
      //std::cout << Q << std::endl;
      double tvalue = tQ(Q);
      double wave = wfC12(p1.P(), L1) * wfC12(p2.P(), L2);
      result1 = M_PI * wave * Fvalue * tvalue * res * k1 * k1;
    }
  }
  if (k2 > 0){
    double res = 1.0 / (- (k2 + l.P() * CA0) / sqrt(k2*k2 + 2.0*k2*l.P()*CA0 + bb) - k2 / sqrt(k2*k2 + Mphi*Mphi));
    TLorentzVector k;
    k.SetXYZM(k2 * sin(kk[0]) * cos(kk[1]), k2 * sin(kk[0]) * sin(kk[1]), k2 * cos(kk[0]), Mphi);
    TLorentzVector p1 = p + k - q;//4-momentum of struck nucleon
    TLorentzVector p2 = pd - k;//4-momentum of reaction nucleon
    p2.SetE(sqrt(p2.P()*p2.P() + MN*MN));//set p2 on-shell energy
    TLorentzVector kp2 = k + p2;//4-momentum of phi-N system
    double s = kp2.M2();//invariant mass square of phi-N system
    double Q2 = (pow(s-Mphi*Mphi-MN*MN, 2) - 4.0*Mphi*Mphi*MN*MN) / (4.0 * s);
    if (Q2 > 0){
      double Q = sqrt(Q2);
      double Fvalue = FQ(Q);
      double tvalue = tQ(Q);
      double wave = wfC12(p1.P(), L1) * wfC12(p2.P(), L2);
      result2 = M_PI * wave * Fvalue * tvalue * res * k2 * k2;
    }
  }
  return result1 + result2;
}
  
double Tint(const double * kk, const double * par){//kk: theta, phi
  //par (every 3): p, pd
  TLorentzVector p, pd;
  p.SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), MN);
  pd.SetXYZM(par[3] * sin(par[4]) * cos(par[5]), par[3] * sin(par[4]) * sin(par[5]), par[3] * cos(par[4]), Md);
  double result = Tk(kk, p, pd) * sin(kk[0]);
  return result;
}

double Tfi(TLorentzVector p, TLorentzVector pd){
  double par[6] = {p.P(), p.Theta(), p.Phi(), pd.P(), pd.Theta(), pd.Phi()};
  double xl[3] = {0, -M_PI};//integration lower bound
  double xu[3] = {M_PI, M_PI};//integration upper bound
  ROOT::Math::WrappedParamFunction<> wf(&Tint, 2, 6);
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
  if (isnan(result)){
    std::cout << "Warning: " << omega[0] << " " << omega[1] << std::endl;
    return 0;
  }
  TLorentzVector q(0, 0, Ebeam, Ebeam);
  TLorentzVector l = q- p - pd;
  //return result * exp(-pow(l.P()/(NA-2)/MN, 2));
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
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  if (isnan(result)){
    std::cout << "Warning: " << Pd[0] << " " << Pd[1] << std::endl;
    return 0;
  }
  return result * pow(2.0 * M_PI, 5) * pow(pd.P(), 2);
}

double dsigmaT(const double * Pd){//combine ss sp ps pp
  double ss, sp, ps, pp;
  L1 = 0; L2 = 0;
  ss = dsigma(Pd);
  L1 = 0; L2 = 1;
  sp = dsigma(Pd);
  L1 = 1; L2 = 0;
  ps = dsigma(Pd);
  L1 = 1; L2 = 1;
  pp = dsigma(Pd);
  return 12.0*ss + 32.0*sp + 32.0*ps + 56.0*pp;
}

double sigmaint(const double * Pd){
  double result = dsigma(Pd) * sin(Pd[1]);
  if (result < 0.0){
    //std::cout << "Warning" << Pd[0] << " " << Pd[1] << " " << result << std::endl;
    result = 0;
  }
  return result;
}

double sigma(){//Total cross section
  double xl[2] = {0.0, 0.0};//integration lower boundary
  double xu[2] = {1.0, M_PI};//integration upper boundary
  ROOT::Math::Functor wf(&sigmaint, 2);
  //ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.01, 1000);
  ROOT::Math::GSLMCIntegrator ig(ROOT::Math::IntegrationMultiDim::kVEGAS, 0.0, 0.01, 500);
  ig.SetFunction(wf);
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
  for (Pd[0] = 0.0; Pd[0] < 1.0; Pd[0] += 0.01){
    for (Pd[1] = 0.0; Pd[1] <= M_PI; Pd[1] += 0.005 * M_PI){
      ds = dsigmaT(Pd);
      std::cout  << Ebeam << "---" << Pd[0] << "   " << Pd[1] << "  " << cos(Pd[1]) << "   " << ds << std::endl;
      fprintf(fp, "%.6E  %.6E  %.6E\n", Pd[0], cos(Pd[1]), ds);
    }
  }
  fclose(fp);
  return 0;
}

int makeDS_Parallel(const char * filename = "ds.dat"){//same as makeDS but parallel
  FILE * fp;
  fp = fopen(filename, "w");
  double Pd[100][201][2];
  double ds[100][201];
  std::cout << "Generating dsigma grid Parallel..." << std::endl;
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic) collapse(2)
    for (int i = 0; i < 100; i++){
      for (int j = 0; j <= 200; j++){
	Pd[i][j][0] = 0.01 * i;
	Pd[i][j][1] = 0.005 * M_PI * j;
	ds[i][j] = dsigma(Pd[i][j]);
	std::cout << Ebeam << "---" << Pd[i][j][0] << "   " << Pd[i][j][1] << "  " << cos(Pd[i][j][1]) << "   " << ds[i][j] << std::endl;
      }
    }
  }
  std::cout << "Writing File ..." << std::endl;
  for (int i = 0; i < 100; i++){
    for (int j = 0; j <= 200; j++){
      fprintf(fp, "%.6E  %.6E  %.6E\n", Pd[i][j][0], cos(Pd[i][j][1]), ds[i][j]);
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
int LoadDS(const char * datfile = "ds.dat"){
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
  std::cout << "Loading ds grid ..." << std::endl;
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
  double xl[2] = {0.0, 0.0};//set lower integration boundary
  double xu[2] = {1.0, 1.0};//set upper integration boundary
  ROOT::Math::Functor wf(&ds, 2);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}

//generator
int genData0(const char * datafile, const Long64_t Nsim = 10){
  TFile * fs = new TFile(datafile, "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double pd[2];//pd, theta
  double phi;
  TLorentzVector Pd;//4-momentum of Pd
  TLorentzVector AP[3];//4-momenta of final particle from NKK channel
  TLorentzVector BP[3];//4-momenta of final particle from LambdaK channel
  double xs;//ds / dPd dcostheta in unit GeV^-3
  Ts->Branch("Ebeam", &Ebeam, "Ebeam/D");
  Ts->Branch("ds", &xs, "ds/D");
  Ts->Branch("Pd", "TLorentzVector", &Pd);
  Ts->Branch("AP0", "TLorentzVector", &AP[0]);
  Ts->Branch("AP1", "TLorentzVector", &AP[1]);
  Ts->Branch("AP2", "TLorentzVector", &AP[2]);
  Ts->Branch("BP0", "TLorentzVector", &BP[0]);
  Ts->Branch("BP1", "TLorentzVector", &BP[1]);
  Ts->Branch("BP2", "TLorentzVector", &BP[2]);
  for (Long64_t i = 0; i < Nsim; i++){
    pd[0] = gRandom->Uniform(0.0, 1.4);//Generate pd
    pd[1] = acos(gRandom->Uniform(0.0, 1.0));//Generate costheta
    phi = gRandom->Uniform(-M_PI, M_PI);//Generate phi
    xs = dsigma(pd);//Calculate ds / dpd dcostheta
    Pd.SetXYZM(pd[0] * sin(pd[1]) * cos(phi), pd[0] * sin(pd[1]) * sin(phi), pd[0] * cos(pd[1]), Md);//Set 4-momentum of bound state
    decay0(Pd, AP);//Generate decay from NKK channel
    decay1(Pd, BP);//Generate decay from LambdaK channel
    Ts->Fill();
  }
  fs->Write();
  return 0;
}


int genData(const char * dsfile, const char * datafile, const Long64_t Nsim = 10){
  LoadDS(dsfile);//Load ds grid
  TFile * fs = new TFile(datafile, "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double pd[2];//pd, costheta
  double phi;
  TLorentzVector Pd;//4-momentum of Pd
  TLorentzVector AP[3];//4-momenta of final particle from NKK channel
  TLorentzVector BP[3];//4-momenta of final particle from LambdaK channel
  double xs;//ds / dPd dcostheta
  Ts->Branch("ds", &xs, "ds/D");
  Ts->Branch("Pd", "TLorentzVector", &Pd);
  Ts->Branch("AP0", "TLorentzVector", &AP[0]);
  Ts->Branch("AP1", "TLorentzVector", &AP[1]);
  Ts->Branch("AP2", "TLorentzVector", &AP[2]);
  Ts->Branch("BP0", "TLorentzVector", &BP[0]);
  Ts->Branch("BP1", "TLorentzVector", &BP[1]);
  Ts->Branch("BP2", "TLorentzVector", &BP[2]);
  for (Long64_t i = 0; i < Nsim; i++){
    pd[0] = gRandom->Uniform(0.0, 1.4);//Generate pd
    pd[1] = gRandom->Uniform(0.0, 1.0);//Generate costheta
    phi = gRandom->Uniform(-M_PI, M_PI);//Generate phi
    xs = ds(pd);//Calculate ds / dpd dcostheta
    Pd.SetXYZM(pd[0] * sqrt(1.0 - pd[1]*pd[1]) * cos(phi), pd[0] * sqrt(1.0 - pd[1]*pd[1]) * sin(phi), pd[0] * pd[1], Md);//Set 4-momentum of bound state
    decay0(Pd, AP);//Generate decay from NKK channel
    decay1(Pd, BP);//Generate decay from LambdaK channel
    Ts->Fill();
  }
  fs->Write();
  return 0;
}
  
  



#endif
