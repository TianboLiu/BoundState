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
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TString.h"

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
double Lambda = 2.56 * GeVfm;//cut off parameter
double Md = MN + Mphi - 0.0077;//bound state mass
double Ebeam = 1.45;//photon beam energy
double b = 1.64 * GeVfm;//C12 shell model harmonic oscillator length
int C00 = 12; int C01 = 32; int C10 = 32; int C11 = 56;
//int C00 = 00; int C01 = 00; int C10 = 00; int C11 = 132;
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
double _Nur_ = 26.568;//wave function renormalization
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
  std::cout << _Nur_ << std::endl;
  return ig.Integral(0.0, 8.5 * GeVfm);
}

double rur2(double r, void * par = 0){//r in unit GeV^-1
  double result = r * pow(ur(r), 2);
  return result;//in unit GeV
}

double r0(){//check ur normalization
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ig.SetFunction(&rur2);
  ig.SetRelTolerance(0.001);
  double result = ig.Integral(0.0, 18.20 * GeVfm);
  std::cout << result / GeVfm << std::endl;
  return ig.Integral(0.0, 8.5 * GeVfm);
}

double Pr(double rmax){//r in unit fm
  //Probability of r < rmax
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ig.SetFunction(&ur2);
  ig.SetRelTolerance(0.001);
  return ig.Integral(0.0, rmax * GeVfm);
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
    result = sqrt(4.0 * M_PI) * r * Veff(r) * ur(r) / pow(2.0 * M_PI, 1.5);
  else
    result = sqrt(4.0 * M_PI) * sin(Q * r) / Q * Veff(r) * ur(r) * exp(-Lambda*Lambda*Q*Q) / pow(2.0 * M_PI, 1.5);
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
double tQ(double Eq, double E1, TLorentzVector k, TLorentzVector p){
  const double M0 = 0.594;//invariant amplitude
  double t = M0 / (4.0 * pow(2.0 * M_PI, 3) * sqrt(Eq * E1 * k.E() * p.E()));//in unit GeV^-2
  TLorentzVector kp = k + p;//Total 4-momentum of p and phi
  double theta = k.Angle(kp.Vect());//polar angle of k w.r.t total P
  if (theta == 0) return t;
  double gamma = kp.Gamma();//Get boost gamma
  double beta = kp.Beta();//Get boost beta
  double c = cos(theta);//cos theta w.r.t P
  double ct0 = gamma * (c / sin(theta) - beta * k.E() / k.P() / sin(theta));
  double c0  = ct0 /sqrt(1 + ct0*ct0);//cos theta in c.m.
  double jac = gamma * (1.0 - beta * k.E()/k.P() * c) * pow((1.0-c0*c0)/(1.0 - c*c), 3.0/2.0);
  return jac * t;
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
double Tk(const double * kk, const TLorentzVector p, const TLorentzVector pd, const int L1, const int L2){//4-momentum of k, p, pd in unit GeV
  //kk: theta_k, phi_k
  if (Ebeam + MN - p.E() <= Mphi) return 0;
  double kpole = sqrt(pow(Ebeam + MN - p.E(), 2) - Mphi * Mphi);
  TLorentzVector k;
  k.SetXYZM(kpole * sin(kk[0]) * cos(kk[1]), kpole * sin(kk[0]) * sin(kk[1]), kpole * cos(kk[0]), Mphi);
  TLorentzVector q(0.0, 0.0, Ebeam, Ebeam);//4-momentum of photon beam
  TLorentzVector p1 = k + p - q;//4-momentum of struck nucleon
  TLorentzVector p2 = pd - k;//4-momentum of reaction nucleon
  double Tvalue = tQ(Ebeam, sqrt(MN*MN + pow(p1.P(), 2)), k, p);
  TLorentzVector Q = 1.0 / (MN + Mphi) * (MN * k - Mphi * p2);//relative momentum of phi N system
  double Fvalue = FQ(Q.P());
  //std::cout << Q.P() << " " << p2.P() << " " << k.P() << " " << p1.P() << " " << p.P() << std::endl;
  //std::cout << Tvalue*Tvalue*p.P()*p.E() << std::endl;
  double result = M_PI * wfC12(p1.P(), L1) * wfC12(p2.P(), L2) * Tvalue * Fvalue * k.P() * k.E();
  return result;
}
  
double Tint(const double * kk, const double * par){//kk: theta, phi
  //par (every 3): p, pd
  TLorentzVector p, pd;
  p.SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), MN);
  pd.SetXYZM(par[3] * sin(par[4]) * cos(par[5]), par[3] * sin(par[4]) * sin(par[5]), par[3] * cos(par[4]), Md);
  int L1 = (int) par[6];
  int L2 = (int) par[7];
  double result = Tk(kk, p, pd, L1, L2) * sin(kk[0]);
  return result;
}

double Tfi(TLorentzVector p, TLorentzVector pd, double L1, double L2){
  double par[8] = {p.P(), p.Theta(), p.Phi(), pd.P(), pd.Theta(), pd.Phi(), L1, L2};
  double xl[3] = {0, -M_PI};//integration lower bound
  double xu[3] = {M_PI, M_PI};//integration upper bound
  ROOT::Math::WrappedParamFunction<> wf(&Tint, 2, 8);
  wf.SetParameters(par);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 100000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}

double T2(const double * omega, const double * par){
  TLorentzVector pd;
  pd.SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), Md);//set 4-momentum of final bound state
  double L1 = par[3];
  double L2 = par[4];
  double pp2 = pow(Ebeam + 2.0 * MN - pd.E(), 2) - MN*MN;
  if (pp2 <= 0.0)
    return 0;
  double pp = sqrt(pp2);//momentum of final nucleon constrained by energy conservation
  TLorentzVector p;
  p.SetXYZM(pp * sin(omega[0]) * cos(omega[1]), pp * sin(omega[0]) * sin(omega[1]), pp * cos(omega[0]), MN);//set 4-momentum of final nucleon
  double result = pow(Tfi(p, pd, L1, L2), 2) * sin(omega[0]) * p.P() * p.E();
  if (isnan(result)){
    std::cout << "Warning: " << omega[0] << " " << omega[1] << std::endl;
    return 0;
  }
  return result;
}

double dsigma(const double * Pd, const double L1, const double L2){//ds / dpd dcostheta
  TLorentzVector pd;
  pd.SetXYZM(Pd[0] * sin(Pd[1]), 0.0, Pd[0] * cos(Pd[1]), Md);
  double par[5] = {pd.P(), pd.Theta(), pd.Phi(), L1, L2};
  double xl[2] = {0.0, -M_PI};//integration lower boundary
  double xu[2] = {M_PI, M_PI};//integration upper boundary
  ROOT::Math::WrappedParamFunction<> wf(&T2, 2, 5);
  wf.SetParameters(par);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  if (isnan(result)){
    std::cout << "Warning: " << Pd[0] << " " << Pd[1] << std::endl;
    return 0;
  }
  return result * pow(2.0 * M_PI, 5) * pow(pd.P(), 2);//in unit GeV^-3
}

double dsigmaT(const double * Pd){//combine ss sp ps pp
  double ss, sp, ps, pp;
  if (C00 != 0) ss = dsigma(Pd, 0, 0);
  else ss = 0.0;
  if (C01 != 0) sp = dsigma(Pd, 0, 1);
  else sp = 0.0;
  if (C10 != 0) ps = dsigma(Pd, 1, 0);
  else ps = 0.0;
  if (C11 != 0) pp = dsigma(Pd, 1, 1);
  else pp = 0.0;
  return C00*ss + C01*sp + C10*ps + C11*pp;//in unit GeV^-3
}

double dsTE(const double * Pd, const double E0){
  double tmp = Ebeam;//same the global Ebeam value
  Ebeam = E0;//change Ebeam to a chosen E0
  double result = dsigmaT(Pd);//calculate the ds with E0
  Ebeam = tmp;//change the global Ebeam value back
  return result;
}

double sigmaTint(const double * Pd){
  double result = dsigmaT(Pd) * sin(Pd[1]);
  if (result < 0.0){
    result = 0;
  }
  return result;
}

double sigmaT(){//Total cross section
  double xl[2] = {0.0, 0.0};//integration lower boundary
  double xu[2] = {1.0, M_PI};//integration upper boundary
  ROOT::Math::Functor wf(&sigmaTint, 2);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 2000);
  //ROOT::Math::GSLMCIntegrator ig(ROOT::Math::IntegrationMultiDim::kVEGAS, 0.0, 0.001, 500);
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
	ds[i][j] = dsigmaT(Pd[i][j]);
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
double decay0(TLorentzVector pd, TLorentzVector * pfs){//NKK decay channel
  double masses[3] = {Mp, Mkaon, Mkaon};
  TGenPhaseSpace event;
  event.SetDecay(pd, 3, masses);
  double weight = event.Generate();
  TLorentzVector * p;
  p = event.GetDecay(0);
  pfs[0].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//p
  p = event.GetDecay(1);
  pfs[1].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K+
  p = event.GetDecay(2);
  pfs[2].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K-
  return weight;
}

double decay1(TLorentzVector pd, TLorentzVector * pfs){//LambdaK decay channel
  TLorentzVector * p;
  double mass1[2] = {MLambda, Mkaon};
  TGenPhaseSpace step1;
  step1.SetDecay(pd, 2, mass1);
  double weight1 = step1.Generate();
  p = step1.GetDecay(1);
  pfs[2].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//K+
  p = step1.GetDecay(0);
  TLorentzVector pL(p->Px(), p->Py(), p->Pz(), p->E());//Lambda
  double mass2[2] = {Mp, Mpion};
  TGenPhaseSpace step2;
  step2.SetDecay(pL, 2, mass2);
  double weight2 = step2.Generate();
  p = step2.GetDecay(0);
  pfs[0].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//p
  p = step2.GetDecay(1);
  pfs[1].SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->E());//pi-
  return weight1 * weight2;
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
double ds_inter(const double * pd){//interpolation of ds/dpd dcostheta
  //pd: pd in GeV, costheta
  return ds2D.Interpolate(pd[0], pd[1]);
}
  
double sigmaT_inter(){//total cross section of bound state production in GeV^-2
  double xl[2] = {0.0, 0.0};//set lower integration boundary
  double xu[2] = {1.0, 1.0};//set upper integration boundary
  ROOT::Math::Functor wf(&ds_inter, 2);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}

//generator
int genData(const char * datafile, const Long64_t Nsim = 10){//Must LoadDS before genData
  TString name = datafile;
  double kmin = 1.2;
  double kmax = 1.8;
  double Ee = 11.0;//Electron beam energy
  TRandom3 r3(4357);//fixed seed
  TF1 * fk = new TF1("fk", "4.0/(3.0*x)-4.0/(3.0*[0])+x/([0]*[0])", kmin, kmax);
  fk->SetParameter(0, Ee);
  fk->SetNpx(500);//Set better resolution for random generator
  TFile * fs = new TFile(name+".root", "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double E0;//photon energy
  double pd[2];//pd, theta
  double phi;
  TLorentzVector Pd;//4-momentum of Pd
  TLorentzVector AP[3];//4-momenta of final particle from NKK channel
  TLorentzVector BP[3];//4-momenta of final particle from LambdaK channel
  double Aw, Bw;//non-normalized weight factors
  double xs;//ds / dPd dcostheta
  Ts->Branch("E0", &E0, "E0/D");
  Ts->Branch("ds", &xs, "ds/D");
  Ts->Branch("Pd", "TLorentzVector", &Pd);
  Ts->Branch("Aw", &Aw, "Aw/D");
  Ts->Branch("AP0", "TLorentzVector", &AP[0]);//p
  Ts->Branch("AP1", "TLorentzVector", &AP[1]);//K+
  Ts->Branch("AP2", "TLorentzVector", &AP[2]);//K-
  Ts->Branch("Bw", &Bw, "Bw/D");
  Ts->Branch("BP0", "TLorentzVector", &BP[0]);//p
  Ts->Branch("BP1", "TLorentzVector", &BP[1]);//pi-
  Ts->Branch("BP2", "TLorentzVector", &BP[2]);//K+
  TH1D * hw = new TH1D("hw", "hw", 2, 0.0, 2.0);
  double Pd_min = 0.0;
  double Pd_max = 1.5;
  double cos_min = 0.0;
  double cos_max = 1.0;
  double volume = (Pd_max - Pd_min) * (cos_max - cos_min);
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%10000 == 9999) std::cout << i + 1 << " in " << Nsim << std::endl;
    E0 = fk->GetRandom();
    pd[0] = gRandom->Uniform(Pd_min, Pd_max);//Generate pd
    pd[1] = acos(gRandom->Uniform(cos_min, cos_max));//Generate theta
    phi = gRandom->Uniform(-M_PI, M_PI);//Generate phi
    xs = dsTE(pd, E0);//Calculate ds / dpd dcostheta
    Pd.SetXYZM(pd[0] * sin(pd[1]) * cos(phi), pd[0] * sin(pd[1]) * sin(phi), pd[0] * cos(pd[1]), Md);//Set 4-momentum of bound state
    Aw = decay0(Pd, AP);//Generate decay from NKK channel
    Bw = decay1(Pd, BP);//Generate decay from LambdaK channel
    hw->Fill(0.5, Aw);
    hw->Fill(1.5, Bw);
    if (xs < 1.0e-9) continue;//skip event with too small ds
    Ts->Fill();
  }
  fs->Write();
  FILE * flog;
  TString logname = name+".log";
  flog = fopen(logname.Data(), "w");
  fprintf(flog, "%d  %.4f  %.1f  %.3f  %.3f  %.6f  %.6f", Nsim, volume, Ee, kmin, kmax, hw->Integral(1,1)/Nsim, hw->Integral(2,2)/Nsim);
  fprintf(flog, "Nsim = %d\n", Nsim);
  fprintf(flog, "Volume = [Pd] * [cos theta] = %.4f GeV\n", volume);
  fprintf(flog, "Electron beam energy = %.1f GeV\n", Ee);
  fprintf(flog, "Photon energy range [kmin, kmax] = [%.3f, %.3f] GeV\n", kmin, kmax);
  fprintf(flog, "Weight normalization W(A) = %.6f, W(B) = %.6f\n", hw->Integral(1,1)/Nsim, hw->Integral(2,2)/Nsim);
  fclose(flog);
  std::cout << "====================================" << std::endl;
  std::cout << "Nsim = " << Nsim << std::endl;
  std::cout << "Volume = [Pd] * [cos theta] = " << volume << std::endl;
  std::cout << "Weight normalization: W(A) = " << hw->Integral(1,1)/Nsim << "  W(B) = " << hw->Integral(2,2)/Nsim << std::endl;
  std::cout << "====================================" << std::endl;
  return 0;
}
  
  



#endif
