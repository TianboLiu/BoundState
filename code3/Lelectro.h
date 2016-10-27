#ifndef _LELECTRO_H_
#define _LELECTRO_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"

TGenPhaseSpace Lphase;//A global variable to generate final state
double _weight_;

/****** Forward tagger cuts ******/
double _Emin_ = 0.5;
double _Emax_ = 4.5;
double _thmin_ = 2.5 / 180.0 * M_PI;
double _thmax_ = 4.5 / 180.0 * M_PI;

int SetTagger(const double Emin, const double Emax, const double thmin = 2.5 / 180.0 * M_PI, const double thmax = 4.5 / 180.0 * M_PI){//Set scattered electron kinematic range
  _Emin_ = Emin;
  _Emax_ = Emax;
  _thmin_ = thmin;
  _thmax_ = thmax;
  return 0;
}
/****** End of tagger cuts ******/

double BreitWigner(const double * M, const double * par){//non-normalized 
  double M0 = par[0];//center mass
  double Gamma = par[1];//decay width
  double E = M[0];//partical mass
  if (E <= 0.0){
    std::cerr << "Unphysical mass in BreitWigner!" << std::endl;
    return -1.0;
  }
  double result = 1.0 / (pow(E*E - M0*M0, 2) + pow(M0*Gamma, 2));
  return result;//non-normalized probability density
}

double BreitWignerNormalization(const double * par){//calculate B-W normalization factor
  double M0 = par[0];
  double Gamma = par[1];
  double Mmin = M0 - 50.0 * Gamma;
  double Mmax = M0 + 50.0 * Gamma;
  if (Mmin < 0) Mmin = 0;
  TF1 f0("B-W function", BreitWigner, Mmin, Mmax, 2);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-4, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(Mmin, Mmax);
  return result;
}

double CarbonMomentum(const double * p0, const double * par){//non-normalized
  double p = p0[0];//nucleon momentum in C12 in unit of GeV
  if (p < 0.0){
    std::cerr << "Unphysical momentum value in CarbonMomentum!" << std::endl;
    return -1.0;
  }
  double A0 = 10.6552;
  double A2 = 32.3105;
  double B2 = 8.43226;
  double result = (A0 + pow(A2 * p, 2)) * exp(-pow(B2 * p, 2));
  return p * p * result / 0.0241519;//Normalized momentum distribution
}

double CarbonEnergy(const double * E0, const double * par){
  double E = E0[0];//nucleon missing energy in C12 in unit of GeV
  if (E <= 0.0){
    std::cerr << "Unphysical energy value in CarbonEnergy!" << std::endl;
    return -1.0;
  }
  double A1 = 0.569161;
  double a1 = 0.0175;
  double b1 = 0.00510586;
  double A2 = 0.0683607;
  double a2 = 0.0369238;
  double b2 = 0.0173035;
  double result = A1 * exp(-pow((E - a1) / b1, 2)) + A2 * exp(-pow((E - a2) / b2, 2));
  return result / 0.00724478;//Normalized missing energy distribution
}

/******* Functions for random sampling *******/
TF1 TF_fp("fp", CarbonMomentum, 0.0, 1.0, 0);//set bound N momentum distri
TF1 TF_fE("fE", CarbonEnergy, 0.0, 0.1, 0);//set bound N missing energy distri
TF1 TF_BWPhi("BWPhi", BreitWigner, 0.819, 1.219, 2);//set phi mass distri
TF1 TF_BWL1520("BWL1520", BreitWigner, 1.43195, 1.600, 2);//set Lambda1520 mass distri
TF1 TF_BWd("BWd", BreitWigner, 1.750, 2.150, 2);//set bound state mass distri
/*** End of Functions for random sampling  ***/

/******* Generator part **********/
int SetFunctions(){
  TF_fp.SetNpx(500);
  TF_fE.SetNpx(500);
  TF_BWPhi.SetParameter(0, 1.019455);//Set phi mass real part
  TF_BWPhi.SetParameter(1, 0.00426);//Set phi width
  TF_BWPhi.SetNpx(2000);
  TF_BWL1520.SetParameter(0, 1.5195);//Set Lambda1520 mass real part
  TF_BWL1520.SetParameter(1, 0.0156);//Set Lambda1520 width
  TF_BWL1520.SetNpx(2000);
  TF_BWd.SetParameter(0, 1.950027);//Set bound state mass real part
  TF_BWd.SetParameter(1, 0.002118);//Set bound state width
  TF_BWd.SetNpx(4000);
  return 0;
}

double Decay2(const TLorentzVector * ki, TLorentzVector * kf, const double * mass){//Generate two body decay event
  double masses[2] = {mass[0], mass[1]};
  TLorentzVector P0 = ki[0];//initial state
  Lphase.SetDecay(P0, 2, masses);
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get 1st particle
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set 1st particle
  tm = Lphase.GetDecay(1);//Get 2nd particle
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set 2nd particle
  return 1;
}

double WeightDistribution3(const double * x, const double * par){//Weight distribution in TGenPhaseSpace
  double E0 = par[0];
  double m0 = par[1];
  double m1 = par[2];
  double m2 = par[3];
  double M0, M1, M2;
  M0 = E0 - m2; M1 = m0; M2 = m1;
  double PMAX0 = sqrt( (M0*M0 - pow(M1+M2, 2)) * (M0*M0 - pow(M1-M2, 2)) ) / (2.0 * M0);
  M0 = E0; M1 = m0 + m1; M2 = m2;
  double PMAX1 = sqrt( (M0*M0 - pow(M1+M2, 2)) * (M0*M0 - pow(M1-M2, 2)) ) / (2.0 * M0);
  M0 = m0 + m1 + x[0] * (E0 - m0 - m1 - m2); M1 = m0; M2 = m1;
  double P0 = sqrt( (M0*M0 - pow(M1+M2, 2)) * (M0*M0 - pow(M1-M2, 2)) ) / (2.0 * M0);
  M0 = E0; M1 = m0 + m1 + x[0] * (E0 - m0 - m1 - m2); M2 = m2;
  double P1 = sqrt( (M0*M0 - pow(M1+M2, 2)) * (M0*M0 - pow(M1-M2, 2)) ) / (2.0 * M0);
  double result = (P0 * P1) / (PMAX0 * PMAX1);
  return result;
} 

double WeightNormalization3(const double M0, const double * mass){//Calculate the weight normalization factor for 3-body decay from TGenPhaseSpace
  double par[4] = {M0, mass[0], mass[1], mass[2]};
  TF1 f0("Weight-3", WeightDistribution3, 0.0, 1.0, 4);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-4, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(0.0, 1.0);
  return result;
}

double Decay3(const TLorentzVector * ki, TLorentzVector * kf, const double * mass){//Generate three body decay event
  double masses[3] = {mass[0], mass[1], mass[2]};
  TLorentzVector P0 = ki[0];//initial state
  Lphase.SetDecay(P0, 3, masses);
  double weight = Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get 1st particle
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set 1st particle
  tm = Lphase.GetDecay(1);//Get 2nd particle
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set 2nd particle
  tm = Lphase.GetDecay(2);//Get 3rd particle
  kf[2].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set 3rd particle
  double Nweight = WeightNormalization3(P0.M(), mass);
  return weight / Nweight;
}

double GetRelativeMomentum(const TLorentzVector * ki){//Get relative momentum of two particles in c.m. frame
  TLorentzVector P = ki[0] + ki[1];//total 4-momentum
  double m1 = ki[0].M();//mass of particle 1
  double m2 = ki[1].M();//mass of particle 2
  double s = P.M2();//Mandelstem variable
  double Q2 = (s - pow(m1 + m2, 2)) * (s - pow(m1 - m2, 2)) / (4.0 * s);
  double Q;
  if (Q2 >= 0.0)
    Q = sqrt(Q2);
  else
    Q = - sqrt(-Q2);
  return Q;
}

double GetTheta0(const TLorentzVector * ki, const TLorentzVector * kf){//Get the polar angle of the produced particle kf[0] in c.m. frame
  TLorentzVector P = ki[0] + ki[1];//total 4-momentum
  double s = P.M2();//Mandelstem variable
  TLorentzVector T = ki[0] - kf[0];//4-momentum transfer
  double t = T.M2();//Mandelstem variable
  double m1 = ki[0].M();//mass of particle 1
  double m2 = ki[1].M();//mass of particle 2
  double p = sqrt( (s - pow(m1 + m2, 2)) * (s - pow(m1 - m2, 2)) / (4.0 * s));//relative momentum of initial state
  double m3 = kf[0].M();//mass of particle 3
  double m4 = kf[1].M();//mass of particle 4
  double Q = sqrt( (s - pow(m3 + m4, 2)) * (s - pow(m3 - m4, 2)) / (4.0 * s));//relative momentum of final state
  double cth = (t - m1*m1 - m3*m3 + 2.0 * sqrt(m1*m1 + p*p) * sqrt(m3*m3 + Q*Q)) / (2.0 * p * Q);//cos theta between p1 and p3 in c.m. frame
  return cth;
}

double GenerateScatteredElectron(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  double Emin = _Emin_;//Set scattered electron energy range
  double Emax = _Emax_;
  double thmin = _thmin_;//Set scattered electron polar angle range
  double thmax = _thmax_;
  double Ee = gRandom->Uniform(Emin, Emax);//Get scattered electron energy
  double theta = acos(gRandom->Uniform(cos(thmax), cos(thmin)));//Get scattered electron polar angle
  double phi = gRandom->Uniform(-M_PI, M_PI);//Get scattered electron azimuthal angle
  kf[0].SetXYZT(Ee * sin(theta) * cos(phi), Ee * sin(theta) * sin(phi), Ee * cos(theta), Ee);//Set scattered electron
  kf[1] = ki[0] - kf[0];//Set virtual photon
  if (kf[1].E() <= 0){//non-physical energy
    weight[0] = 0;
    return 0;
  }
  double alpha_em = 1.0 / 137.0;//fine structure constant
  double Amp = 16.0 * M_PI * alpha_em * (ki[0] * kf[0]);
  double vol = 2.0 * M_PI * (cos(thmin) - cos(thmax)) * (Emax - Emin);
  double Lips = Ee / (2.0 * pow(2.0 * M_PI, 3));//phase space factor in lab frame
  weight[0] = Amp * Lips * vol;
  return weight[0];
}

double AmplitudePhotoproductionPhi(const TLorentzVector * ki, const TLorentzVector * kf){//Calculate phi-meson photoproduction amplitude square
  double par[5] = {0.232612, 1.95038, 4.02454, 1.52884, 0.525636};
  double x = GetTheta0(ki, kf);//theta in c.m. frame
  double M = kf[0].M() + kf[1].M();//threshold
  TLorentzVector Pout = kf[0] + kf[1];//total 4-momentum
  double sr = Pout.M();//energy in c.m. frame
  if (sr < M) return 0.0;//below threshold
  double b1 = par[3] * pow(sr * sr - M * M, par[4]);
  double a2 = par[2] * (sr - M);
  double a0 = par[0] * atan(par[1] * par[1]* (sr * sr - M * M));
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
  double Q = GetRelativeMomentum(kf);//Relative momentum of final state in c.m.
  double Flux = 4.0 * (ki[0] * ki[1]);//Lorentz-invariant relative velosity
  double Lips = Q / (16.0 * M_PI * M_PI * Pout.M());
  double Amp2 = ds * Flux / Lips;//invariant amplitude square in unit 1
  return Amp2;
}

double GeneratePhotoproductionPhiNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  TLorentzVector Pout = ki[0] + ki[1];//Total 4-momentum
  double Mp = 0.938272;//proton mass
  double Mphi = TF_BWPhi.GetRandom();//phi-meson mass
  if (Pout.M() < Mp + Mphi) {
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mphi, Mp};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate event uniform in Omega_k(c.m.)
  kf[0] = *Lphase.GetDecay(0);//Get phi 4-momentum
  kf[1] = *Lphase.GetDecay(1);//Get proton 4-momentum
  double Amp = AmplitudePhotoproductionPhi(ki, kf);//Get invariant amplitude square
  double Flux = 4.0 * (ki[0] * ki[1]);//Lorentz-invariant relative velosity
  double Q = GetRelativeMomentum(kf);
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double vol = 4.0 * M_PI;
  weight[0] = Amp * Lips * vol / Flux;//total amplitude weight
  return weight[0];
}

double GenerateElectroproductionPhiNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  TLorentzVector kf1[2];
  double w1 = GenerateScatteredElectron(ki, kf1);//Get amplitude2 * Lips weight
  kf[0] = kf1[0];//Set final scattered electron
  if (w1 == 0){
    weight[0] = 0;
    return 0;
  }
  double prop = 1.0 / (pow(kf1[1].M2(), 2));//photon propagator square
  TLorentzVector ki2[2] = {kf1[1], ki[1]};//Set virtual photon, proton 4-momentum
  const double Mp = 0.938272;//Proton mass
  double Mphi = TF_BWPhi.GetRandom();//Get phi meson mass
  TLorentzVector Pout = ki2[0] + ki2[1];//total 4-momentum of phi proton final state
  if (Pout.M() < Mphi + Mp){
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mphi, Mp};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate phi proton uniform in Omega_k(c.m.)
  kf[1] = *Lphase.GetDecay(0);//Set phi 4-momentum
  kf[2] = *Lphase.GetDecay(1);//Set proton 4-momentum
  double Amp2 = AmplitudePhotoproductionPhi(ki2, &kf[1]);//virtual photoproducton part amplitude square
  double Q = GetRelativeMomentum(&kf[1]);
  double Lips2 = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double vol2 = 4.0 * M_PI;
  double Flux = 4.0 * (ki[0] * ki[1]);
  weight[0] = w1 * prop * Amp2 * Lips2 * vol2 / Flux;
  return weight[0];
}

double GenerateNucleonInCarbon(TLorentzVector * P){//Generate a bound nucleon in the carbon nucleus
  const double Mp = 0.938272;//proton mass
  double p = TF_fp.GetRandom();//Generate momentum
  double E = sqrt(Mp * Mp + p * p) - TF_fE.GetRandom();//Generate off-shell energy
  double theta = acos(gRandom->Uniform(-1, 1));//Generate polar angle
  double phi = gRandom->Uniform(-M_PI, M_PI);//Generate azimuthal angle
  P[0].SetXYZT(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), E);//Set 4-momentum
  return 1.0;
}

double GeneratePhotoproductionPhiCarbon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  const double Mp = 0.938272;
  const double MA = 12.0 * Mp;
  TLorentzVector ki1[2];
  ki1[0] = ki[0];//Set photon 4-momentum
  GenerateNucleonInCarbon(&ki1[1]);//Set off-shell nucleon 4-momentum
  double Normal = MA / ki1[1].E();//invariant wave function normalization
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of phi proton final state
  double Mphi = TF_BWPhi.GetRandom();//Get phi meson mass
  if (Pout.M() < Mp + Mphi){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mphi, Mp};//Set final state masses
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate phi p uniform in Omega_k(c.m.)
  kf[0] = *Lphase.GetDecay(0);//Get phi 4-momentum
  kf[1] = *Lphase.GetDecay(1);//Get proton 4-momentum
  double Amp = AmplitudePhotoproductionPhi(ki1, kf);//Get invariant amplitude square
  double Q = GetRelativeMomentum(kf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Normal * Amp * Lips * vol / Flux;//total weight
  return weight[0];
}
				  
double tt(const double * x){
  //x[0] = 1;
  return 0;
}




#endif
