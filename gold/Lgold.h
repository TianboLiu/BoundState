#ifndef _LGOLD_H_
#define _LGOLD_H_

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
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
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

/******* Functions for random sampling *******/
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

double GoldMomentum(const double * p0, const double * par){//non-normalized
  double p = p0[0];//nucleon momentum in Au197 in unit of GeV
  if (p < 0.0){
    std::cerr << "Unphysical momentum value in GoldMomentum!" << std::endl;
    return -1.0;
  }
  double A0 = 58.3382;
  double A2 = 69.2938;
  double B2 = 7.82756;
  double result = (A0 + pow(A2 * p, 2)) * exp(-pow(B2 * p, 2));
  return p * p * result / 0.162508;//Normalized momentum distribution
}

double GoldEnergy(const double * E0, const double * par){
  double E = E0[0];//nucleon missing energy in Au197 in unit of GeV
  if (E <= 0.0){
    std::cerr << "Unphysical energy value in GoldEnergy!" << std::endl;
    return -1.0;
  }
  double A1 = 1.73622;
  double a1 = 3.07375;
  double b1 = 0.645561;
  double A2 = 14.1433;
  double a2 = 0.795058;
  double result = A1 * atan(A2 * pow(E/0.01, a1)) * exp(-b1 * pow(E/0.01, a2));
  return result / 0.0433967;//Normalized missing energy distribution
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

TF1 TF_fp("fp", GoldMomentum, 0.0, 1.0, 0);//set bound N momentum distri
TF1 TF_fE("fE", GoldEnergy, 0.0, 0.3, 0);//set bound N missing energy distri
TF1 TF_BWPhi("BWPhi", BreitWigner, 0.819, 1.219, 2);//set phi mass distri
TF1 TF_BWL1520("BWL1520", BreitWigner, 1.43195, 1.600, 2);//set Lambda1520 mass distri
TF1 TF_BWd("BWd", BreitWigner, 1.900, 2.000, 2);//set bound state mass distri

int SetFunctions(){
  TF_fp.SetNpx(500);
  TF_fE.SetNpx(1000);
  TF_BWPhi.SetParameter(0, 1.019455);//Set phi mass real part
  TF_BWPhi.SetParameter(1, 0.00426);//Set phi width
  TF_BWPhi.SetNpx(2000);
  TF_BWL1520.SetParameter(0, 1.5195);//Set Lambda1520 mass real part
  TF_BWL1520.SetParameter(1, 0.0156);//Set Lambda1520 width
  TF_BWL1520.SetNpx(2000);
  TF_BWd.SetParameter(0, 1.950027);//Set bound state mass real part
  TF_BWd.SetParameter(1, 0.004094);//Set bound state width
  TF_BWd.SetNpx(4000);
  return 0;
}
/*** End of Functions for random sampling  ***/

/****** Functions for general use ******/
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

double GetTheta0(const TLorentzVector * ki, const TLorentzVector * kf){//Get the polar angle of the produced particle kf[0] in c.m. frame w.r.t ki[0]
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
  const double Me = 0.000510998;//electron mass
  kf[0].SetXYZM(Ee * sin(theta) * cos(phi), Ee * sin(theta) * sin(phi), Ee * cos(theta), Me);//Set scattered electron
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


double GenerateNucleonInGold(TLorentzVector * P){//Generate a bound nucleon in the gold nucleus
  const double Mp = 0.938272;//proton mass
  double p = TF_fp.GetRandom();//Generate momentum
  double E = sqrt(Mp * Mp + p * p) - TF_fE.GetRandom();//Generate off-shell energy
  double theta = acos(gRandom->Uniform(-1, 1));//Generate polar angle
  double phi = gRandom->Uniform(-M_PI, M_PI);//Generate azimuthal angle
  P[0].SetXYZT(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), E);//Set 4-momentum
  return 1.0;
}

/****** Production channels ******/

//phi meson production part
double AmplitudePhotoproductionPhi(const TLorentzVector * ki, const TLorentzVector * kf){//Calculate phi-meson photoproduction amplitude square
  double par[5] = {0.232612, 1.95038, 4.02454, 1.52884, 0.525636};
  double x = GetTheta0(ki, kf);//theta in c.m. frame
  double M = kf[0].M() + kf[1].M();//threshold
  TLorentzVector Pout = kf[0] + kf[1];//total 4-momentum
  double sr = Pout.M();//energy in c.m. frame
  if (sr < M) return 0;//below threshold
  double b1 = par[3] * pow(sr * sr - M * M, par[4]);
  double a2 = par[2] * (sr - M);
  double a0 = par[0] * atan(par[1] * par[1]* (sr * sr - M * M));
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
  double Q = GetRelativeMomentum(kf);//Relative momentum of final state in c.m.
  if (Q < 1.0e-30) return 0;
  double Flux = 4.0 * (ki[0] * ki[1]);//Lorentz-invariant relative velosity
  double Lips = Q / (16.0 * M_PI * M_PI * Pout.M());
  double Amp2 = ds * Flux / Lips;//invariant amplitude square in unit 1
  return Amp2;
}

double GeneratePhotoproductionPhiNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: photon, p;  kf: phi, p
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
  //ki: e, p;  kf: e', phi, p
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

double GeneratePhotoproductionPhiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: photon;  kf: phi, p
  const double Mp = 0.938272;
  const double MA = 197.0 * Mp;
  TLorentzVector ki1[2];
  ki1[0] = ki[0];//Set photon 4-momentum
  GenerateNucleonInGold(&ki1[1]);//Set off-shell nucleon 4-momentum
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

double GenerateElectroproductionPhiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: e;  kf: e', phi, p
  TLorentzVector kf1[2];
  double w1 = GenerateScatteredElectron(ki, kf1);//Get amplitude2 * Lips weight
  kf[0] = kf1[0];//Set final scattered electron
  if (w1 == 0){//unphysical virtual photon energy
    weight[0] = 0;
    return 0;
  }
  double prop = 1.0 / pow(kf1[1].M2(), 2);//virtual photon propagator square
  TLorentzVector ki2[2];
  ki2[0] = kf1[1];//Set intial virtual photon for phi-production step
  GenerateNucleonInGold(&ki2[1]);//Set off-shell nucleon 4-momentum
  const double Mp = 0.938272;//proton mass
  const double Mphi = TF_BWPhi.GetRandom();//Get phi meson mass
  TLorentzVector Pout = ki2[0] + ki2[1];//total 4-momentum of phi proton final state
  if (Pout.M() < Mphi + Mp){//below threshold
    weight[0] = 0;
    return 0;
  }
  const double MA = 197.0 * Mp;
  double Normal = MA / ki2[1].E();//invariant wave function normalization
  double mass[2] = {Mphi, Mp};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate phi proton uniform in Omega_k(c.m.);
  kf[1] = *Lphase.GetDecay(0);//Set phi meson 4-momentum
  kf[2] = *Lphase.GetDecay(1);//Set proton 4-momentum
  double Amp2 = AmplitudePhotoproductionPhi(ki2, &kf[1]);//Get amplitude square of phi production part
  double Q = GetRelativeMomentum(&kf[1]);//Get relative momentum of phi proton in c.m. frame
  double vol2 = 4.0 * M_PI;
  double Lips2 = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Normal * w1 * prop * Amp2 * Lips2 * vol2 / Flux;//total weight
  return weight[0];
}

//direct K+K- production part
double GeneratePhotoproductionKKNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: photon, p1;  kf: p, K+, K-
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  TLorentzVector Pout = ki[0] + ki[1];//total 4-momentum of pKK
  double mass[3] = {Mp, Mkaon, Mkaon};
  if (Pout.M() < Mp + Mkaon + Mkaon){//below threshold
    weight[0] = 0;
    return 0;
  }
  double wd = Decay3(&Pout, kf, mass);//Get phase space weight factor
  TLorentzVector tf[2] = {kf[0], kf[1] + kf[2]};//group K+K-
  double Amp = AmplitudePhotoproductionPhi(ki, tf);//Get amplitude square
  double Q = GetRelativeMomentum(tf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * (ki[0] * ki[1]);
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = wd * Br * Amp * vol * Lips / Flux;//total weight
  return weight[0];
}

double GenerateElectroproductionKKNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: e, p1;  kf: e', p, K+, K-
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  TLorentzVector kf0[2];
  double w0 = GenerateScatteredElectron(ki, kf0);//Get amp * Lips weight
  if (w0 == 0){//unphysical
    weight[0] = 0;
    return 0;
  }
  kf[0] = kf0[0];//Set scattered electron
  double prop = 1.0 / pow(kf0[1].M2(), 2);//virtual photon propagator square
  TLorentzVector ki1[2] = {kf0[1], ki[1]};
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of pKK
  double mass[3] = {Mp, Mkaon, Mkaon};
  if (Pout.M() < Mp + Mkaon + Mkaon){//below threshold
    weight[0] = 0;
    return 0;
  }
  double wd = Decay3(&Pout, &kf[1], mass);//Get phase space weight factor
  TLorentzVector tf[2] = {kf[1], kf[2] + kf[3]};//group K+K-
  double Amp = AmplitudePhotoproductionPhi(ki1, tf);//Get amplitude square
  double Q = GetRelativeMomentum(tf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * (ki[0] * ki[1]);
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = w0 * prop * wd * Br * Amp * vol * Lips / Flux;//total weight
  return weight[0];
}

double GeneratePhotoproductionKKGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: photon;  kf: p, K+, K-
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  const double MA = 197.0 * Mp;
  TLorentzVector ki1[2];
  ki1[0] = ki[0];//photon
  GenerateNucleonInGold(&ki1[1]);//reaction nucleon
  double Normal = MA / ki1[1].E();//invariant wave function normalization
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of pKK
  double mass[3] = {Mp, Mkaon, Mkaon};
  if (Pout.M() < Mp + Mkaon + Mkaon){//below threshold
    weight[0] = 0;
    return 0;
  }
  double wd = Decay3(&Pout, kf, mass);//Get phase space weight factor
  TLorentzVector tf[2] = {kf[0], kf[1] + kf[2]};//group K+K-
  double Amp = AmplitudePhotoproductionPhi(ki1, tf);//Get amplitude square
  double Q = GetRelativeMomentum(tf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * (ki[0].E() * MA);
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = Normal * wd * Br * Amp * vol * Lips / Flux;//total weight
  return weight[0];
}

double GenerateElectroproductionKKGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: e;  kf: e', p, K+, K-
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  const double MA = 197.0 * Mp;
  TLorentzVector kf0[2];
  double w0 = GenerateScatteredElectron(ki, kf0);//Get amp * Lips weight
  if (w0 == 0){//unphysical
    weight[0] = 0;
    return 0;
  }
  kf[0] = kf0[0];//Set scattered electron
  double prop = 1.0 / pow(kf0[1].M2(), 2);//virtual photon propagator square
  TLorentzVector ki1[2];
  ki1[0] = kf0[1];//virtual photon
  GenerateNucleonInGold(&ki1[1]);
  double Normal = MA / ki1[1].E();//invariant wave function normalization
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of pKK
  double mass[3] = {Mp, Mkaon, Mkaon};
  if (Pout.M() < Mp + Mkaon + Mkaon){//below threshold
    weight[0] = 0;
    return 0;
  }
  double wd = Decay3(&Pout, &kf[1], mass);//Get phase space weight factor
  TLorentzVector tf[2] = {kf[1], kf[2] + kf[3]};//group K+K-
  double Amp = AmplitudePhotoproductionPhi(ki1, tf);//Get amplitude square
  double Q = GetRelativeMomentum(tf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * (ki[0].E() * MA);
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = Normal * w0 * prop * wd * Br * Amp * vol * Lips / Flux;//total weight
  return weight[0];
}

//Lambda1520 production part
double AmplitudePhotoproductionLambda1520(const TLorentzVector * ki, const TLorentzVector * kf){//Calculate K+ Lambda1520 photoproduction amplitude square
  //ki: q, p1; K, Lambda
  double par[5] = {11.299, 4.60959, 0.835621, 0.54681, 1.827941};
  double x = GetTheta0(ki, kf);//theta in c.m. frame
  double M = kf[0].M() + kf[1].M();//threshold
  TLorentzVector Pout = kf[0] + kf[1];//total 4-momentum
  double sr = Pout.M();//energy in c.m. frame
  if (sr < M) return 0;//below threshold
  double b1 = par[3] * (sr - M + par[4]);
  double a2 = 0.25;
  double a0 = par[0] * (sr - M) * exp(-par[1] * pow(sr - M, par[2]));;//total cross section in unit mb
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
  double Q = GetRelativeMomentum(kf);//Relative momentum of final state in c.m.
  if (Q < 1.0e-30) return 0;
  double Flux = 4.0 * (ki[0] * ki[1]);//Lorentz-invariant relative velosity
  double Lips = Q / (16.0 * M_PI * M_PI * Pout.M());
  double Amp2 = ds * Flux / Lips;//invariant amplitude square in unit 1
  return Amp2;
}

double GeneratePhotoproductionLambda1520Nucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: photon, p;  kf: K+, Lambda
  TLorentzVector Pout = ki[0] + ki[1];//Total 4-momentum
  const double MLambda = TF_BWL1520.GetRandom();//Lambda mass
  const double Mkaon = 0.493677;//kaon mass 
  if (Pout.M() < Mkaon + MLambda) {//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mkaon, MLambda};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate event uniform in Omega_k(c.m.)
  kf[0] = *Lphase.GetDecay(0);//Get K+ 4-momentum
  kf[1] = *Lphase.GetDecay(1);//Get Lambda1520 4-momentum
  double Amp = AmplitudePhotoproductionLambda1520(ki, kf);//Get invariant amplitude square
  double Flux = 4.0 * (ki[0] * ki[1]);//Lorentz-invariant relative velosity
  double Q = GetRelativeMomentum(kf);
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double vol = 4.0 * M_PI;
  weight[0] = Amp * Lips * vol / Flux;//total amplitude weight
  return weight[0];
}

double GenerateElectroproductionLambda1520Nucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: e, p;  kf: e', K+, Lambda
  TLorentzVector kf1[2];
  double w1 = GenerateScatteredElectron(ki, kf1);//Get amplitude2 * Lips weight
  kf[0] = kf1[0];//Set final scattered electron
  if (w1 == 0){
    weight[0] = 0;
    return 0;
  }
  double prop = 1.0 / (pow(kf1[1].M2(), 2));//photon propagator square
  TLorentzVector ki2[2] = {kf1[1], ki[1]};//Set virtual photon, proton 4-momentum
  const double Mkaon = 0.493677;//kaon mass
  double MLambda = TF_BWL1520.GetRandom();//Get Lambda mass
  TLorentzVector Pout = ki2[0] + ki2[1];//total 4-momentum of K+ L1520 final state
  if (Pout.M() < Mkaon + MLambda){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mkaon, MLambda};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate K+ Lambda uniform in Omega_k(c.m.)
  kf[1] = *Lphase.GetDecay(0);//Set K+ 4-momentum
  kf[2] = *Lphase.GetDecay(1);//Set Lambda1520 4-momentum
  double Amp2 = AmplitudePhotoproductionLambda1520(ki2, &kf[1]);//virtual photoproducton part amplitude square
  double Q = GetRelativeMomentum(&kf[1]);
  double Lips2 = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double vol2 = 4.0 * M_PI;
  double Flux = 4.0 * (ki[0] * ki[1]);
  weight[0] = w1 * prop * Amp2 * Lips2 * vol2 / Flux;
  return weight[0];
}

double GeneratePhotoproductionLambda1520Gold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: photon;  kf: K+, Lambda
  const double Mp = 0.938272;
  const double Mkaon = 0.493677;
  const double MA = 197.0 * Mp;
  TLorentzVector ki1[2];
  ki1[0] = ki[0];//Set photon 4-momentum
  GenerateNucleonInGold(&ki1[1]);//Set off-shell nucleon 4-momentum
  double Normal = MA / ki1[1].E();//invariant wave function normalization
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of K+ L1520 final state
  double MLambda = TF_BWL1520.GetRandom();//Get Lambda1520 mass
  if (Pout.M() < Mkaon + MLambda){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mkaon, MLambda};//Set final state masses
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate K+ Lambda uniform in Omega_k(c.m.)
  kf[0] = *Lphase.GetDecay(0);//Get K+ 4-momentum
  kf[1] = *Lphase.GetDecay(1);//Get Lambda1520 4-momentum
  double Amp = AmplitudePhotoproductionLambda1520(ki1, kf);//Get invariant amplitude square
  double Q = GetRelativeMomentum(kf);
  double vol = 4.0 * M_PI;
  double Lips = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Normal * Amp * Lips * vol / Flux;//total weight
  return weight[0];
}

double GenerateElectroproductionLambda1520Gold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: e;  kf: e', K, Lambda
  TLorentzVector kf1[2];
  double w1 = GenerateScatteredElectron(ki, kf1);//Get amplitude2 * Lips weight
  kf[0] = kf1[0];//Set final scattered electron
  if (w1 == 0){//unphysical virtual photon energy
    weight[0] = 0;
    return 0;
  }
  double prop = 1.0 / pow(kf1[1].M2(), 2);//virtual photon propagator square
  TLorentzVector ki2[2];
  ki2[0] = kf1[1];//Set intial virtual photon for K+ L1520 production step
  GenerateNucleonInGold(&ki2[1]);//Set off-shell nucleon 4-momentum
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  const double MLambda = TF_BWL1520.GetRandom();//Get Lambda mass
  TLorentzVector Pout = ki2[0] + ki2[1];//total 4-momentum of Lambda final state
  if (Pout.M() < Mkaon + MLambda){//below threshold
    weight[0] = 0;
    return 0;
  }
  const double MA = 197.0 * Mp;
  double Normal = MA / ki2[1].E();//invariant wave function normalization
  double mass[2] = {Mkaon, MLambda};
  Lphase.SetDecay(Pout, 2, mass);
  Lphase.Generate();//Generate K+ Lambda uniform in Omega_k(c.m.);
  kf[1] = *Lphase.GetDecay(0);//Set K+ meson 4-momentum
  kf[2] = *Lphase.GetDecay(1);//Set Lambda 4-momentum
  double Amp2 = AmplitudePhotoproductionLambda1520(ki2, &kf[1]);//Get amplitude square of K+ L1520 production part
  double Q = GetRelativeMomentum(&kf[1]);//Get relative momentum of K+ in c.m. frame
  double vol2 = 4.0 * M_PI;
  double Lips2 = Q / (4.0 * Pout.M() * pow(2.0 * M_PI, 2));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Normal * w1 * prop * Amp2 * Lips2 * vol2 / Flux;//total weight
  return weight[0];
}

//Bound state production part
double AmplitudeBoundStateFormation(const TLorentzVector * ki, const TLorentzVector * kf){//phi M -> d invariant amplitude (not square)
  //ki: phi, p;  kf: d
  double Q = GetRelativeMomentum(ki);//Get phi N relative momentum in c.m. frame
  if (Q < 0){//unphysical momentum
    return 0;
  }
  double A0 = 0.112356;//Fit to FQ0.dat with poly*gaussian
  double A2 = 0.609078;
  double B2 = 2.75841;
  double FQ = (A0 + A2 * Q * Q) * exp(-B2 * Q * Q);//without cutoff
  double Ep = sqrt(ki[1].M2() + Q*Q);
  double Normal = 2.0 * Ep * sqrt(2.0 * kf[0].M());
  return Normal * FQ;//invariant amplitude in unit GeV
}

double GeneratePhotoproductionBoundStateGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: q; kf: p, d
  const double Mp = 0.938272;
  const double MA = 197.0 * Mp;
  //step 1
  TLorentzVector p1;
  GenerateNucleonInGold(&p1);//Get reaction nucleon-1
  TLorentzVector ki1[2] = {ki[0], p1};//Set phi production part initial state
  const double Mphi = TF_BWPhi.GetRandom(0.976, 1.062);//Get phi mass
  const double Normal_phi = 352.455;//phi meson propagator normalization
  double mass[2] = {Mphi, Mp};
  TLorentzVector Pout1 = ki1[0] + ki1[1];//total 4-momentum of phi p system
  if (Pout1.M() < Mphi + Mp){//below phi threshold
    weight[0] = 0;
    return 0;
  }
  Lphase.SetDecay(Pout1, 2, mass);
  Lphase.Generate();
  kf[0] = *Lphase.GetDecay(1);//Set final proton
  TLorentzVector k = *Lphase.GetDecay(0);//Get intermediate phi meson
  double Jac_phi = k.M() / k.E();//jacobian dE_phi / dM_phi
  TLorentzVector kf1[2] = {k, kf[0]};//phi production part final state
  double Amp1 = AmplitudePhotoproductionPhi(ki1, kf1);//Get invariant amplitude square
  double w1 = Amp1 * Jac_phi * Normal_phi;//total weight of step 1
  //step 2
  TLorentzVector p2;
  GenerateNucleonInGold(&p2);//Get reaction nucleon-2
  TLorentzVector ki2[2] = {k, p2};//Set bound state formation part initial state
  TLorentzVector Pout2 = ki2[0] + ki2[1];//total 4-momentum of final bound state
  if (Pout2.M() <= TF_BWd.GetXmin() || Pout2.M() >= TF_BWd.GetXmax()){//fail to match bound state mass
    weight[0] = 0;
    return 0;
  }
  kf[1] = Pout2;//Set final bound state 4-momentum
  double Amp2 = pow(AmplitudeBoundStateFormation(ki2, &Pout2), 2);//Get invariant amplitude square of bound state formation
  const double Normal_d = 98.2728;//bound state mass distribution normalization
  double distri_Md = TF_BWd.Eval(Pout2.M()) / Normal_d;//mass matching weight
  double Q = GetRelativeMomentum(kf1);//Get relative momentum of phi p system in c.m. frame
  double Ep0 = sqrt(Mp * Mp + Q * Q);
  TLorentzVector p0 = kf[0];
  p0.Boost(-Pout1.BoostVector());//Get final p in phi p cm frame
  double dEpdQ = Pout1.Gamma() * ( Q / Ep0 + Pout1.Beta() * cos(p0.Angle(Pout1.Vect())) );
  double dEddMd = Pout2.M() / Pout2.E();
  double Jac_d = std::abs(1.0 / dEpdQ / dEddMd);
  double w2 = Amp2 * distri_Md * Jac_d;//total weight of step 2
  double Nor = MA * (MA - p1.E() - p2.E());
  double vol = (4.0 * M_PI) / ( pow(2.0 * M_PI, 3) * 2.0 * Ep0) / ( pow(2.0 * M_PI, 3) * 2.0 * Pout2.E()) / ( pow(2.0 * M_PI, 3) * 2.0 * (MA - p1.E() - p2.E()));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Nor * w1 * w2 * vol / Flux;//total weight
  return weight[0];
}
  
double GenerateElectroproductionBoundStateGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//
  //ki: e;  kf: e', p, d
  const double Mp = 0.938272;
  const double MA = 197.0 * Mp;
  //step 0
  TLorentzVector kf0[2];
  double w0 = GenerateScatteredElectron(ki, kf0);//Get weight * Lips factor of virtual photon production part
  if (w0 == 0){//unphysical event
    weight[0] = 0;
    return 0;
  }
  kf[0] = kf0[0];//Set final electron
  double prop = 1.0 / pow(kf0[1].M2(), 2);//Get virtual photon propagater square
  //step 1
  TLorentzVector p1;
  GenerateNucleonInGold(&p1);//Get reaction nucleon-1
  TLorentzVector ki1[2] = {kf0[1], p1};//Set phi production part initial state
  const double Mphi = TF_BWPhi.GetRandom(0.976, 1.062);//Get phi mass
  const double Normal_phi = 352.455;//phi meson propagator normalization
  double mass[2] = {Mphi, Mp};
  TLorentzVector Pout1 = ki1[0] + ki1[1];//total 4-momentum of phi p system
  if (Pout1.M() < Mphi + Mp){//below phi threshold
    weight[0] = 0;
    return 0;
  }
  Lphase.SetDecay(Pout1, 2, mass);
  Lphase.Generate();
  kf[1] = *Lphase.GetDecay(1);//Set final proton
  TLorentzVector k = *Lphase.GetDecay(0);//Get intermediate phi meson
  double Jac_phi = k.M() / k.E();//jacobian dE_phi / dM_phi
  TLorentzVector kf1[2] = {k, kf[0]};//phi production part final state
  double Amp1 = AmplitudePhotoproductionPhi(ki1, kf1);//Get invariant amplitude square
  double w1 = Amp1 * Jac_phi * Normal_phi;//total weight of step 1
  //step 2
  TLorentzVector p2;
  GenerateNucleonInGold(&p2);//Get reaction nucleon-2
  TLorentzVector ki2[2] = {k, p2};//Set bound state formation part initial state
  TLorentzVector Pout2 = ki2[0] + ki2[1];//total 4-momentum of final bound state
  if (Pout2.M() <= TF_BWd.GetXmin() || Pout2.M() >= TF_BWd.GetXmax()){//fail to match bound state mass
    weight[0] = 0;
    return 0;
  }
  kf[2] = Pout2;//Set final bound state 4-momentum
  double Amp2 = pow(AmplitudeBoundStateFormation(ki2, &Pout2), 2);//Get invariant amplitude square of bound state formation
  const double Normal_d = 98.2728;//bound state mass distribution normalization
  double distri_Md = TF_BWd.Eval(Pout2.M()) / Normal_d;//mass matching weight
  double Q = GetRelativeMomentum(kf1);//Get relative momentum of phi p system in c.m. frame
  double Ep0 = sqrt(Mp * Mp + Q * Q);
  TLorentzVector p0 = kf[0];
  p0.Boost(-Pout1.BoostVector());//Get final p in phi p cm frame
  double dEpdQ = Pout1.Gamma() * ( Q / Ep0 + Pout1.Beta() * cos(p0.Angle(Pout1.Vect())) );
  double dEddMd = Pout2.M() / Pout2.E();
  double Jac_d = std::abs(1.0 / dEpdQ / dEddMd);
  double w2 = Amp2 * distri_Md * Jac_d;//total weight of step 2
  double Nor = MA * (MA - p1.E() - p2.E());
  double vol = (4.0 * M_PI) / ( pow(2.0 * M_PI, 3) * 2.0 * Ep0) / ( pow(2.0 * M_PI, 3) * 2.0 * Pout2.E()) / ( pow(2.0 * M_PI, 3) * 2.0 * (MA - p1.E() - p2.E()));
  double Flux = 4.0 * ki[0].E() * MA;
  weight[0] = Nor * w0 * prop * w1 * w2 * vol / Flux;//total weight
  return weight[0];
}


/****** Event generator part ******/

double GenerateEvent_NKK_withPhotoproductionPhiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with NKK in final state with phi-meson production channel
  //ki: q; kf: N, K+, K-
  const double Mkaon = 0.493677;//kaon mass
  const double scale = 197.0;//nucleon number scale
  TLorentzVector kf1[2];//intermediate state: phi, N
  double w1 = GeneratePhotoproductionPhiGold(ki, kf1);//Get phi-meson electroproduction weight
  if (kf1[0].M() < Mkaon + Mkaon){//below K+K- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mkaon, Mkaon};
  TLorentzVector kf2[2];
  double w2 = Decay2(&kf1[0], kf2, mass);//Get K+K- decayed from phi-meson
  kf[0] = kf1[1];//Set proton
  kf[1] = kf2[0];//Set K+
  kf[2] = kf2[1];//Set K-
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}

double GenerateEvent_NKK_withPhotoproductionKKGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with NKK in final state with direct K+K- production channel
  //ki: q; kf: N, K+, K-
  const double scale = 197.0;//nucleon number scale
  double w1 = GeneratePhotoproductionKKGold(ki, kf);//Get electroproduction KK weight
  weight[0] = scale * w1;
  return weight[0];
}

double GenerateEvent_NKK_withPhotoproductionLambda1520Gold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with NKK in final state with Lambda1520 production channel
  //ki: q; kf: N, K+, K-
  const double scale = 197.0;//nucleon number scale
  const double Mkaon = 0.493677;//kaon mass
  const double Mp = 0.938272;//proton mass
  TLorentzVector kf1[2];//intermediate state: K+, Lambda1520
  double w1 = GeneratePhotoproductionLambda1520Gold(ki, kf1);//Get electroproduction Lambda1520 weight
  if (kf1[1].M() < Mkaon + Mp){//below pK- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mp, Mkaon};
  TLorentzVector kf2[2];
  double w2 = Decay2(&kf1[1], kf2, mass);//Get p K- decayed from Lambda1520
  kf[0] = kf2[0];//Set proton
  kf[1] = kf1[0];//Set K+
  kf[2] = kf2[1];//Set K-
  double Br = 0.45 / 2.0;//Branch ratio of Lambda1520 to pK-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}

double GenerateEvent_NKKN_withPhotoproductionBoundStateGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with NKKN in final state with bound state production channel
  //ki: 1; kf: N, (K+, K-, N)
  const double scale = 197.0 * 196.0;//nucleon number scale
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  TLorentzVector kf1[2];//intermediate state: p, d
  double w1 = GeneratePhotoproductionBoundStateGold(ki, kf1);//Get electroproduction bound state weight
  if (kf1[1].M() < Mp + Mkaon + Mkaon){//below N K+ K- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mkaon, Mkaon};
  TLorentzVector kf2[3];
  double w2 = Decay3(&kf1[1], kf2, mass);//Get N K+ K- decayed from bound state
  kf[0] = kf1[0];//Set produced proton
  kf[1] = kf2[1];//Set K+
  kf[2] = kf2[2];//Set K-
  kf[3] = kf2[0];//Set proton decayed from bound state
  double Br = 0.951 * 0.489;//Branch ratio of bound state to NK+K-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}

double GenerateEvent_eNKK_withElectroproductionPhiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with eNKK in final state with phi-meson production channel
  //ki: e; kf: e', N, K+, K-
  const double Mkaon = 0.493677;//kaon mass
  const double scale = 197.0;//nucleon number scale
  TLorentzVector kf1[3];//intermediate state: e', phi, N
  double w1 = GenerateElectroproductionPhiGold(ki, kf1);//Get phi-meson electroproduction weight
  if (kf1[1].M() < Mkaon + Mkaon){//below K+K- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mkaon, Mkaon};
  TLorentzVector kf2[2];
  double w2 = Decay2(&kf1[1], kf2, mass);//Get K+K- decayed from phi-meson
  kf[0] = kf1[0];//Set scattered electron
  kf[1] = kf1[2];//Set proton
  kf[2] = kf2[0];//Set K+
  kf[3] = kf2[1];//Set K-
  double Br = 0.489;//Branch ratio of phi to K+K-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}

double GenerateEvent_eNKK_withElectroproductionKKGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with eNKK in final state with direct K+K- production channel
  //ki: e; kf: e', N, K+, K-
  const double scale = 197.0;//nucleon number scale
  double w1 = GenerateElectroproductionKKGold(ki, kf);//Get electroproduction KK weight
  weight[0] = scale * w1;
  return weight[0];
}

double GenerateEvent_eNKK_withElectroproductionLambda1520Gold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with eNKK in final state with Lambda1520 production channel
  //ki: e; kf: e', N, K+, K-
  const double scale = 197.0;//nucleon number scale
  const double Mkaon = 0.493677;//kaon mass
  const double Mp = 0.938272;//proton mass
  TLorentzVector kf1[3];//intermediate state: e', K+, Lambda1520
  double w1 = GenerateElectroproductionLambda1520Gold(ki, kf1);//Get electroproduction Lambda1520 weight
  if (kf1[2].M() < Mkaon + Mp){//below pK- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[2] = {Mp, Mkaon};
  TLorentzVector kf2[2];
  double w2 = Decay2(&kf1[2], kf2, mass);//Get p K- decayed from Lambda1520
  kf[0] = kf1[0];//Set scattered electron
  kf[1] = kf2[0];//Set proton
  kf[2] = kf1[1];//Set K+
  kf[3] = kf2[1];//Set K-
  double Br = 0.45 / 2.0;//Branch ratio of Lambda1520 to pK-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}

double GenerateEvent_eNKKN_withElectroproductionBoundStateGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with eNKKN in final state with bound state production channel
  //ki: e; kf: e', N, (K+, K-, N)
  const double scale = 197.0 * 196.0;//nucleon number scale
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//kaon mass
  TLorentzVector kf1[3];//intermediate state: e', p, d
  double w1 = GenerateElectroproductionBoundStateGold(ki, kf1);//Get electroproduction bound state weight
  if (kf1[2].M() < Mp + Mkaon + Mkaon){//below N K+ K- threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mkaon, Mkaon};
  TLorentzVector kf2[3];
  double w2 = Decay3(&kf1[2], kf2, mass);//Get N K+ K- decayed from bound state
  kf[0] = kf1[0];//Set scattered electron
  kf[1] = kf1[1];//Set produced proton
  kf[2] = kf2[1];//Set K+
  kf[3] = kf2[2];//Set K-
  kf[4] = kf2[0];//Set proton decayed from bound state
  double Br = 0.951 * 0.489;//Branch ratio of bound state to NK+K-
  weight[0] = scale * w1 * w2 * Br;
  return weight[0];
}
				    
/****** detector part ******/

int DetectorResolutionSmear(TLorentzVector * P, double * res){//Smearing the three momentum with gaussian, fixed mass
  double M = P->M();
  double p = P->P() * (1.0 + gRandom->Gaus(0.0, res[0]));//smear p
  if (p < 0.0) p = 0.0;//Set negative p to 0
  double theta = P->Theta() + gRandom->Gaus(0.0, res[1]);//smear theta
  double phi = P->Phi() + gRandom->Gaus(0.0, res[2]);//smear phi
  P[0].SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), M);//Set smeared three momentum with fixed mass
  return 0;
}

double LongLifetimeDecayFactor(const TLorentzVector * P, const double ct, const double dis = 2.0){//Calculate the remaining factor when arriving at the detector due to decay
  double beta = P->Beta();
  double df = exp(- dis / beta / ct);
  return df;
}

bool MomentumCut(const TLorentzVector * P, const double pmin, const double pmax){//check the momentum in a selected region
  double p = P->P();//in unit GeV
  if (p < pmin || p > pmax)
    return false;
  else
    return true;
}

bool PolarAngleCut(const TLorentzVector * P, const double thmin, const double thmax, const char * unit = "deg"){//check the polar angle in a selected region
  double th = P->Theta();
  if (strcmp(unit, "deg") == 0)
    th = th * 180.0 / M_PI;
  if (th < thmin || th > thmax)
    return false;
  else
    return true;
}

bool MassCut(const TLorentzVector * P, const double Mmin, const double Mmax){//check the invariant mass in a selected region
  double M = P->M();//in unit GeV
  if (M < Mmin || M > Mmax)
    return false;
  else
    return true;
}

const double target_length = 0.01383;//in unit cm
const double target_radius = 0.3;//in unit cm
const int NPT_P_GOLD = 132;
double EK_P_GOLD[132], dEK_P_GOLD[132];
ROOT::Math::Interpolator SP_P_GOLD(NPT_P_GOLD, ROOT::Math::Interpolation::kCSPLINE);
const int NPT_K_GOLD = 132;
double EK_K_GOLD[132], dEK_K_GOLD[132];
ROOT::Math::Interpolator SP_K_GOLD(NPT_K_GOLD, ROOT::Math::Interpolation::kCSPLINE);

int SetStoppingPower(){
  char tmp[256];
  std::ifstream fpgold("stopping_power_p_gold.dat");
  for (int i = 0; i < 8; i++)
    fpgold.getline(tmp, 256);
  for (int i = 0; i < NPT_P_GOLD; i++){
    fpgold >> EK_P_GOLD[i] >> dEK_P_GOLD[i];
    EK_P_GOLD[i] = EK_P_GOLD[i] / 1000.0;//convert to GeV
    dEK_P_GOLD[i] = dEK_P_GOLD[i] * 19.32 / 1000.0;//GeV / cm
  }
  SP_P_GOLD.SetData(NPT_P_GOLD, EK_P_GOLD, dEK_P_GOLD);
  fpgold.close();
  std::ifstream fkgold("stopping_power_p_gold.dat");
  for (int i = 0; i < 8; i++)
    fkgold.getline(tmp, 256);
  for (int i = 0; i < NPT_K_GOLD; i++){
    fkgold >> EK_K_GOLD[i] >> dEK_K_GOLD[i];
    EK_K_GOLD[i] = EK_K_GOLD[i] / 1000.0;//convert to GeV
    dEK_K_GOLD[i] = dEK_K_GOLD[i] * 19.32 / 1000.0;//GeV / cm
  }
  SP_K_GOLD.SetData(NPT_K_GOLD, EK_K_GOLD, dEK_K_GOLD);
  fkgold.close();
  return 0;
}

double CheckBONUS(const TLorentzVector * P, const char * particle = "K+", const double z = 0.5){//Get the probability of the detection by BONUS
  if (MomentumCut(P, 0.0, 0.06))//momentum below 60 MeV
    return 0.0;
  if (!PolarAngleCut(P, 20.0, 160.0, "deg"))//angle outside BONUS coverage
    return 0.0;
  double Ek0 = P->E() - P->M();//Get kinematical energy
  double dEk = 0;
  if (strcmp(particle, "K+") == 0 || strcmp(particle, "K-") == 0)
    dEk = SP_K_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else if (strcmp(particle, "p") == 0)
    dEk = SP_P_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else {
    std::cerr << "Unknown particle option in CheckBONUS!" << std::endl;
    return 0;
  }
  double dz = target_length * (1.0 - z) / cos(P->Theta());
  if (dz < 0.0) dz = dz / (1.0 - z) * (-z);
  double dx = std::min(dz, target_radius / sin(P->Theta()));//travel distance in target
  double Ek = Ek0 - dEk * dx;//Kinematical energy outside target
  if (Ek < 0.0)//fail to get out off target
    return 0.0;
  double p = sqrt(pow(Ek+P->M(), 2) - P->M2());//Get momentum
  if (p < 0.06 || p > 0.25)//momentum out off BONUS coverage
    return 0.0;
  double wd = 1.0;
  return wd;
}

double CheckCLAS12FA(const TLorentzVector * P, const char * particle = "K+", const double z = 0.5){//Get the probability of the detection by CLAS12 forward detector
  if (MomentumCut(P, 0.0, 0.25))//momentum below 250 MeV
    return 0.0;
  if (!PolarAngleCut(P, 5.0, 35.0, "deg"))//angle outside CLAS12 FA coverage
    return 0.0;
  double Ek0 = P->E() - P->M();//Get kinematical energy
  double dEk = 0;
  if (strcmp(particle, "K+") == 0 || strcmp(particle, "K-") == 0)
    dEk = SP_K_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else if (strcmp(particle, "p") == 0)
    dEk = SP_P_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else {
    std::cerr << "Unknown particle option in CheckBONUS!" << std::endl;
    return 0;
  }
  double dz = target_length * (1.0 - z) / cos(P->Theta());
  if (dz < 0.0) dz = dz / (1.0 - z) * (-z);
  double dx = std::min(dz, target_radius / sin(P->Theta()));//travel distance in target
  double Ek = Ek0 - dEk * dx;//Kinematical energy outside target
  if (Ek < 0.0)//fail to get out off target
    return 0.0;
  double p = sqrt(pow(Ek+P->M(), 2) - P->M2());//Get momentum
  if (p < 0.25)//momentum out off BONUS coverage
    return 0.0;
  double wd = 1.0;
  if (strcmp(particle, "K+") == 0 || strcmp(particle, "K-") == 0){
    double th = P->Theta();
    double azi = P->Phi();
    TLorentzVector PP;
    PP.SetXYZM(p * sin(th) * cos(azi), p * sin(th) * sin(azi), p * cos(th), 0.493677);
    wd = LongLifetimeDecayFactor(&PP, 3.712, 2.0);
  }
  return wd;
}

double CheckCLAS12LA(const TLorentzVector * P, const char * particle = "K+", const double z = 0.5){//Get the probability of the detection by CLAS12 central detector
  if (MomentumCut(P, 0.0, 0.25))//momentum below 250 MeV
    return 0.0;
  if (!PolarAngleCut(P, 35.0, 125.0, "deg"))//angle outside CLAS12 FA coverage
    return 0.0;
  double Ek0 = P->E() - P->M();//Get kinematical energy
  double dEk = 0;
  if (strcmp(particle, "K+") == 0 || strcmp(particle, "K-") == 0)
    dEk = SP_K_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else if (strcmp(particle, "p") == 0)
    dEk = SP_P_GOLD.Eval(Ek0);//Get dE/dx in GeV/cm
  else {
    std::cerr << "Unknown particle option in CheckBONUS!" << std::endl;
    return 0;
  }
  double dz = target_length * (1.0 - z) / cos(P->Theta());
  if (dz < 0.0) dz = dz / (1.0 - z) * (-z);
  double dx = std::min(dz, target_radius / sin(P->Theta()));//travel distance in target
  double Ek = Ek0 - dEk * dx;//Kinematical energy outside target
  if (Ek < 0.0)//fail to get out off target
    return 0.0;
  double p = sqrt(pow(Ek+P->M(), 2) - P->M2());//Get momentum
  if (p < 0.25)//momentum out off BONUS coverage
    return 0.0;
  double wd = 1.0;
  if (strcmp(particle, "K+") == 0 || strcmp(particle, "K-") == 0){
    double th = P->Theta();
    double azi = P->Phi();
    TLorentzVector PP;
    PP.SetXYZM(p * sin(th) * cos(azi), p * sin(th) * sin(azi), p * cos(th), 0.493677);
    wd = LongLifetimeDecayFactor(&PP, 3.712, 2.0);
  }
  return wd;
}


#endif
