#ifndef _LSIMULATION_H_
#define _LSIMULATIOM_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "TROOT.h"
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

double Bremsstrahlung(const double * k, const double * par){//non-normalized
  //E0: electron beam energy, k: photon energy
  double E0 = par[0];
  double y = k[0] / E0;
  if (y < 0.01) {
    std::cerr << "Out of range in Bremsstrahlung!" << std::endl;
    return -1.0;
  }
  double result = 1.0 / (y * E0) * (4.0 / 3.0 - 4.0 / 3.0 * y + y * y);
  return result;//non-normalized probability density
}

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
  return result;//non-normalized momentum distribution
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
  return result;//non-normalized missing energy distribution
}

double SigmaPhiProduction(const double * mom, const double * par){//total cross section of phi meson photoproduction from a bound nucleon
  double scalar = mom[0];//scalar product of income p1 and p2
  double Q = mom[1];//relative momentum of produced phi meson in c.m. frame
  double Mphi = par[0];//phi meson mass
  const double Amp = 0.59405;//invariant amplitude
  const double Mp = 0.938272;//proton mass
  double result = Amp * Amp / (16.0 * M_PI * scalar) * Q / (sqrt(Q*Q+Mphi*Mphi) + sqrt(Q*Q+Mp*Mp));
  return result;//cross section in unit GeV^-2
}

double AmplitudeBoundState(const double * mom, const double * par = 0){//Amplitude of phi N to bound state
  double Q = mom[0];//relative momentum of N phi system in unit GeV
  if (Q < 0.0) {
    std::cerr << "Unphysical relative momentum in AmplitudeBoundState!" << std::endl;
    return -1.0;
  }
  double A0 = 0.112356;//Fit to FQ0.dat with poly*gaussian
  double A2 = 0.609078;
  double B2 = 2.75841;
  double result = (A0 + A2 * Q * Q) * exp(-B2 * Q * Q);
  return result;//amplitude in unit GeV^-1/2;
}

double ProbabilityBoundState(const double * mom, const double * par = 0){//Probability density w.r.t Q
  double result = pow(AmplitudeBoundState(mom, par), 2);
  return result;//in unit GeV^-1
}

/******* Functions for random sampling *******/
TF1 TF_fq("fq", Bremsstrahlung, 1.2, 1.8, 1);//set photon energy distri
TF1 TF_fp("fp", CarbonMomentum, 0.0, 0.5, 0);//set bound N momentum distri
TF1 TF_fE("fE", CarbonEnergy, 0.0, 0.1, 0);//set bound N missing energy distri
TF1 TF_BWPhi("BWPhi", BreitWigner, 0.819, 1.219, 2);//set phi mass distri
TF1 TF_BWd("BWd", BreitWigner, 1.750, 2.150, 2);//set bound state mass distri
TGenPhaseSpace Lphase;//A global variable to generate final state
/*** End of Functions for random sampling  ***/

int SetFunctions(){
  TF_fq.SetParameter(0, 11.0);
  TF_fq.SetNpx(500);
  TF_fp.SetNpx(500);
  TF_fE.SetNpx(500);
  TF_BWPhi.SetParameter(0, 1.019455);
  TF_BWPhi.SetParameter(1, 0.00426);
  TF_BWPhi.SetNpx(2000);
  TF_BWd.SetParameter(0, 1.950027);
  TF_BWd.SetParameter(1, 0.001789);
  TF_BWd.SetNpx(4000);
  return 0;
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

double GeneratePhiProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate a phi meson photoproduction event
  TLorentzVector q = ki[0];//4-momentum of incoming photon
  TLorentzVector p1 = ki[1];//4-momentum of incoming bound nucleon
  const double Mp = 0.938272;
  const double Mphi = TF_BWPhi.GetRandom();
  TLorentzVector Pout = q + p1;//Total 4-momentum
  const double s = Pout.M2();//mandelstem variable s
  if (s <= pow(Mp + Mphi, 2)){
    return 0.0;//below the threshold
  }
  double mom[2];
  mom[0] = q * p1;//scalar product of incoming 2 particles
  mom[1] = sqrt( (s - pow(Mp + Mphi, 2)) * (s - pow(Mphi - Mp, 2)) / (4.0 * s) );//Relative momentum Q
  double sigma = SigmaPhiProduction(mom, &Mphi);//total cross section
  double mass[2] = {Mphi, Mp};//Set final particle masses
  Lphase.SetDecay(Pout, 2, mass);//Set kinematics
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get phi 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set phi 4-momentum
  tm = Lphase.GetDecay(1);//Get recoil proton 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set recoil proton 4-momentum
  weight[0] = sigma;
  return sigma;//return total cross section with particular p1 and Mphi, in unit GeV^-2
}

double GenerateBoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator a bound state event
  TLorentzVector k = ki[0];//4-momentum of phi
  TLorentzVector p2 = ki[1];//4-momentum of bound Nucleon in carbon
  TLorentzVector Pout = k + p2;//Total 4-momentum
  double Md = Pout.M();
  const double Mmin = TF_BWd.GetXmin();
  const double Mmax = TF_BWd.GetXmax();
  if (Md < Mmin || Md >= Mmax)
    weight[0] = 0.0;
  else
    weight[0] = TF_BWd.Eval(Md) / TF_BWd.Integral(Mmin, Mmax);
  double s = Pout.M2();
  double Q = sqrt( (s - pow(k.M() + p2.M(), 2)) * (s - pow(k.M() - p2.M(), 2)) / (4.0 * s));
  weight[1] = ProbabilityBoundState(&Q);
  kf[0] = Pout;
  return weight[0] * weight[1];//in unit GeV^-2
}

double GenerateEventBoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator an event of bound photoproduction
  TLorentzVector q = ki[0];//4-momentum of photon
  TLorentzVector p1, p2;
  GenerateNucleonInCarbon(&p1);//Generate first bound N in carbon
  GenerateNucleonInCarbon(&p2);//Generate second bound N in carbon
  TLorentzVector ki1[2] = {q, p1};//Set incoming particles in step 1
  TLorentzVector kf1[2];//outgoing particles in step 1
  double weight1;
  GeneratePhiProduction(ki1, kf1, &weight1);//Generate step 1
  TLorentzVector k = kf1[0];//Get produced phi
  kf[0] = kf1[1];//Get recoil proton
  TLorentzVector ki2[2] = {k, p2};//Set incoming particles in step 2
  TLorentzVector Pd;//4-momentum of bound state
  double weight2[2];
  GenerateBoundState(ki2, &Pd, weight2);//Generate step 2
  kf[1] = Pd;//Get bound state 4-momentum
  weight[0] = weight1 * weight2[0] * weight2[1];
  weight[1] = weight1;
  weight[2] = weight2[0];
  weight[3] = weight2[1];
  return weight[0];
}
  
  







#endif
