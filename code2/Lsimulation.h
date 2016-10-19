#ifndef _LSIMULATION_H_
#define _LSIMULATIOM_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

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
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"

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

double BremsstrahlungPhotonLumi(const double * k, const double * par){//Calculate the ratio of effective gamma A lumi to eN lumi
  double E0 = k[0];//Beam energy
  double kmin = k[1];//photon energy minimum
  double kmax = k[2];//photon energy maximum
  double fL = par[0];//radiation length fraction
  double A = par[1];//nucleon number
  double Y = 4.0 / 3.0 * log(kmax / kmin) - 4.0 / (3.0 * E0) * (kmax - kmin) + 1.0 / (2.0 * E0 * E0) * (kmax * kmax - kmin * kmin);
  double result = Y * fL / (2.0 * A);
  return result;
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

double dSigmaPhiProduction(const TLorentzVector * ki, const TLorentzVector * kf){//differential cross section of phi photoproduction in c.m. frame
  double par[5] = {0.224551, 2.00386, 3.75764, 1.38537, 0.909071};
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
  double result = a0 * (r0 + a2 * r2) / (1.0 + a2);//ds in unit mb
  return result / 389.379;//ds in unit GeV^-2
}

double dSigmaVirtualPhiProduction(const TLorentzVector * ki, const TLorentzVector * kf){//differential cross section of phi photoproduction in c.m. frame
  double par[5] = {0.224551, 2.00386, 3.75764, 1.38537, 0.909071};
  double x = GetTheta0(ki, kf);//theta in c.m. frame
  double M = kf[0].M() + kf[1].M();//threshold
  TLorentzVector Pout = kf[0] + kf[1];//total 4-momentum
  double sr = Pout.M();//energy in c.m. frame
  if (sr < M) return 0.0;//below threshold
  double b1 = par[3] * pow(sr * sr - M * M, par[4]);
  double a2 = par[2] * (sr - M);
  double a0 = 0.33;//constant total cross section in unit mb
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double result = a0 * (r0 + a2 * r2) / (1.0 + a2);//ds in unit mb
  return result / 389.379;//ds in unit GeV^-2
}

double dSigmaLambda1520Production(const TLorentzVector * ki, const TLorentzVector * kf){//differential cross section of Lambda 1520 + K photoproduction in c.m. frame
  double par[5] = {11.299, 4.60959, 0.835621, 0.54681, 1.827941};
  double x = GetTheta0(ki, kf);//theta in c.m. frame
  double M = kf[0].M() + kf[1].M();//threshold
  TLorentzVector Pout = kf[0] + kf[1];//total 4-momentum
  double sr = Pout.M();//energy in c.m. frame
  if (sr < M) return 0.0;//below threshold
  double b1 = par[3] * (sr - M + par[4]);
  double a2 = 0.25;
  double a0 = par[0] * (sr - M) * exp(-par[1] * pow(sr - M, par[2]));;//total cross section in unit mb
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double result = a0 * (r0 + a2 * r2) / (1.0 + a2);//ds in unit mb
  return result / 389.379;//ds in unit GeV^-2
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
  double result1 = (A0 + A2 * Q * Q) * exp(-B2 * Q * Q);//without cutoff
  //double result2 = 0.114843 * exp(-Q * Q / (2.0 * pow(0.2, 2)));//with cutoff
  return result1;//amplitude in unit GeV^-1/2;
}

double ProbabilityBoundState(const double * mom, const double * par = 0){//Probability density w.r.t Q
  double result = pow(AmplitudeBoundState(mom, par), 2);
  return result;//in unit GeV^-1
}

int DetectorResolutionSmear(TLorentzVector * P, double * res){//Smearing the three momentum with gaussian, fixed mass
  double M = P->M();
  double p = P->P() * (1.0 + gRandom->Gaus(0.0, res[0]));//smear p
  if (p < 0.0) p = 0.0;//Set negative p to 0
  double theta = P->Theta() + gRandom->Gaus(0.0, res[1]);//smear theta
  double phi = P->Phi() + gRandom->Gaus(0.0, res[2]);//smear phi
  P[0].SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), M);//Set smeared three momentum with fixed mass
  //P[0].SetRho(p); P[0].SetTheta(theta); P[0].SetPhi(phi);
  return 0;
} 

/******* Functions for random sampling *******/
TF1 TF_fq("fq", Bremsstrahlung, 1.1, 1.8, 1);//set photon energy distri
TF1 TF_fp("fp", CarbonMomentum, 0.0, 1.0, 0);//set bound N momentum distri
TF1 TF_fE("fE", CarbonEnergy, 0.0, 0.1, 0);//set bound N missing energy distri
TF1 TF_BWPhi("BWPhi", BreitWigner, 0.819, 1.219, 2);//set phi mass distri
TF1 TF_BWL1520("BWL1520", BreitWigner, 1.43195, 1.600, 2);//set Lambda1520 mass distri
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
  TF_BWL1520.SetParameter(0, 1.5195);
  TF_BWL1520.SetParameter(1, 0.0156);
  TF_BWL1520.SetNpx(2000);
  TF_BWd.SetParameter(0, 1.950027);
  TF_BWd.SetParameter(1, 0.002118);
  TF_BWd.SetNpx(4000);
  return 0;
}

double GenerateBremsstrahlungPhoton(TLorentzVector * q, const double * k){//Generate a photon from Bremsstrahlung
  double E0 = TF_fq.GetRandom(k[0], k[1]);
  q[0].SetXYZT(0.0, 0.0, E0, E0);
  return 1.0;
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
  const double Mphi = TF_BWPhi.GetRandom();//Get a phi meson mass
  TLorentzVector Pout = q + p1;//Total 4-momentum
  const double s = Pout.M2();//mandelstem variable s
  if (s <= pow(Mp + Mphi, 2)){
    weight[0] = 0.0;
    return 0.0;//below the threshold
  }
  double mass[2] = {Mphi, Mp};//Set final particle masses
  Lphase.SetDecay(Pout, 2, mass);//Set kinematics
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get phi 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set phi 4-momentum
  tm = Lphase.GetDecay(1);//Get recoil proton 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set recoil proton 4-momentum
  double ds = dSigmaPhiProduction(ki, kf);//ds / dcos theta
  weight[0] = ds * 2;
  return weight[0];//cross section weighting factor in unit GeV^-2
}

double GenerateVirtualPhiProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate a phi meson photoproduction event
  TLorentzVector q = ki[0];//4-momentum of incoming photon
  TLorentzVector p1 = ki[1];//4-momentum of incoming bound nucleon
  const double Mp = 0.938272;
  const double Mphi = TF_BWPhi.GetRandom();//off-shell supression
  TLorentzVector Pout = q + p1;//Total 4-momentum
  const double s = Pout.M2();//mandelstem variable s
  if (s <= pow(Mp + Mphi, 2)){
    weight[0] = 0.0;
    return 0.0;//below the threshold
  }
  double mass[2] = {Mphi, Mp};//Set final particle masses
  Lphase.SetDecay(Pout, 2, mass);//Set kinematics
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get phi 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set phi 4-momentum
  tm = Lphase.GetDecay(1);//Get recoil proton 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set recoil proton 4-momentum
  double ds = dSigmaVirtualPhiProduction(ki, kf);//ds / dcos theta
  weight[0] = ds * 2;
  return weight[0];//cross section weighting factor in unit GeV^-2
}

double GenerateLambda1520Production(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate a Lambda 1520 photoproduction event
  TLorentzVector q = ki[0];//4-momentum of incoming photon
  TLorentzVector p1 = ki[1];//4-momentum of incoming bound nucleon
  const double Mkaon = 0.493677;
  const double MLambda = TF_BWL1520.GetRandom();
  TLorentzVector Pout = q + p1;//Total 4-momentum
  const double s = Pout.M2();//mandelstem variable s
  if (s <= pow(Mkaon + MLambda, 2)){
    weight[0] = 0.0;
    return 0.0;//below the threshold
  }
  double mass[2] = {Mkaon, MLambda};//Set final particle masses
  Lphase.SetDecay(Pout, 2, mass);//Set kinematics
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get phi 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set kaon 4-momentum
  tm = Lphase.GetDecay(1);//Get recoil proton 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set Lambda1520 4-momentum
  double ds = dSigmaLambda1520Production(ki, kf);//ds / dcos theta
  weight[0] = ds * 2.0;
  return weight[0];//cross section weighting factor of particular kaon and Lambda1520, in unit GeV^-2
}

double GenerateBoundStateFormation(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator a bound state event with 4-momentum matching
  TLorentzVector Pout = ki[0] + ki[1];//Total 4-momentum
  double Q = GetRelativeMomentum(ki);//Get relative momentum in c.m.
  if (Q <= 0.0){
    weight[0] = 0.0;
    return 0.0;
  }
  weight[0] = ProbabilityBoundState(&Q);//Momentum matching probability
  kf[0] = Pout;//Set bound state momentum
  return weight[0];//in unit GeV^-2
}

double GenerateBoundStateProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator an event of bound state photoproduction
  TLorentzVector q = ki[0];//4-momentum of photon
  TLorentzVector p1, p2;
  GenerateNucleonInCarbon(&p1);//Generate first bound N in carbon
  GenerateNucleonInCarbon(&p2);//Generate second bound N in carbon
  TLorentzVector ki1[2] = {q, p1};//Set incoming particles in step 1
  TLorentzVector kf1[2];//outgoing particles in step 1
  double weight1;
  GenerateVirtualPhiProduction(ki1, kf1, &weight1);//Generate phi production
  if (weight1 < 1.0e-85){
    weight[0] = 0.0;
    return 0.0;
  }
  TLorentzVector k = kf1[0];//Get produced phi
  kf[0] = kf1[1];//Get recoil proton
  TLorentzVector ki2[2] = {k, p2};//Set incoming particles in step 2
  double weight2;
  GenerateBoundStateFormation(ki2, &kf[1], &weight2);//Generate bound state
  double weight3;//Energy matching weight
  if (kf[1].M() <= TF_BWd.GetMinimum() || kf[1].M() >= TF_BWd.GetMaximum())
    weight3 = 0.0;
  else
    weight3 = TF_BWd.Eval(kf[1].M()) * kf[1].E() / kf[1].M() / 194.844;//!!! mass width par dependent normalization
  //double weight3 = TF_BWd.Eval(kf[1].M()) * kf[1].E() / kf[1].M() / TF_BWd.Integral(TF_BWd.GetMinimum(), TF_BWd.GetMaximum());//normalized weight
  double Lt = 0.01864;//in unit GeV^2 for C-12(11)
  weight[0] = weight1 * weight2 * weight3 * Lt;
  if (isnan(weight[0])) std::cout << "Warning!" << std::endl;
  return weight[0];
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

double GenerateKKProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){
  /* need update */
//Generate a direct KK photoproduction event
  TLorentzVector q = ki[0];//4-momentum of incoming photon
  TLorentzVector p1 = ki[1];//4-momentum of incoming bound nucleon
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//K meson mass
  TLorentzVector Pout = q + p1;//Total 4-momentum
  double s = Pout.M2();//Total invariant mass
  if (s <= pow(Mp + 2.0 * Mkaon, 2)){
    weight[0] = 0.0;
    return 0.0;
  }
  double mass[3] = {Mp, Mkaon, Mkaon};//Set final particle masses
  double weight1 = Decay3(&Pout, kf, mass);//phase space weight
  TLorentzVector Ptoy = kf[1] + kf[2];;//virtual intermediate
  TLorentzVector Pf[2] = {Ptoy, kf[0]};//Set intermediate state
  double ds = dSigmaPhiProduction(ki, Pf);
  double weight2 = ds * 2;//cross section weight
  double Br = 0.489;//Branch ratio to K+K-
  weight[0] = weight1 * weight2 * Br;
  return weight[0];
}

double GenerateEventNNKKwithBoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NNKK with bound state production
  const double Mp = 0.938272;
  const double Mkaon = 0.493677;
  TLorentzVector kf1[2];
  double weight1;
  GenerateBoundStateProduction(ki, kf1, &weight1);//Generate a bound state
  kf[0] = kf1[0];//recoil proton
  TLorentzVector Pd = kf1[1];//4-momentum of bound state
  if (weight1 < 1.0e-85){//bound state not produced
    weight[0] = 0.0;
    return 0.0;
  }
  if (Pd.M() <= Mp + 2.0 * Mkaon){//below NKK threshold
    weight[0] = 0.0;
    return 0.0;
  }
  double mass[3] = {Mp, Mkaon, Mkaon};
  TLorentzVector kf2[3];
  double weight2 = Decay3(&Pd, kf2, mass);//Generate bound state decay
  kf[1] = kf2[0];//decayed proton
  kf[2] = kf2[1];//decayed K+
  kf[3] = kf2[2];//decayed K-
  double Br = 0.919 * 0.489;//Branching ratio to NK+K- channel
  weight[0] = weight1 * weight2 * Br;//in unit GeV^-2
  return weight[0];
}

double GenerateEventNKKwithPhi(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NKK with phi production
  TLorentzVector p1;
  GenerateNucleonInCarbon(&p1);//Get a nucleon in carbon
  TLorentzVector ki1[2] = {ki[0], p1};//Set incoming particles for phi production
  TLorentzVector kf1[2];
  double weight1;
  GeneratePhiProduction(ki1, kf1, &weight1);
  if (weight1 < 1.0e-85){//phi meson not produced
    weight[0] = 0.0;
    return 0.0;
  }
  const double Mkaon = 0.493677;
  if (kf1[0].M() <= 2.0 * Mkaon){//below KK threshold
    weight[0] = 0.0;
    return 0.0;
  }  
  kf[0] = kf1[1];//recoil proton
  TLorentzVector kphi = kf1[0];//produced phi meson
  double mass[2] = {Mkaon, Mkaon};
  TLorentzVector kf2[2];
  Decay2(&kphi, kf2, mass);
  kf[1] = kf2[0];//decayed K+
  kf[2] = kf2[1];//decayed K-
  double Br = 0.489;//Branch ratio of phi -> K+K- channel
  weight[0] = weight1 * Br;//in unit GeV^-2
  return weight[0];
}

double GenerateEventNKKwithLambda1520(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NKK with Lambda1520 production
  TLorentzVector p1;
  GenerateNucleonInCarbon(&p1);//Get a nucleon in carbon
  TLorentzVector ki1[2] = {ki[0], p1};//Set incoming particles for Lambda1520 production
  TLorentzVector kf1[2];
  double weight1;
  GenerateLambda1520Production(ki1, kf1, &weight1);
  kf[1] = kf1[0];//produced K+
  TLorentzVector kLambda = kf1[1];//produced Lambda1520
  double mass[2] = {0.938272, 0.493677};
  TLorentzVector kf2[2];
  Decay2(&kLambda, kf2, mass);
  kf[0] = kf2[0];//decayed p
  kf[2] = kf2[1];//decayed K-
  double Br = 0.45 / 2.0;//Branch ratio of Lambda1520 -> NK channel divided by 2 to select pK- only
  weight[0] = weight1 * Br;
  return weight[0];
}

double GenerateEventNKKwithout(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NKK direct production
  TLorentzVector p1;
  GenerateNucleonInCarbon(&p1);//Get a nucleon in carbon
  TLorentzVector ki1[2] = {ki[0], p1};//Set incoming particles
  double weight1;
  GenerateKKProduction(ki1, kf, &weight1);
  weight[0] = weight1;
  return weight[0];
}



#endif
