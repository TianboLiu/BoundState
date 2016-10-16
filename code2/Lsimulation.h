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
  return result / 1.89142;//Normalized momentum distribution
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

double SigmaPhiProduction(const double * mom, const double * par){//total cross section of phi meson photoproduction from a bound nucleon
  double scalar = mom[0];//scalar product of income p1 and p2
  double Q = mom[1];//relative momentum of produced phi meson in c.m. frame
  double Mphi = par[0];//phi meson mass
  const double Amp = 0.594054;//invariant amplitude
  const double Mp = 0.938272;//proton mass
  double result = Amp * Amp / (16.0 * M_PI * scalar) * Q / (sqrt(Q*Q+Mphi*Mphi) + sqrt(Q*Q+Mp*Mp));
  return 0.17/0.389379e3;//cross section in unit GeV^-2
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
  double cth = GetTheta0(ki, kf);
  double result = exp(cth);
  return result;
}

double SigmaLambda1520Production(const double * mom, const double * par){//total cross section of Lambda 1520 + K photoproduction
  double scalar = mom[0];//scalar product of income p1 and p2
  double Q = mom[1];//relative momentum of produced K meson in c.m. frame
  double MLambda = par[0];//Lambda 1520 mass
  const double Amp = 0.801197;//invariant amplitude
  const double Mkaon = 0.493677;//kaon mass
  double result = Amp * Amp / (16.0 * M_PI * scalar) * Q / (sqrt(Q*Q+MLambda*MLambda) + sqrt(Q*Q+Mkaon*Mkaon));
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
  double result1 = (A0 + A2 * Q * Q) * exp(-B2 * Q * Q);
  double result2 = 0.114843 * exp(-100.15 * Q * Q);
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
TF1 TF_fp("fp", CarbonMomentum, 0.0, 0.5, 0);//set bound N momentum distri
TF1 TF_fE("fE", CarbonEnergy, 0.0, 0.1, 0);//set bound N missing energy distri
TF1 TF_BWPhi("BWPhi", BreitWigner, 0.987354, 1.219, 2);//set phi mass distri
TF1 TF_BWL1520("BWL1520", BreitWigner, 1.43195, 1.600, 2);//set Lambda1520 mass distri
TF1 TF_BWd("BWd", BreitWigner, 1.925626, 2.150, 2);//set bound state mass distri
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
  const double Mphi = TF_BWPhi.GetRandom();
  TLorentzVector Pout = q + p1;//Total 4-momentum
  const double s = Pout.M2();//mandelstem variable s
  if (s <= pow(Mp + Mphi, 2)){
    weight[0] = 0.0;
    return 0.0;//below the threshold
  }
  double mom[2];
  mom[0] = std::abs(q * p1);//scalar product of incoming 2 particles
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
  sigma = dSigmaPhiProduction(ki, kf) * 2.0;//test
  weight[0] = sigma;
  return sigma;//return total cross section with particular p1 and Mphi, in unit GeV^-2
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
  double mom[2];
  mom[0] = q * p1;//scalar product of incoming 2 particles
  mom[1] = sqrt( (s - pow(Mkaon + MLambda, 2)) * (s - pow(Mkaon - MLambda, 2)) / (4.0 * s) );//Relative momentum Q
  double sigma = SigmaLambda1520Production(mom, &MLambda);//total cross section
  double mass[2] = {Mkaon, MLambda};//Set final particle masses
  Lphase.SetDecay(Pout, 2, mass);//Set kinematics
  Lphase.Generate();
  TLorentzVector * tm;//tmp 4-momentum pointer
  tm = Lphase.GetDecay(0);//Get phi 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set kaon 4-momentum
  tm = Lphase.GetDecay(1);//Get recoil proton 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set Lambda1520 4-momentum
  weight[0] = sigma;
  return sigma;//return total cross section with particular kaon and Lambda1520, in unit GeV^-2
}

double GenerateKKProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate a direct KK photoproduction event
  TLorentzVector q = ki[0];//4-momentum of incoming photon
  TLorentzVector p1 = ki[1];//4-momentum of incoming bound nucleon
  double scalar = std::abs(q * p1);//scalar product of incoming paricles
  const double Mp = 0.938272;//proton mass
  const double Mkaon = 0.493677;//K meson mass
  TLorentzVector Pout = q + p1;//Total 4-momentum
  double s = Pout.M2();//Total invariant mass
  if (s <= pow(Mp + 2.0 * Mkaon, 2)){
    weight[0] = 0.0;
    weight[1] = 0.0;
    weight[2] = 0.0;
    weight[3] = 0.0;
    return 0.0;
  }
  double mass[3] = {Mp, Mkaon, Mkaon};//Set final particle masses
  Lphase.SetDecay(Pout, 3, mass);//
  Lphase.Generate();
  TLorentzVector * tm;
  tm = Lphase.GetDecay(0);//Get p 4-momentum
  kf[0].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set p 4-momentum
  tm = Lphase.GetDecay(1);//Get K+ 4-momentum
  kf[1].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set p 4-momentum
  tm = Lphase.GetDecay(2);//Get K- 4-momentum
  kf[2].SetXYZT(tm->X(), tm->Y(), tm->Z(), tm->T());//Set p 4-momentum
  TLorentzVector Ptoy;//virtual intermediate
  double Mtoy;
  double mom[2] = {scalar, 0.0};
  Ptoy = kf[1] + kf[2];//4-momentum of (K+K-)
  Mtoy = Ptoy.M();//virtual "phi"(K+K-) mass
  mom[1] = sqrt( (s - pow(Mp + Mtoy, 2)) * (s - pow(Mp - Mtoy, 2)) / (4.0 * s) );//Relative momentum Q
  weight[1] = SigmaPhiProduction(mom, &Mtoy);
  Ptoy = kf[0] + kf[1];//4-momentum of (pK+)
  Mtoy = Ptoy.M();//virtual (pK+)
  mom[1] = sqrt( (s - pow(Mkaon + Mtoy, 2)) * (s - pow(Mkaon - Mtoy, 2)) / (4.0 * s) );//Relative momentum Q
  weight[2] = SigmaLambda1520Production(mom, &Mtoy);
  Ptoy = kf[0] + kf[2];//4-momentum of (pK-)
  Mtoy = Ptoy.M();//virtual (pK-)
  mom[1] = sqrt( (s - pow(Mkaon + Mtoy, 2)) * (s - pow(Mkaon - Mtoy, 2)) / (4.0 * s) );//Relative momentum Q
  weight[3] = SigmaLambda1520Production(mom, &Mtoy);
  weight[0] = (weight[1] + weight[2] + weight[3]) / 3.0;
  return weight[0];
}

double GenerateBoundStateFormation(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator a bound state event with 4-momentum matching
  TLorentzVector k = ki[0];//4-momentum of phi
  TLorentzVector p2 = ki[1];//4-momentum of bound Nucleon in carbon
  TLorentzVector Pout = k + p2;//Total 4-momentum
  double s = Pout.M2();
  double Q = sqrt( (s - pow(k.M() + p2.M(), 2)) * (s - pow(k.M() - p2.M(), 2)) / (4.0 * s));
  weight[0] = ProbabilityBoundState(&Q);
  kf[0] = Pout;
  return weight[0];//in unit GeV^-2
}

double GenerateBoundStateProduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generator an event of bound photoproduction
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
  double Md = TF_BWd.GetRandom();//Get bound state mass according to B-W
  TLorentzVector Pout = k + p2;//Get total momentum
  Pout.SetE(sqrt(Md * Md + Pout.P() * Pout.P()));//Set on-shell energy
  const double Mp = 0.938272;
  p2.SetE(Pout.E() - k.E());//Set required nucleon energy
  double Em = sqrt(Mp * Mp + p2.P() * p2.P()) - p2.E();//Calculate missing energy
  double weight2;//energy match weight
  if (Em < TF_fE.GetXmin() || Em >= TF_fE.GetXmax()){
    weight2 = 0.0;
    weight[0] = 0.0;
    return 0.0;
  }
  else
    weight2 = 1.0;//TF_fE.Eval(Em);//Get the probability of find such nucleon
  TLorentzVector ki2[2] = {k, p2};//Set incoming particles in step 2
  TLorentzVector Pd;//4-momentum of bound state
  double weight3;//momentum match weight
  GenerateBoundStateFormation(ki2, &Pd, &weight3);//Generate step 2
  kf[1] = Pd;//Get bound state 4-momentum
  weight[0] = weight1 * weight2 * weight3;
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

double GenerateEventNNKKwithBoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NNKK with bound state production
  TLorentzVector kf1[2];
  double weight1[4];
  GenerateBoundStateProduction(ki, kf1, weight1);
  kf[0] = kf1[0];//recoil proton
  TLorentzVector Pd = kf1[1];//4-momentum of bound state
  if (weight1[0] == 0.0) {
    weight[0] = 0.0;
    return 0.0;
  }
  double mass[3] = {0.938272, 0.493677, 0.493677};
  TLorentzVector kf2[3];
  double weight2 = Decay3(&Pd, kf2, mass);
  //std::cout << weight2 << std::endl;
  kf[1] = kf2[0];//decayed proton
  kf[2] = kf2[1];//decayed K+
  kf[3] = kf2[2];//decayed K-
  double Br = 0.919 * 0.489;//Branching ratio of NKK channel
  weight[0] = weight1[0] * weight2 * Br;
  return weight[0];
}

double GenerateEventNKKwithPhi(const TLorentzVector * ki, TLorentzVector * kf, double * weight){//Generate an event of NKK with phi production
  TLorentzVector p1;
  GenerateNucleonInCarbon(&p1);//Get a nucleon in carbon
  TLorentzVector ki1[2] = {ki[0], p1};//Set incoming particles for phi production
  TLorentzVector kf1[2];
  double weight1;
  GeneratePhiProduction(ki1, kf1, &weight1);
  kf[0] = kf1[1];//recoil proton
  TLorentzVector kphi = kf1[0];//produced phi meson
  double mass[2] = {0.493677, 0.493677};
  TLorentzVector kf2[2];
  Decay2(&kphi, kf2, mass);
  kf[1] = kf2[0];//decayed K+
  kf[2] = kf2[1];//decayed K-
  double Br = 0.489;//Branch ratio of phi -> K+K- channel
  weight[0] = weight1 * Br;
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
  double weight1[4];
  GenerateKKProduction(ki1, kf, weight1);
  weight[0] = weight1[1] * 0.489;
  //weight[0] = (weight1[1] * 0.489 + weight1[2] * 0.45 + weight1[2] * 0.45) / 3.0;
  //weight[0] = weight[0] * 0.5;
  return weight[0];
}



#endif
