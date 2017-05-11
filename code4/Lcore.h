#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "Math/Functor.h"
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
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"

#include "Lparticle.h"

using namespace std;

const double Mp = PARTICLE::proton.M();

namespace MODEL{
  
  const double Mass = 1.950;//bound state mass
  const double Eb = Mp + PARTICLE::phi.M() - Mass;//binding energy
  const double Width = 0.004;//bound state width
  const double FractionNphi = 1.0;//Nphi component fraction
  const double BrNphi = 0.951;//branch ratio of decay from NJpsi component

  ROOT::Math::Interpolator Ur_INTER(161, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator Veff_INTER(60, ROOT::Math::Interpolation::kCSPLINE);

  int SetVeff(){//
    double x[60], y[60];
    ifstream infile("wave/Veff.dat");
    for (int i = 0; i < 60; i++)
      infile >> x[i] >> y[i];
    Veff_INTER.SetData(60, x, y);
    infile.close();
    return 0;
  }

  int SetUr(){//
    double x[161], y[161];
    ifstream infile("wave/wf.dat");
    for (int i = 0; i < 161; i++)
      infile >> x[i] >> y[i];
    Ur_INTER.SetData(161, x, y);
    infile.close();
    return 0;
  }

  double Veff(const double r){//effective potential
    double rfm = r * Phys::hbar;//convert GeV^-1 to fm
    if (rfm < 0 || rfm > 5.0)
      return 0;
    return Veff_INTER.Eval(rfm) / 1000.0;//in GeV
  }

  double Ur(const double r){//radial wave function
    double N = 26.568;//Normalization factor
    double rfm = r * Phys::hbar;//convert GeV^-1 to fm
    if (rfm < 0 || rfm > 7.95)
      return 0;
    return N * Ur_INTER.Eval(rfm);//in GeV^1/2
  }

  double FQ_integrand(const double r, void * par){
    double * k = (double *) par;
    if (k[0] == 0)
      return -r * Veff(r) * Ur(r);
    return -sin(k[0] * r) / k[0] * Veff(r) * Ur(r);
  }

  double FQk(double k){
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-4);
    ig.SetFunction(&FQ_integrand, &k);
    double result = ig.Integral(0.0, 30.0);
    return result / (4.0 * M_PI);
  }

  int CalculateFQ(){
    FILE * fp = fopen("wave/FQ0.dat", "w");
    double k, fk;
    for (int i = 0; i < 500; i++){
      k = i * 0.004;
      fk = FQk(k);
      cout << k << "   " << fk << endl;
      fprintf(fp, "%.6E\t%.6E\n", k, fk);
    }
    fclose(fp);
    return 0;
  }
      
  ROOT::Math::Interpolator FQ_INTER(500, ROOT::Math::Interpolation::kCSPLINE);
  int SetFQ(){
    ifstream infile("wave/FQ0.dat");
    double x[500], y[500];
    for (int i = 0; i < 500; i++)
      infile >> x[i] >> y[i];
    FQ_INTER.SetData(500, x, y);
    infile.close();
    return 0;
  }

  double FQ(const double k){
    if (k < 0.0 || k > 1.85)
      return 0;
    return FQ_INTER.Eval(k);
  }
  
  double BreitWigner(const double * E, const double * par){
    return 1.0 / (pow(E[0] * E[0] - Mass * Mass, 2) + E[0] * E[0] * Width * Width) / 96.7216;
  }

  TF1 TF_fMass("fM", BreitWigner, Mass - 5.0 * Width, Mass + 5.0 * Width, 0);

  int SetMODEL(){
    SetVeff();
    SetUr();
    SetFQ();
    TF_fMass.SetNpx(1000);
    return 0;
  }

}


namespace GOLD{

  const double NA = 197.0;
  const double ProtonDensity = 79.0 / (4.0 * M_PI * pow(7.3, 3) / 3.0) * pow(Phys::hbar, 3);//GeV^3

  double fMomentum(const double * p0, const double * par = 0){//non-normalized
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
  
  double fEnergy(const double * E0, const double * par = 0){
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

  TF1 TF_fMomentum("fp", fMomentum, 0.0, 1.0, 0);
  TF1 TF_fEnergy("fE", fEnergy, 0.0, 0.3, 0);
 
  int SetGOLD(){
    TF_fMomentum.SetNpx(1000);
    TF_fEnergy.SetNpx(1000);
    gRandom->SetSeed(0);
    return 0;
  }

}



namespace GENERATE{

  TRandom3 random(0);
  TGenPhaseSpace GenPhase;
  double Weight = 0.0;
  
  int NucleonGold(TLorentzVector * P){
    double p = GOLD::TF_fMomentum.GetRandom();
    double cth = random.Uniform(-1.0, 1.0);
    double phi = random.Uniform(-M_PI, M_PI);
    double dE = GOLD::TF_fEnergy.GetRandom();
    P->SetXYZT(p * sqrt(1.0 - cth * cth) * cos(phi), p * sqrt(1.0 - cth * cth) * sin(phi), p * cth, sqrt(p * p + Mp * Mp) - dE);
    return 0;
  }

  double dSigmaPhi(const double W0, const double W, const double cth){
    if (W <= W0) return 0;
    double par[5] = {0.232612, 1.95038, 4.02454, 1.52884, 0.525636};
    double b1 = par[3] * pow(W * W - W0 * W0, par[4]);
    double a2 = par[2] * (W - W0);
    double a0 = par[0] * atan(par[1] * par[1] * (W * W - W0 * W0));
    double r0 = exp(b1 * cth) * b1 / (2.0 * sinh(b1));
    double r2 = cth * cth * exp(b1 * cth) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
    double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
    return ds;//GeV^-2
  }

  double PhiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma, N; kf: phi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    const double s = Pout.M2();//c.m. energy square
    const double MK = PARTICLE::K.M();
    const double Mphi = PARTICLE::phi.RandomM(MK * 2.0 + 1.0e-20, 1.1);//random phi meson mass 
    if (s <= pow(Mphi + Mp, 2)){
      weight[0] = 0;
      return 0;
    }
    double mass[2] = {Mphi, Mp};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);
    kf[1] = *GenPhase.GetDecay(1);
    TLorentzVector ka = ki[0];
    TLorentzVector kc = kf[0];
    ka.Boost(-Pout.BoostVector());
    kc.Boost(-Pout.BoostVector());
    double cth = cos(ka.Angle(kc.Vect()));
    weight[0] = dSigmaPhi(Mphi + Mp, Pout.M(), cth) * 4.0 * M_PI;
    return weight[0];//GeV^-2
  }

  double PhiPhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: phi, N'
    TLorentzVector ki1[2];
    ki1[0] = ki[0];//photon
    NucleonGold(&ki1[1]);//off-shell nucleon
    PhiPhotoproduction(ki1, kf, weight);
    weight[0] *= GOLD::NA;
    return weight[0];//GeV^-2 
  }
  
  double BoundStateFormationGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){//
    //ki: phi; kf: d
    const double Md = MODEL::TF_fMass.GetRandom();//bound state mass
    const double dE = GOLD::TF_fEnergy.GetRandom();//missing energy
    const double Et = ki[0].E() - dE;
    const double k = ki[0].P();//phi momentum
    const double MM = Md * Md + k * k - Et * Et - Mp * Mp;
    const double cth = random.Uniform(-1.0, 1.0);
    const double phi = random.Uniform(-M_PI, M_PI);
    const double a = Et * Et - k * k * cth * cth;
    const double b = -MM * k * cth;
    const double c = Et * Et * Mp * Mp - MM * MM / 4.0;
    const double DD = b * b - 4.0 * a * c;
    if (DD < 0){
      weight[0] = 0;
      return 0;
    }
    double p2;
    if (a * c < 0)
      p2 = (-b + sqrt(DD)) / (2.0 * a);
    else
      p2 = (-b - sqrt(DD)) / (2.0 * a);
    if (p2 < 0){
      weight[0] = 0;
      return 0;
    }
    weight[0] = GOLD::fMomentum(&p2);
    TLorentzVector P2;
    P2.SetXYZT(p2 * sqrt(1.0 - cth * cth) * cos(phi), p2 * sqrt(1.0 - cth * cth) * sin(phi), p2 * cth, sqrt(Mp * Mp + p2 * p2) - dE);
    P2.RotateY(ki[0].Theta());
    P2.RotateZ(ki[0].Phi());
    const double Q = sqrt( (Md * Md - pow(ki[0].M() + P2.M(), 2)) * (Md * Md - pow(ki[0].M() - P2.M(), 2))) / (2.0 * Md);
    weight[0] *= pow(MODEL::FQ(Q), 2) * MODEL::FractionNphi;
    kf[0] = ki[0] + P2;
    weight[0] *= GOLD::ProtonDensity * ki[0].Gamma() / PARTICLE::phi.Gamma();
    //weight[0] *= GOLD::ProtonDensity * (7.3 / Phys::hbar / ki[0].Beta());
    return weight[0];
  }

  double BoundStatePhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: N', d
    TLorentzVector kf1[2];
    double weight1;
    PhiPhotoproductionGold(ki, kf1, &weight1);//produce phi
    if (weight1 == 0){
      weight[0] = 0;
      return 0;
    }
    kf[0] = kf1[1];//N'
    double weight2;
    BoundStateFormationGold(&kf1[0], &kf[1], &weight2);//form bound state
    weight[0] = weight1 * weight2 * GOLD::NA;
    return weight[0];
  }

  double ScatteredElectron(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', gamma
    const double Pmin = 0.5;
    const double Pmax = ki[0].P();
    const double thmin = 2.5 / 180.0 * M_PI;
    const double thmax = 4.5 / 180.0 * M_PI;
    const double Pe = random.Uniform(Pmin, Pmax);
    const double theta = acos(random.Uniform(cos(thmax), cos(thmin)));
    const double phi = random.Uniform(-M_PI, M_PI);
    const double Me = PARTICLE::e.M();
    const double Ee = sqrt(Pe * Pe - Me * Me);
    kf[0].SetXYZT(Pe * sin(theta) * cos(phi), Pe * sin(theta) * sin(phi), Pe * cos(theta), Ee);//e'
    kf[1] = ki[0] - kf[0];//photon
    const double alpha_em = 1.0 / 137.0;
    weight[0] = 16.0 * M_PI * alpha_em * (kf[0] * ki[0]);//vertex M square
    weight[0] *= 1.0 / pow(kf[1].M2(), 2);//propagator square
    weight[0] *= 2.0 * M_PI * (cos(thmin) - cos(thmax)) * (Pmax - Pmin) * Pe / (2.0 * pow(2.0 * M_PI, 3));
    return weight[0];
  }

  double BoundStateElectroproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', d
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1 = kf[1];//virtual photon
    double weight2;
    BoundStatePhotoproductionGold(&ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    TLorentzVector PA(0, 0, 0, GOLD::NA * Mp);
    weight[0] *= sqrt(pow(ki1 * PA, 2) - ki1.M2() * PA.M2()) / sqrt(pow(ki[0] * PA, 2) - ki[0].M2() * PA.M2());
    return weight[0];
  }

  double PhiElectroproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', phi, N'
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1 = kf[1];//virtual photon
    double weight2;
    PhiPhotoproductionGold(&ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    TLorentzVector PA(0, 0, 0, GOLD::NA * Mp);
    weight[0] *= sqrt(pow(ki1 * PA, 2) - ki1.M2() * PA.M2()) / sqrt(pow(ki[0] * PA, 2) - ki[0].M2() * PA.M2());
    return weight[0];
  }

  int NKKWeight(){
    FILE * fp = fopen("wave/NKK.dat", "w");
    double x, y;
    double mass[3] = {Mp, PARTICLE::K$p.M(), PARTICLE::K$m.M()};
    TLorentzVector PP(0, 0, 0, 0);
    double sum = 0;
    for (int i =0; i < 100; i++){
      cout << i << endl;
      x = Mp + PARTICLE::K.M() * 2.0 + 0.05 * i;
      PP.SetE(x);
      GenPhase.SetDecay(PP, 3, mass);
      sum = 0;
      for (int j = 0; j < 10000000; j++){
	sum += GenPhase.Generate();
      }
      y = sum / 10000000;
      fprintf(fp, "%.6E\t%.6E\n", x, 1.0 / y);
    }
    fclose(fp);
    return 0;
  }

  ROOT::Math::Interpolator fNKK_INTER(100, ROOT::Math::Interpolation::kCSPLINE);
  int SetfNKK(){
    ifstream infile("wave/NKK.dat");
    double x[100], y[100];
    for (int i = 0; i < 100; i++)
      infile >> x[i] >> y[i];
    fNKK_INTER.SetData(100, x, y);
    infile.close();
    return 0;
  }

  double fNKK(const double M){
    if (M <= 1.92563)
      return 0;
    return fNKK_INTER.Eval(M);
  }

  double KKPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma, N; kf: N', K+, K-
    TLorentzVector Pout = ki[0] + ki[1];
    const double MK = PARTICLE::K.M();
    if (Pout.M() <= Mp + MK + MK){
      weight[0] = 0;
      return weight[0];
    }
    double mass[3] = {Mp, MK, MK};
    GenPhase.SetDecay(Pout, 3, mass);
    weight[0] = GenPhase.Generate() * fNKK(Pout.M());//phase space weight
    kf[0] = *GenPhase.GetDecay(0);//N'
    kf[1] = *GenPhase.GetDecay(1);//K+
    kf[2] = *GenPhase.GetDecay(2);//K-
    TLorentzVector kk = kf[1] + kf[2];
    double W0 = kk.M() + Mp;
    TLorentzVector kp = ki[0];
    kp.Boost(-Pout.BoostVector());
    kk.Boost(-Pout.BoostVector());
    double cth = cos(kk.Angle(kp.Vect()));
    weight[0] *= dSigmaPhi(W0, Pout.M(), cth) * 4.0 * M_PI;//phi(KK) cross section
    weight[0] *= 0.489;//Branch ratio to K+K-
    return weight[0];
  }

  double KKPhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: N', K+, K-
    TLorentzVector ki1[2];
    ki1[0] = ki[0];//photon
    NucleonGold(&ki1[1]);//off-shell nucleon
    KKPhotoproduction(ki1, kf, weight);
    weight[0] *= GOLD::NA;
    return weight[0];//GeV^-2
  }

  double KKElectroproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', K+, K-
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1 = kf[1];//virtual photon
    double weight2;
    KKPhotoproductionGold(&ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    TLorentzVector PA(0, 0, 0, GOLD::NA * Mp);
    weight[0] *= sqrt(pow(ki1 * PA, 2) - ki1.M2() * PA.M2()) / sqrt(pow(ki[0] * PA, 2) - ki[0].M2() * PA.M2());
    return weight[0];
  }

  double dSigmaL1520(const double W0, const double W, const double cth){
    if (W <= W0) return 0;
    double par[5] = {11.299, 4.60959, 0.835621, 0.54681, 1.827941};
    double b1 = par[3] * (W - W0 + par[4]);
    double a2 = 0.25;
    double a0 = par[0] * (W - W0) * exp(-par[1] * pow(W - W0, par[2]));;//total cross section in unit mub
    double r0 = exp(b1 * cth) * b1 / (2.0 * sinh(b1));
    double r2 = cth * cth * exp(b1 * cth) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
    double ds = a0 * (r0 + a2 * r2) / (1.0 + a2) / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
    return ds;
  }

  double L1520Photoproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma, N; kf: K+, L1520
    TLorentzVector Pout = ki[0] + ki[1];
    const double MK = PARTICLE::K.M();
    const double ML = PARTICLE::Lambda1520.RandomM(MK + Mp + 1.0e-20, 1.7);
    if (Pout.M() <= MK + ML){
      weight[0] = 0;
      return weight[0];
    }
    double mass[2] = {MK, ML};
    GenPhase.SetDecay(Pout, 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);
    kf[1] = *GenPhase.GetDecay(1);
    TLorentzVector ka = ki[0];
    TLorentzVector kc = kf[0];
    ka.Boost(-Pout.BoostVector());
    kc.Boost(-Pout.BoostVector());
    double cth = cos(ka.Angle(kc.Vect()));
    weight[0] = dSigmaL1520(MK + ML, Pout.M(), cth) * 4.0 * M_PI;
    return weight[0];
  }

  double L1520PhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: K+, L1520
    TLorentzVector ki1[2];
    ki1[0] = ki[0];//photon
    NucleonGold(&ki1[1]);//off-shell nucleon
    L1520Photoproduction(ki1, kf, weight);
    weight[0] *= GOLD::NA;
    return weight[0];//GeV^-2
  }

  double L1520ElectroproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', K+, L1520
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1 = kf[1];//virtual photon
    double weight2;
    L1520PhotoproductionGold(&ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    TLorentzVector PA(0, 0, 0, GOLD::NA * Mp);
    weight[0] *= sqrt(pow(ki1 * PA, 2) - ki1.M2() * PA.M2()) / sqrt(pow(ki[0] * PA, 2) - ki[0].M2() * PA.M2());
    return weight[0];
  }

  double Event_eNKKN_BoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', [K+, K-, p]
    TLorentzVector kk[3];
    BoundStateElectroproductionGold(ki, kk, weight);
    kf[0] = kk[0];//e'
    kf[1] = kk[1];//N'
    const double MK = PARTICLE::K.M();
    double mass[3] = {Mp, MK, MK};
    GenPhase.SetDecay(kk[2], 3, mass);
    weight[0] *= GenPhase.Generate() * fNKK(kk[2].M());
    kf[2] = *GenPhase.GetDecay(1);//K+
    kf[3] = *GenPhase.GetDecay(2);//K-
    kf[4] = *GenPhase.GetDecay(0);//p
    weight[0] *= 0.465;//Branch ratio to pK+K-
    return weight[0];
  }

  double Event_eNKK_Phi(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', [K+, K-]
    TLorentzVector kk[3];
    PhiElectroproductionGold(ki, kk, weight);
    kf[0] = kk[0];//e'
    kf[1] = kk[2];//N'
    const double MK = PARTICLE::K.M();
    double mass[2] = {MK, MK};
    GenPhase.SetDecay(kk[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//K+
    kf[3] = *GenPhase.GetDecay(1);//K-
    weight[0] *= 0.489;//Branch ratio to K+K-
    weight[0] *= 79.0 / 197.0;//require proton
    return weight[0];
  }

  double Event_eNKK_KK(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', K+, K-
    KKElectroproductionGold(ki, kf, weight);
    weight[0] *= 79.0 / 197.0;
    return weight[0];
  }

  double Event_eNKK_L1520(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', K+, K-
    TLorentzVector kk[3];
    L1520ElectroproductionGold(ki, kk, weight);
    kf[0] = kk[0];//e'
    kf[2] = kk[1];//K+
    const double MK = PARTICLE::K.M();
    double mass[2] = {Mp, MK};
    GenPhase.SetDecay(kk[2], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//p
    kf[3] = *GenPhase.GetDecay(1);//K-
    weight[0] *= 0.45 / 2.0;//Branch ratio to pK-
    weight[0] *= 79.0 / 197.0;//require proton
    return weight[0];
  }


}


namespace DETECTOR{

  TFile * facc1;
  TFile * facc2;
  TH3F * acc_pip;
  TH3F * acc_pim;
  TH3F * acc_ele;

  int SetDETECTOR(){
    facc1 = new TFile("acceptance/clasev_acceptance.root", "r");
    facc2 = new TFile("acceptance/acceptance_ele_vertex_cP3375.root", "r");
    acc_pip = (TH3F *) facc1->Get("acceptance_PThetaPhi_pip");
    acc_pim = (TH3F *) facc1->Get("acceptance_PThetaPhi_pim");
    acc_ele = (TH3F *) facc1->Get("acceptance_PThetaPhi_ele");
    return 0;
  }

  double Acceptance(const TLorentzVector P, const char * part){
    double p = P.P();
    double theta = P.Theta() * 180.0 / M_PI;
    double phi = P.Phi() * 180.0 / M_PI;
    if (phi < 0) phi = phi + 360.0;
    TH3F * acc;
    if (strcmp(part, "p") == 0) acc = acc_pip;
    else if (strcmp(part, "e") == 0) acc = acc_ele;
    else if (strcmp(part, "K+") == 0) acc = acc_pip;
    else if (strcmp(part, "K-") == 0) acc = acc_pim;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(phi);
    int biny = acc->GetYaxis()->FindBin(theta);
    int binz = acc->GetZaxis()->FindBin(p);
    return acc->GetBinContent(binx, biny, binz);
  }




}



int Initialize(){
  MODEL::SetMODEL();
  GOLD::SetGOLD();
  GENERATE::SetfNKK();
  DETECTOR::SetDETECTOR();
  return 0;
}


#endif
