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
    const double Mphi = PARTICLE::phi.RandomM();//random phi meson mass 
    if (s <= pow(Mphi + Mp, 2)){
      weight[0] = 0;
      return 0;
    }
    const double cth = random.Uniform(-1.0, 1.0);
    const double phi = random.Uniform(-M_PI, M_PI);
    const double Q = sqrt( (s - pow(Mphi + Mp, 2)) * (s - pow(Mphi - Mp, 2))) / (2.0 * Pout.M());
    weight[0] = dSigmaPhi(Mphi + Mp, Pout.M(), cth) * 4.0 * M_PI;
    kf[0].SetXYZM(Q * sqrt(1.0 - cth * cth) * cos(phi), Q * sqrt(1.0 - cth * cth) * sin(phi), Q * cth, Mphi);//phi in c.m. along gamma
    TLorentzVector kk = ki[0];
    kk.Boost(-Pout.BoostVector());//gamma in c.m.
    kf[0].RotateY(kk.Theta());
    kf[0].RotateZ(kk.Phi());//phi in c.m.
    kf[0].Boost(Pout.BoostVector());//phi
    kf[1] = Pout - kf[0];//N'
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
    return weight[0];
  }

  double BoundStatePhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: N', d
    TLorentzVector kf1[2];
    double weight1;
    PhiPhotoproductionGold(ki, kf1, &weight1);//produce Jpsi
    if (weight1 == 0){
      weight[0] = 0;
      return 0;
    }
    kf[0] = kf1[1];//N'
    double weight2;
    BoundStateFormationGold(&kf1[0], &kf[1], &weight2);//form bound state
    weight[0] = weight1 * weight2;
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


    
}




#endif
