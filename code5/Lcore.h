#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TFile.h"
//#include "TTree.h"
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/SpecFuncMathCore.h"
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

  TF1 TF_fMass("fM", BreitWigner, PARTICLE::K.M() * 2.0 + Mp + 1.0e-20, Mass + 10.0 * Width, 0);

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

namespace PHIMODEL{
  ROOT::Math::Interpolator CLASm08(65, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm07(78, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm06(80, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm05(78, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm04(81, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm03(86, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm02(86, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASm01(85, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp00(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp01(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp02(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp03(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0375(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0425(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0475(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0525(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0575(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0625(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0675(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0725(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0775(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0825(87, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0875(86, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator CLASp0925(80, ROOT::Math::Interpolation::kCSPLINE);
  
  int SetInterpolation(){
    double x[100], y[100], tmp;
    ifstream fm08("data/dsc-0.8.dat");
    int i = 0;
    while (fm08 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm08.SetData(65, x, y);
    fm08.close();
    ifstream fm07("data/dsc-0.7.dat");
    i = 0;
    while (fm07 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm07.SetData(78, x, y);
    fm07.close();
    ifstream fm06("data/dsc-0.6.dat");
    i = 0;
    while (fm06 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm06.SetData(80, x, y);
    fm06.close();
    ifstream fm05("data/dsc-0.5.dat");
    i = 0;
    while (fm05 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm05.SetData(78, x, y);
    fm05.close();
    ifstream fm04("data/dsc-0.4.dat");
    i = 0;
    while (fm04 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm04.SetData(81, x, y);
    fm04.close();
    ifstream fm03("data/dsc-0.3.dat");
    i = 0;
    while (fm03 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm03.SetData(86, x, y);
    fm03.close();
    ifstream fm02("data/dsc-0.2.dat");
    i = 0;
    while (fm02 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm02.SetData(86, x, y);
    fm02.close();
    ifstream fm01("data/dsc-0.1.dat");
    i = 0;
    while (fm01 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASm01.SetData(85, x, y);
    fm01.close();
    ifstream fp00("data/dsc0.0.dat");
    i = 0;
    while (fp00 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp00.SetData(87, x, y);
    fp00.close();
    ifstream fp01("data/dsc0.1.dat");
    i = 0;
    while (fp01 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp01.SetData(87, x, y);
    fp01.close();
    ifstream fp02("data/dsc0.2.dat");
    i = 0;
    while (fp02 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp02.SetData(87, x, y);
    fp02.close();
    ifstream fp03("data/dsc0.3.dat");
    i = 0;
    while (fp03 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp03.SetData(87, x, y);
    fp03.close();
    ifstream fp0375("data/dsc0.375.dat");
    i = 0;
    while (fp0375 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0375.SetData(87, x, y);
    fp0375.close();
    ifstream fp0425("data/dsc0.425.dat");
    i = 0;
    while (fp0425 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0425.SetData(87, x, y);
    fp0425.close();
    ifstream fp0475("data/dsc0.475.dat");
    i = 0;
    while (fp0475 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0475.SetData(87, x, y);
    fp0475.close();
    ifstream fp0525("data/dsc0.525.dat");
    i = 0;
    while (fp0525 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0525.SetData(87, x, y);
    fp0525.close();
    ifstream fp0575("data/dsc0.575.dat");
    i = 0;
    while (fp0575 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0575.SetData(87, x, y);
    fp0575.close();
    ifstream fp0625("data/dsc0.625.dat");
    i = 0;
    while (fp0625 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0625.SetData(87, x, y);
    fp0625.close();
    ifstream fp0675("data/dsc0.675.dat");
    i = 0;
    while (fp0675 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0675.SetData(87, x, y);
    fp0675.close();
    ifstream fp0725("data/dsc0.725.dat");
    i = 0;
    while (fp0725 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0725.SetData(87, x, y);
    fp0725.close();
    ifstream fp0775("data/dsc0.775.dat");
    i = 0;
    while (fp0775 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0775.SetData(87, x, y);
    fp0775.close();
    ifstream fp0825("data/dsc0.825.dat");
    i = 0;
    while (fp0825 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0825.SetData(87, x, y);
    fp0825.close();
    ifstream fp0875("data/dsc0.875.dat");
    i = 0;
    while (fp0875 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0875.SetData(86, x, y);
    fp0875.close();
    ifstream fp0925("data/dsc0.925.dat");
    i = 0;
    while (fp0925 >> x[i] >> tmp >> y[i] >> tmp) i++;
    CLASp0925.SetData(80, x, y);
    fp0925.close();
    return 0;
  }

  double dSigmaPhi_clas(const double W0, const double W, const double cth){
    if (W <= W0) return 0;
    double result = 0.0;
    double dE = W - W0;
    if (dE > 0.85) dE = 0.85;
    if (cth < -0.75) result = CLASm08.Eval(dE);
    else if (cth >= -0.75 && cth < -0.65) result = CLASm07.Eval(dE);
    else if (cth >= -0.65 && cth < -0.55) result = CLASm06.Eval(dE);
    else if (cth >= -0.55 && cth < -0.45) result = CLASm05.Eval(dE);
    else if (cth >= -0.45 && cth < -0.35) result = CLASm04.Eval(dE);
    else if (cth >= -0.35 && cth < -0.25) result = CLASm03.Eval(dE);
    else if (cth >= -0.25 && cth < -0.15) result = CLASm02.Eval(dE);
    else if (cth >= -0.15 && cth < -0.05) result = CLASm01.Eval(dE);
    else if (cth >= -0.05 && cth < 0.05) result = CLASp00.Eval(dE);
    else if (cth >= 0.05 && cth < 0.15) result = CLASp01.Eval(dE);
    else if (cth >= 0.15 && cth < 0.25) result = CLASp02.Eval(dE);
    else if (cth >= 0.25 && cth < 0.35) result = CLASp03.Eval(dE);
    else if (cth >= 0.35 && cth < 0.40) result = CLASp0375.Eval(dE);
    else if (cth >= 0.40 && cth < 0.45) result = CLASp0425.Eval(dE);
    else if (cth >= 0.45 && cth < 0.50) result = CLASp0475.Eval(dE);
    else if (cth >= 0.50 && cth < 0.55) result = CLASp0525.Eval(dE);
    else if (cth >= 0.55 && cth < 0.60) result = CLASp0575.Eval(dE);
    else if (cth >= 0.60 && cth < 0.65) result = CLASp0625.Eval(dE);
    else if (cth >= 0.65 && cth < 0.70) result = CLASp0675.Eval(dE);
    else if (cth >= 0.70 && cth < 0.75) result = CLASp0725.Eval(dE);
    else if (cth >= 0.75 && cth < 0.80) result = CLASp0775.Eval(dE);
    else if (cth >= 0.80 && cth < 0.85) result = CLASp0825.Eval(dE);
    else if (cth >= 0.85 && cth < 0.90) result = CLASp0875.Eval(dE);
    else if (cth >= 0.90) result = CLASp0925.Eval(dE);
    else result = 0.0;
    return result / (2.0 * M_PI) / 389.379;//ds/dOmega in unit GeV^-2
  }
   
    
    
}

namespace GENERATE{

  TRandom3 random(0);
  TGenPhaseSpace GenPhase;
  double Weight = 0.0;

  double Bremsstrahlung(const double * k, const double * par){//non-normalized bremsstrahlung photon
    //E0: electron beam energy; k: photon energy
    double E0 = par[0];
    double y = k[0] / E0;
    if (y < 0.01) {
      std::cerr << "Out of range in Bremsstrahlung!" << std::endl;
      return -1.0;
    }
    double result = 1.0 / (y * E0) * (4.0 / 3.0 - 4.0 / 3.0 * y + y * y);
    return result;
  }
  
  TF1 TF_fBremsstrahlung("fBremsstrahlung", Bremsstrahlung, 1.0, 2.0, 1);//set photon energy distribution

  double BremsstrahlungPhoton(TLorentzVector * q, const double * k){//Generate a Bremsstrahlung photon
    //q: photon; k: Emin, Emax
    double E0 = TF_fBremsstrahlung.GetRandom(k[0], k[1]);
    q[0].SetXYZT(0.0, 0.0, E0, E0);
    return 1.0;
  }
  
  int NucleonGold(TLorentzVector * P){
    double p = GOLD::TF_fMomentum.GetRandom();
    double cth = random.Uniform(-1.0, 1.0);
    double phi = random.Uniform(-M_PI, M_PI);
    double dE = GOLD::TF_fEnergy.GetRandom();
    P->SetXYZT(p * sqrt(1.0 - cth * cth) * cos(phi), p * sqrt(1.0 - cth * cth) * sin(phi), p * cth, sqrt(p * p + Mp * Mp) - dE);
    return 0;
  }

  double FormFactor(const double tau, const double Q2){//twist-tau form factor
    if (Q2 < 0){
      cout << "Timelike region Q2!" << endl;
    }
    return ROOT::Math::beta(tau - 1.0, 0.5 + Q2 / (4.0 * pow(0.5337, 2))) / ROOT::Math::beta(tau - 1.0, 0.5);
  }

  double GEp(const double Q2){//proton electric form factor model
    double F1p = 1.1612 * FormFactor(3.0, Q2) - 0.1612 * FormFactor(4.0, Q2) + 0.0011 * FormFactor(5.0, Q2) - 0.0011 * FormFactor(6.0, Q2);
    double F2p = 1.4400 * FormFactor(4.0, Q2) + 0.3528 * FormFactor(6.0, Q2);
    return F1p - Q2 / (4.0 * Mp * Mp) * F2p;
  }

  double GMp(const double Q2){//proton magnetic form factor model
    double F1p = 1.1612 * FormFactor(3.0, Q2) - 0.1612 * FormFactor(4.0, Q2) + 0.0011 * FormFactor(5.0, Q2) - 0.0011 * FormFactor(6.0, Q2);
    double F2p = 1.4400 * FormFactor(4.0, Q2) + 0.3528 * FormFactor(6.0, Q2);
    return F1p + F2p;
  }

  double Event_ep_QE(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){//quasi-elastic ep event from gold target
    //ki: e; kf: e', p
    TLorentzVector P;
    NucleonGold(&P);//off-shell proton
    TLorentzVector l = ki[0];//incoming electron
    l.Boost(-P.BoostVector());//boost to proton rest frame
    TLorentzVector lp(0,0,1,1);//scattered electron, temporary value
    lp.SetTheta(acos(random.Uniform(-1.0, 1.0)));//polar angle
    lp.SetPhi(random.Uniform(-M_PI, M_PI));//azimuthal angle
    double th = lp.Angle(l.Vect());//angle w.r.t. the incoming electron in the proton rest frame
    double El = (pow(P.M(), 2) + 2.0 * l.E() * P.M() - pow(Mp, 2)) / (2.0 * (P.M() + l.E() * (1.0 - cos(th))));//energy of scattered electron
    lp.SetRho(El);//set momentum, electron mass neglected
    lp.SetE(El);
    TLorentzVector q = l - lp;//virtual photon
    double Q2 = -q * q;
    double alpha_em = 1.0 / 137.0;
    double sigma0 = 4.0 * pow(alpha_em, 2) * pow(cos(th/2.0), 2) * pow(lp.E(), 3) / l.E() / pow(Q2, 2);//point-like cross section, GeV^-2
    double tau = Q2 / (4.0 * Mp * Mp);
    double epsilon = 1.0 / (1.0 + 2.0 * (1.0 + tau) * pow(tan(th/2.0), 2));//photon polarization
    double sigma = sigma0 * (epsilon * pow(GEp(Q2), 2) + tau * pow(GMp(Q2), 2)) / (epsilon * (1.0 + tau));//cross section, GeV^-2
    double volume = 4.0 * M_PI;//phase space
    weight[0] = sigma * volume;//weighting factor
    weight[0] *= 79.0 / 197.0;//proton fraction in gold
    lp.Boost(P.BoostVector());//boost to lab frame
    kf[0] = lp;//scattered electron
    kf[1] = ki[0] + P - kf[0];//outgoing proton
    return weight[0];
  }

  double (*dSigmaPhi)(const double, const double, const double);

  double dSigmaPhi_old(const double W0, const double W, const double cth){
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
    const double Pmin = 0.05;
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

  double PhiElectroproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', phi, N'
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1[2];
    ki1[0] = kf[1];//virtual photon
    ki1[1].SetXYZT(0, 0, 0, Mp);//rest nucleon
    double weight2;
    PhiPhotoproduction(ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    weight[0] *= sqrt(pow(ki1[0] * ki1[1], 2) - ki1[0].M2() * ki1[1].M2()) / sqrt(pow(ki[0] * ki1[1], 2) - ki[0].M2() * ki1[1].M2());
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

  ROOT::Math::Interpolator dSigmaPiPi_INTER(117, ROOT::Math::Interpolation::kCSPLINE);
  int SetdSigmaPiPi(){
    ifstream infile("wave/protonpippimtotal.dat");
    double x[117], y[117], w[117];
    for (int i = 0; i < 117; i++){
      infile >> x[i] >> y[i];
      w[i] = sqrt(Mp * Mp + 2.0 * Mp * x[i]);
      y[i] *= 0.1 / pow(Phys::hbar, 2);//convert to GeV^-2
    }
    dSigmaPiPi_INTER.SetData(117, w, y);
    infile.close();
    return 0;
  }

  double dSigmaPiPi(const double W){
    if (W <= Mp + PARTICLE::pi.M() * 2.0) return 0;
    return dSigmaPiPi_INTER.Eval(W);
  }

  int NPiPiWeight(){
    FILE * fp = fopen("wave/NPiPi.dat", "w");
    double x, y;
    double mass[3] = {Mp, PARTICLE::pi$p.M(), PARTICLE::pi$m.M()};
    TLorentzVector PP(0, 0, 0, 0);
    double sum = 0;
    for (int i =0; i < 100; i++){
      cout << i << endl;
      x = Mp + PARTICLE::pi.M() * 2.0 + 0.05 * i;
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

  ROOT::Math::Interpolator fNPiPi_INTER(100, ROOT::Math::Interpolation::kCSPLINE);
  int SetfNPiPi(){
    ifstream infile("wave/NPiPi.dat");
    double x[100], y[100];
    for (int i = 0; i < 100; i++)
      infile >> x[i] >> y[i];
    fNPiPi_INTER.SetData(100, x, y);
    infile.close();
    return 0;
  }

  double fNPiPi(const double M){
    if (M <= Mp + PARTICLE::pi.M() * 2.0)
      return 0;
    return fNPiPi_INTER.Eval(M);
  }

  double PiPiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma, N; kf: N', pi+, pi-
    TLorentzVector Pout = ki[0] + ki[1];
    const double Mpi = PARTICLE::pi.M();
    if (Pout.M() <= Mp + Mpi + Mpi){
      weight[0] = 0;
      return weight[0];
    }
    double mass[3] = {Mp, Mpi, Mpi};
    GenPhase.SetDecay(Pout, 3, mass);
    weight[0] = GenPhase.Generate() * fNPiPi(Pout.M());//phase space weight
    kf[0] = *GenPhase.GetDecay(0);//N'
    kf[1] = *GenPhase.GetDecay(1);//pi+
    kf[2] = *GenPhase.GetDecay(2);//pi-
    weight[0] *= dSigmaPiPi(Pout.M());//phi(KK) cross section
    return weight[0];
  }

  double PiPiPhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: N', pi+, pi-
    TLorentzVector ki1[2];
    ki1[0] = ki[0];//photon
    NucleonGold(&ki1[1]);//off-shell nucleon
    PiPiPhotoproduction(ki1, kf, weight);
    weight[0] *= GOLD::NA;
    return weight[0];//GeV^-2
  }

  double PiPiElectroproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', pi+, pi-
    double weight1;
    ScatteredElectron(ki, kf, &weight1);
    TLorentzVector ki1 = kf[1];//virtual photon
    double weight2;
    PiPiPhotoproductionGold(&ki1, &kf[1], &weight2);
    weight[0] = weight1 * weight2;
    TLorentzVector PA(0, 0, 0, GOLD::NA * Mp);
    weight[0] *= sqrt(pow(ki1 * PA, 2) - ki1.M2() * PA.M2()) / sqrt(pow(ki[0] * PA, 2) - ki[0].M2() * PA.M2());
    return weight[0];
  }

  double Event_NKKN_BoundState(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: q; kf: N', [K+, K-, p]
    TLorentzVector kk[2];
    BoundStatePhotoproductionGold(ki, kk, weight);
    kf[0] = kk[0];//N'
    const double MK = PARTICLE::K.M();
    double mass[3] = {Mp, MK, MK};
    GenPhase.SetDecay(kk[1], 3, mass);
    weight[0] *= GenPhase.Generate() * fNKK(kk[1].M());
    kf[1] = *GenPhase.GetDecay(1);//K+
    kf[2] = *GenPhase.GetDecay(2);//K-
    kf[3] = *GenPhase.GetDecay(0);//p
    weight[0] *= 0.465;//Branch ratio to pK+K-
    return weight[0];
  }

  double Event_NKK_Phi(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: q; kf: N', [K+, K-]
    TLorentzVector kk[2];
    PhiPhotoproductionGold(ki, kk, weight);
    kf[0] = kk[1];//N'
    const double MK = PARTICLE::K.M();
    double mass[2] = {MK, MK};
    GenPhase.SetDecay(kk[0], 2, mass);
    GenPhase.Generate();
    kf[1] = *GenPhase.GetDecay(0);//K+
    kf[2] = *GenPhase.GetDecay(1);//K-
    weight[0] *= 0.489;//Branch ratio to K+K-
    weight[0] *= 79.0 / 197.0;//require proton
    return weight[0];
  }

  double Event_NKK_KK(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: q; kf: N', K+, K-
    KKPhotoproductionGold(ki, kf, weight);
    weight[0] *= 79.0 / 197.0;//require proton
    return weight[0];
  }

  double Event_NKK_L1520(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: q; kf: N', K+, K-
    TLorentzVector kk[2];
    L1520PhotoproductionGold(ki, kk, weight);
    kf[1] = kk[0];//K+
    const double MK = PARTICLE::K.M();
    double mass[2] = {Mp, MK};
    GenPhase.SetDecay(kk[1], 2, mass);
    GenPhase.Generate();
    kf[0] = *GenPhase.GetDecay(0);//p
    kf[2] = *GenPhase.GetDecay(1);//K-
    weight[0] *= 0.45 / 2.0;//Branch ratio to pK-
    weight[0] *= 79.0 / 197.0;//require proton
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

  double Event_eNKK_Phi_Nucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', [K+, K-]
    TLorentzVector kk[3];
    PhiElectroproduction(ki, kk, weight);
    kf[0] = kk[0];//e'
    kf[1] = kk[2];//N'
    const double MK = PARTICLE::K.M();
    double mass[2] = {MK, MK};
    GenPhase.SetDecay(kk[1], 2, mass);
    GenPhase.Generate();
    kf[2] = *GenPhase.GetDecay(0);//K+
    kf[3] = *GenPhase.GetDecay(1);//K-
    weight[0] *= 0.489;//Branch ratio to K+K-
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

  double Event_eNPiPi_PiPi(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: e; kf: e', N', pi+, pi-
    PiPiElectroproductionGold(ki, kf, weight);
    weight[0] *= 79.0 / 197.0;
    return weight[0];
  }

  int SetGENERATE(){
    TF_fBremsstrahlung.SetParameter(0, 11.0);//Set electron energy
    TF_fBremsstrahlung.SetNpx(500);
    dSigmaPhi = &PHIMODEL::dSigmaPhi_clas;
    //dSigmaPhi = &dSigmaPhi_old;
    SetfNKK();
    SetdSigmaPiPi();
    SetfNPiPi();
    return 0;
  }


}


namespace DETECTOR{

  TRandom3 random(0);

  TFile * facc1, * facc2, * facc3;
  TFile * fres1, * fres2;
  TH3F * acc_ele_clas;
  TH3F * acc_pip_clas;
  TH3F * acc_pim_clas;
  TH3F * acc_Kp_clas;
  TH3F * acc_Km_clas;
  TH3F * acc_proton_clas;
  TH2D * acc_proton_alert;
  TH2D * acc_Kp_alert;
  TH2D * acc_pip_alert;
  TH2D * acc_Km_alert;
  TH2D * acc_pim_alert;
  TH2D * res_Kp_alert_p;
  TH2D * res_Kp_alert_theta;
  TH2D * res_Kp_alert_phi;
  TH2D * res_Km_alert_p;
  TH2D * res_Km_alert_theta;
  TH2D * res_Km_alert_phi;
  TH2D * res_pip_alert_p;
  TH2D * res_pip_alert_theta;
  TH2D * res_pip_alert_phi;
  TH2D * res_pim_alert_p;
  TH2D * res_pim_alert_theta;
  TH2D * res_pim_alert_phi;
  TH2D * res_proton_alert_p;
  TH2D * res_proton_alert_theta;
  TH2D * res_proton_alert_phi;

  int SetDETECTOR(){
    facc1 = new TFile("acceptance/clasev_acceptance_binP20MeVTheta1degPhi1deg.root", "r");
    facc2 = new TFile("acceptance/acceptance_ele_vertex_cP3375.root", "r");
    facc3 = new TFile("acceptance/acc_alert_20190427.root", "r");
    fres1 = new TFile("acceptance/res_kp_20190429.root", "r");
    fres2 = new TFile("acceptance/res_proton_20190429.root", "r");
    acc_pip_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_pip");
    acc_pim_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_pim");
    acc_ele_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_ele");
    acc_proton_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_pip");
    acc_Kp_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_pip");
    acc_Km_clas = (TH3F *) facc1->Get("acceptance_PThetaPhi_pim");
    acc_proton_alert = (TH2D *) facc3->Get("h0");
    acc_Kp_alert = (TH2D *) facc3->Get("h1");
    acc_pip_alert = (TH2D *) facc3->Get("h2");
    acc_Km_alert = (TH2D *) facc3->Get("h3");
    acc_pim_alert = (TH2D *) facc3->Get("h4");
    res_Kp_alert_p = (TH2D *) fres1->Get("h1");
    res_Kp_alert_theta = (TH2D *) fres1->Get("h2");
    res_Kp_alert_phi = (TH2D *) fres1->Get("h3");
    res_Km_alert_p = (TH2D *) fres1->Get("h1");
    res_Km_alert_theta = (TH2D *) fres1->Get("h2");
    res_Km_alert_phi = (TH2D *) fres1->Get("h3");
    res_pip_alert_p = (TH2D *) fres1->Get("h1");
    res_pip_alert_theta = (TH2D *) fres1->Get("h2");
    res_pip_alert_phi = (TH2D *) fres1->Get("h3");
    res_pim_alert_p = (TH2D *) fres1->Get("h1");
    res_pim_alert_theta = (TH2D *) fres1->Get("h2");
    res_pim_alert_phi = (TH2D *) fres1->Get("h3");
    res_proton_alert_p = (TH2D *) fres2->Get("h1");
    res_proton_alert_theta = (TH2D *) fres2->Get("h2");
    res_proton_alert_phi = (TH2D *) fres2->Get("h3");
    return 0;
  }

  double AcceptanceCLAS12(const TLorentzVector P, const char * part){
    double p = P.P();
    double theta = P.Theta() * 180.0 / M_PI;
    if (theta > 35.0) return 0;//only forward detector
    double phi = P.Phi() * 180.0 / M_PI;
    if (phi < 0) phi = phi + 360.0;
    TH3F * acc;
    if (strcmp(part, "e") == 0) acc = acc_ele_clas;
    else if (strcmp(part, "p") == 0) acc = acc_proton_clas;
    else if (strcmp(part, "K+") == 0) acc = acc_Kp_clas;
    else if (strcmp(part, "K-") == 0) acc = acc_Km_clas;
    else if (strcmp(part, "pi+") == 0) acc = acc_pip_clas;
    else if (strcmp(part, "pi-") == 0) acc = acc_pim_clas;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(phi);
    int biny = acc->GetYaxis()->FindBin(theta);
    int binz = acc->GetZaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny, binz);
    if (strcmp(part, "K+") == 0 || strcmp(part, "K-") == 0) result *= exp(-6.5 / Phys::c / PARTICLE::K.Tau() / P.Beta() / P.Gamma());//kaon decay
    return result;
  }

  double AcceptanceALERT(const TLorentzVector P, const char * part){
    double p = P.P() * 1000.0;//MeV
    if (p > 350.0) return 0;//sharp cut on momentum
    double theta = P.Theta() * 180.0 / M_PI;
    TH2D * acc;
    if (strcmp(part, "p") == 0) acc = acc_pip_alert;
    else if (strcmp(part, "K+") == 0) acc = acc_Kp_alert;
    else if (strcmp(part, "K-") == 0) acc = acc_Km_alert;
    else if (strcmp(part, "pi+") == 0) acc = acc_pip_alert;
    else if (strcmp(part, "pi-") == 0) acc = acc_pim_alert;
    else return 0;
    int binx = acc->GetXaxis()->FindBin(theta);
    int biny = acc->GetYaxis()->FindBin(p);
    double result = acc->GetBinContent(binx, biny);
    return result;
  }

  double Acceptance(const TLorentzVector P, const char * part){
    double acc_clas = AcceptanceCLAS12(P, part);
    double acc_alert = AcceptanceALERT(P, part);
    return 1.0 - (1.0 - acc_clas) * (1.0 - acc_alert);
  }

  double Smear(TLorentzVector * P, const char * part){
    double m = P->M();
    double p = P->P();
    double theta = P->Theta();
    double phi = P->Phi();
    double res[3];
    double acc = AcceptanceALERT(P[0], part);
    TH2D * resp, * restheta, * resphi;
    if (acc > 0){
      if (strcmp(part, "p") == 0){
	resp = res_proton_alert_p;
	restheta = res_proton_alert_theta;
	resphi = res_proton_alert_phi;
      }
      else if (strcmp(part, "K+") == 0){
	resp = res_Kp_alert_p;
	restheta = res_Kp_alert_theta;
	resphi = res_Kp_alert_phi;
      }
      else if (strcmp(part, "K-") == 0){
	resp = res_Km_alert_p;
	restheta = res_Km_alert_theta;
	resphi = res_Km_alert_phi;
      }
      else if (strcmp(part, "pi+") == 0){
        resp = res_pip_alert_p;
        restheta = res_pip_alert_theta;
        resphi = res_pip_alert_phi;
      }
      else if (strcmp(part, "pi-") == 0){
        resp = res_pim_alert_p;
        restheta = res_pim_alert_theta;
        resphi = res_pim_alert_phi;
      }
      else
	return 0;
      res[0] = resp->GetBinContent(resp->GetXaxis()->FindBin(theta / M_PI * 180.0), resp->GetYaxis()->FindBin(p * 1000.0)) / 100.0;
      res[1] = restheta->GetBinContent(restheta->GetXaxis()->FindBin(theta / M_PI * 180.0), restheta->GetYaxis()->FindBin(p * 1000.0));
      res[2] = resphi->GetBinContent(resphi->GetXaxis()->FindBin(theta / M_PI * 180.0), resphi->GetYaxis()->FindBin(p * 1000.0));
      p = p * abs(random.Gaus(1, res[0]));
      theta = random.Gaus(theta, res[1]);
      phi = random.Gaus(phi, res[2]);
      P->SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), m);
      return acc;
    }
    acc = AcceptanceCLAS12(P[0], part);
    if (acc > 0){
      res[0] = 0.01;
      res[1] = 0.001;
      res[2] = 0.004;
      p = p * abs(random.Gaus(1, res[0]));
      theta = random.Gaus(theta, res[1]);
      phi = random.Gaus(phi, res[2]);
      P->SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), m);
      return acc;
    }
    return 0;
  }


}



int Initialize(){
  MODEL::SetMODEL();
  GOLD::SetGOLD();
  PHIMODEL::SetInterpolation();
  GENERATE::SetGENERATE();
  DETECTOR::SetDETECTOR();
  return 0;
}


#endif
