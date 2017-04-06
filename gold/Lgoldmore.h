#ifndef _LGOLDMORE_H_
#define _LGOLDMORE_H_

#include "Lgold.h"

ROOT::Math::Interpolator Sigma_ProtonPionpPionm(117, ROOT::Math::Interpolation::kCSPLINE);
int SetSigmaProtonPionpPionm(){
  std::ifstream fppipi("protonpippimtotal.dat");
  double Ephoton[117], sigma[117];
  for (int i = 0; i < 117; i++){
    fppipi >> Ephoton[i] >> sigma[i];
    sigma[i] = sigma[i] / 0.389379;//convert to GeV^-2
  }
  Sigma_ProtonPionpPionm.SetData(117, Ephoton, sigma);
  fppipi.close();
  return 0;
}

double GeneratePhotoproductionpipiNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: photon, p; kf: p, pi+, pi-
  TLorentzVector Pout = ki[0] + ki[1];//Total 4-momentum
  const double Mp = 0.938272;//proton mass
  const double Mpi = 0.13957;//charged pion mass
  if (Pout.M() < Mp + Mpi + Mpi){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mpi, Mpi};
  double wd = Decay3(&Pout, kf, mass);//Get phase space weight factor
  double Eeff = (Pout.M2() - Mp * Mp) / (2.0 * Mp);//equivalent photon energy
  weight[0] = wd * Sigma_ProtonPionpPionm.Eval(Eeff);//total weight
  return weight[0];
}

double GenerateElectroproductionpipiNucleon(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: e, p; kf: e', p, pi+, pi-
  const double Mp = 0.938272;
  const double Mpi = 0.13957;
  TLorentzVector kf0[2];
  double w0 = GenerateScatteredElectron(ki, kf0);//Get amp * Lips weight
  if (w0 == 0){//unphysical
    weight[0] = 0;
    return 0;
  }
  kf[0] = kf0[0];//Set scattered electron
  double prop = 1.0 / pow(kf0[1].M2(), 2);//virtual photon propagator square
  TLorentzVector ki1[2] = {kf0[1], ki[1]};//virtual photon, nucleon
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of Ppipi
  if (Pout.M() < Mp + Mpi + Mpi){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mpi, Mpi};
  double wd = Decay3(&Pout, &kf[1], mass);//Get phase space weight factor
  double virtualFlux = 4.0 * sqrt(pow(ki1[0] * ki1[1], 2) - ki1[0].M2() * ki1[1].M2());
  double Flux = 4.0 * (ki[0] * ki[1]);
  double Eeff = (Pout.M2() - Mp * Mp) / (2.0 * Mp);//equivalent photon energy
  weight[0] = w0 * prop * wd * Sigma_ProtonPionpPionm.Eval(Eeff) * virtualFlux / Flux;//total weight
  return weight[0];
}

double GeneratePhotoproductionpipiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: photon; kf: p, pi+, pi-
  const double Mp = 0.938272;
  const double Mpi = 0.13957;
  const double MA = 197.0 * Mp;
  TLorentzVector ki1[2];
  ki1[0] = ki[0];//photon
  GenerateNucleonInGold(&ki1[1]);//reaction nucleon
  double Normal = MA / ki1[1].E();//invariant wave function normalization
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of Ppipi
  if (Pout.M() < Mp + Mpi + Mpi){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mpi, Mpi};
  double wd = Decay3(&Pout, kf, mass);//Get phase space weight factor
  double virtualFlux = 4.0 * (ki1[0] * ki1[1]);
  double Eeff = (Pout.M2() - Mp * Mp) / (2.0 * Mp);//equivalent photon energy
  double Flux = 4.0 * (ki[0].E() * MA);
  weight[0] = Normal * wd * Sigma_ProtonPionpPionm.Eval(Eeff) * virtualFlux / Flux;
  return weight[0];
}

double GenerateElectroproductionpipiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){
  //ki: e; kf: e', p, pi+, pi-
  const double Mp = 0.938272;
  const double Mpi = 0.13957;
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
  TLorentzVector Pout = ki1[0] + ki1[1];//total 4-momentum of Ppipi
  if (Pout.M() < Mp + Mpi + Mpi){//below threshold
    weight[0] = 0;
    return 0;
  }
  double mass[3] = {Mp, Mpi, Mpi};
  double wd = Decay3(&Pout, &kf[1], mass);//Get phase space weight factor
  double virtualFlux = 4.0 * sqrt(pow(ki1[0] * ki1[1], 2) - ki1[0].M2() * ki1[1].M2());
  double Eeff = (Pout.M2() - Mp * Mp) / (2.0 * Mp);//equivalent photon energy
  double Flux = 4.0 * (ki[0].E() * MA);
  weight[0] = Normal * w0 * prop * wd * Sigma_ProtonPionpPionm.Eval(Eeff) * virtualFlux / Flux;//total weight
  return weight[0];
}

double GenerateEvent_Npipi_withPhotoproductionpipiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with Npipi in final state with pipi production
  //ki: q; kf: N, pi+, pi-
  const double scale = 197.0;//nucleon number scale
  double w1 = GeneratePhotoproductionpipiGold(ki, kf);//Get photoproduction pipi weight
  weight[0] = scale * w1;
  return weight[0];
}

double GenerateEvent_eNpipi_withElectroproductionpipiGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &_weight_){//Generate an event with eNpipi in final state with pipi electroproduction
  //ki: e; kf: e', N, pi+, pi-
  const double scale = 197.0;//nucleon number scale
  double w1 = GenerateElectroproductionpipiGold(ki, kf);//Get electroproduction pipi weight
  weight[0] = scale * w1;
  return weight[0];
}  

int MisPIDtoKaon(TLorentzVector * P){
  double Px = P->X();
  double Py = P->Y();
  double Pz = P->Z();
  const double Mkaon = 0.493677;
  P->SetXYZM(Px, Py, Pz, Mkaon);
  return 0;
}







#endif
