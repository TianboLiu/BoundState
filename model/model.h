#ifndef _MODEL_H_
#define _MODEL_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>

#include "TF1.h"
#include "TF3.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

/************** Constants ***************/
const double GeVfm = 1.0 / 0.1973269718;//GeV times fm
const double RadToDeg = 180.0 / M_PI;
const double DegToRad = M_PI / 180.0;
const double Mp = 0.9382720;//proton mass
const double Mn = 0.9395654;//neutron mass
const double MN = (Mp + Mn) / 2.0;//average nucleon mass
const double MDelta = 1.232;//Delta baryon mass
const double Mphi = 1.019455;//phi meson mass
const double Metac = 2.9810;//eta_c meson mass
const double Mpsi = 3.096916;//J/psi meson mass
/*********** End of Constants ***********/

/************** Parameters **************/
double beta;//wavefunction parameter
double NA;//nucleon number in nucleus
double Lambda;//cut off parameters
double Mnuclear;//nucleus mass
double Mmeson;//meson mass
double Mreduced;//reduced mass
double w;//harmonic oscillator parameter
double Md;
int L1, L2;//nucleon orbital number
/********** End of Parameters ***********/

/************** Functions ***************/
const int _MAX0_ = 200;
double _ix0[_MAX0_], _iy0[_MAX0_];//Veff grid
gsl_interp_accel * accV;
gsl_spline * splineV;

const int _MAX1_ = 300;
double _ix1[_MAX1_], _iy1[_MAX1_], _iz1[_MAX1_];//wf1 grid
gsl_interp_accel * accR;
gsl_spline * splineR;
gsl_interp_accel * accU;
gsl_spline * splineU;

const int _MAX2_ = 200;
double _ix2[_MAX2_], _iy2[_MAX2_];//du_r_2 grid
gsl_interp_accel * accDU2;
gsl_spline * splineDU2;


int readgrid(){
  std::ifstream f0("dwf1.dat");
  if (!f0.is_open()){
    std::cerr << "File Veff.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX0_; i++){
    f0 >> _ix0[i] >> _iy0[i];
  }
  f0.close();
  accV = gsl_interp_accel_alloc();
  splineV = gsl_spline_alloc(gsl_interp_cspline, _MAX0_);
  gsl_spline_init(splineV, _ix0, _iy0, _MAX0_);

  std::ifstream f1("wf1.dat");
  if (!f1.is_open()){
    std::cerr << "File wf1.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX1_; i++){
    f1 >> _ix1[i] >> _iy1[i] >> _iz1[i];
  }
  f1.close();
  accR = gsl_interp_accel_alloc();
  splineR = gsl_spline_alloc(gsl_interp_cspline, _MAX1_);
  gsl_spline_init(splineR, _ix1, _iy1, _MAX1_);
  accU = gsl_interp_accel_alloc();
  splineU = gsl_spline_alloc(gsl_interp_cspline, _MAX1_);
  gsl_spline_init(splineU, _ix1, _iz1, _MAX1_);

  std::ifstream f2("dwf1.dat");
  if (!f2.is_open()){
    std::cerr << "File dwf1.dat does not exist!" << std::endl;
    return 1;
  }
  for (int i = 0; i < _MAX2_; i++){
    f2 >> _ix2[i] >> _iy2[i];
  }
  f2.close();
  accDU2 = gsl_interp_accel_alloc();
  splineDU2 = gsl_spline_alloc(gsl_interp_cspline, _MAX2_);
  gsl_spline_init(splineDU2, _ix2, _iy2, _MAX2_);

  return 0;
}

/* double Veff(double r){ */
/*   double rf = r / GeVfm; */
/*   //double rf = r; */
/*   if (rf > 3.0) return 0; */
/*   double y = gsl_spline_eval(splineV, rf, accV); */
/*   double result  = (y - _iy0[_MAX0_-1]) / 1000.0; */
/*   return result; */
/* } */

double Rr(double r, void * params = 0){
  double rf = r / GeVfm;
  if (rf > 15.0) return 0;
  double y = gsl_spline_eval(splineR, rf, accR);
  double result = y;
  return result;
}

double ur(double r, void * params = 0){
  double rf = r / GeVfm;
  if (rf > 15.0) return 0;
  double y = gsl_spline_eval(splineU, rf, accU);
  double result = y * GeVfm;
  return result;
}

double ur_r(double r, void * params = 0){
  double res, err;
  gsl_function F;
  F.function = &ur;
  F.params = 0;
  gsl_deriv_forward(&F, r, 1e-3, &res, &err);
  return res;
}

double ur_r_2(double r, void * params = 0){
  double rf = r / GeVfm;
  if (rf > 10.0) return 0;
  double y = gsl_spline_eval(splineDU2, rf, accDU2);
  double result = y / GeVfm;
  return result;
}

double integrand(double * rp, double * pars){
  double r = rp[0];
  //double result = -ur(r) * ur_r_2(r) / (2.0 * Mreduced) + ur(r) * ur(r) * Veff(r);
  //double result = ur(r) * ur(r);
  //double result = ur(r) * ur(r) * Veff(r);
  double result = -ur(r) * ur_r_2(r) / (2.0 * Mreduced);
  return result;
}

double one(){
  TF1 f("Integrand", integrand, 0.0, 100.0, 0);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GaussLegendreIntegrator ig;
  ig.SetFunction(wf1);
  ig.SetRelTolerance(0.001);
  ig.SetNumberPoints(60);
  double result = ig.Integral(0.0, 30.0);
  return result;
}

/* double shooting(double E, double rmax, bool wr = false){ */
/*   const int points = 6000; */
/*   double h = rmax / points; */
/*   double hh = h * h; */
/*   double x[points + 1], y[points + 1]; */
/*   double f[points + 1]; */
/*   for(int i = 0; i <= points; i++){ */
/*     x[i] = h * i; */
/*     f[i] = 2.0 * Mreduced * (E - Veff(x[i])); */
/*   } */
/*   y[0] = 0; */
/*   y[1] = h; */
/*   for (int i = 2; i <= points; i++){ */
/*     y[i] = ( 2.0 * y[i-1] * (1.0 - 5.0 * hh / 12.0 * f[i-1]) - y[i-2] * (1.0 + hh/12.0 * f[i-2]) ) / (1.0 + hh / 12.0 * f[i]); */
/*   } */
/*   printf("%8.6f   %18.6E   %18.6E   %18.6E\n", E, y[points/2], y[points], pow((y[points] - y[points-1])/h, 2)); */
/*   if (wr){ */
/*     FILE * fw; */
/*     fw = fopen("ur.dat","w"); */
/*     for (int i = 0; i <= points; i++) */
/*       fprintf(fw, "%E  %E\n", x[i], y[i]); */
/*     fclose(fw); */
/*   } */
/*   return y[points]; */
/* } */

int allfree(){
  gsl_spline_free(splineV);
  gsl_interp_accel_free(accV);
  gsl_spline_free(splineR);
  gsl_interp_accel_free(accR);
  gsl_spline_free(splineU);
  gsl_interp_accel_free(accU);
  gsl_spline_free(splineDU2);
  gsl_interp_accel_free(accDU2);
  return 0;
}


/*
double Psi(const double r){
  //double result = exp(-r*r / (2.0 * beta * beta)) / pow(M_PI * beta * beta, 3.0/4.0);
  double result = sqrt(pow(beta, 3) / M_PI) * exp(-beta * r);
  return result;
}

double EnergyT(const double b){
  beta = b;
  double m = Mreduced;
  //double result = 3.0 / (4.0 * m * beta * beta);
  double result = beta * beta / (2.0 * m); 
  return result;
}

double EVr(const double * rp, const double * foo){
  double r = rp[0];
  double wf2 = pow(Psi(r), 2);
  return 4.0 * M_PI * r * r * wf2 * Veff(r);
}

double EnergyV(const double b){
  beta = b;
  TF1 f("Integrand", EVr, 0.0, 18.0, 0);
  ROOT::Math::WrappedTF1 wf1(f);
  ROOT::Math::GaussLegendreIntegrator ig;
  ig.SetFunction(wf1);
  ig.SetRelTolerance(0.001);
  ig.SetNumberPoints(40);
  double result = ig.Integral(0.0, 18.0);
  return result;
}

double Ebind(const double * b){
  return EnergyT(b[0]) + EnergyV(b[0]);
}

double findbeta(){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
  min->SetMaxFunctionCalls(1000);
  min->SetMaxIterations(100);
  min->SetTolerance(1.0e-5);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Ebind, 1);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "beta", 2.0, 0.001, 0.5, 10.0);
  min->Minimize();
  const double * xs = min->X();
  beta = xs[0];
  return xs[0];
}



/*
double Ebind(const double x){
  //binding energy of Nuclear-meson in GeV for variation
  double m = Mreduced;
  return x * x * (1 / (2.0 * m) - 4.0 * alpha * x / (pow(2.0 * x + mu, 2)));
}

double vari_g(const double * x){
  double cc = x[0];
  return Ebind(cc);
}
double findgamma(){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
  min->SetMaxFunctionCalls(1000);
  min->SetMaxIterations(100);
  min->SetTolerance(1.0e-5);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&vari_g, 1);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "gamma", 1.0, 0.001, 0.0, 5.0);
  min->Minimize();
  const double * xs = min->X();
  gam = xs[0];
  return xs[0];
}

double ampF(const double Q){
  //phi-N -> d amplitude
  double num = - alpha * pow(gam, 3.0 / 2.0) * exp(-Q*Q/Lambda/Lambda); 
  double den = sqrt(2.0) * pow(M_PI, 3.0/2.0) * (Q*Q + pow(gam + mu, 2.0));
  return num / den;
}

double wf_N_iso(const double p, const int L){
  //nucleon wavefunction spherical symmetric
  double m = MN;
  double wf;
  if (L == 0)
    wf = exp(-p*p / (2.0 * m * w)) / pow(M_PI * m * w, 3.0/4.0);
  else if (L == 1)
    wf = sqrt(2.0 / (3.0 * m * w)) * p * exp(-p*p / (2.0 * m * w)) / pow(M_PI * m * w, 3.0/4.0);
  else wf = 0;
  return wf;
}

double ampT(const double qc){
  //gamma N -> phi N amplitude
  double t0 = 0.005;//GeV^-2
  double t = t0 / (4.0 * M_PI);
  return t;
}

double E0(const double m, const double p){
  return sqrt(m*m + p*p);
}

double EnergyDenominator(const double q, const double p1, const double p, const double k){
  double result = q + NA * MN - E0((NA - 1) * MN, p1) - E0(Mmeson, k) - E0(MN, p);
  return result;
}

double getrelativep(const double mp, const double * p, const double mk, const double * k){
  //p: p, theta, phi
  double pk = p[0] * k[0] * (cos(p[1]) * cos(k[1]) + sin(p[1]) * sin(k[1]) * cos(p[2] - k[2]));
  double K2 = sqrt( (mp*mp + p[0]*p[0]) * (mk*mk + k[0]*k[0])) - pk;
  return sqrt( (K2 * K2 - mp*mp * mk*mk) / (2.0 * K2 + mp*mp + mk*mk));
}

int vector_add(const double * p, const double * k, double * r){
  //get vector p + k
  double pk = p[0] * k[0] * (cos(p[1]) * cos(k[1]) + sin(p[1]) * sin(k[1]) * cos(p[2] - k[2]));
  r[0] = sqrt(p[0]*p[0] + k[0]*k[0] + 2.0 * pk);
  r[1] = acos( (p[0] * cos(p[1]) + k[0] * cos(k[1])) / r[0]);
  double ry = p[0] * sin(p[1]) * sin(p[2]) + k[0] * sin(k[1]) * sin(k[2]);
  double rx = p[0] * sin(p[1]) * cos(p[2]) + k[0] * sin(k[1]) * cos(k[2]);
  double tr = 0;
  if (rx == 0 && ry >=0) r[2] = M_PI / 2.0;
  else if (rx == 0 && ry < 0) r[2] = -M_PI / 2.0;
  else {
    tr = ry / rx;
    if (rx > 0) r[2] = atan(tr);
    else r[2] = atan(tr) + M_PI;
  }
  if (r[2] > M_PI) r[2] = r[2] - 2.0 * M_PI;
  return 0;
}

int vector_sub(const double * p, const double * k, double * r){
  //get vector p - k
  double t[3];
  t[0] = k[0];
  t[1] = M_PI - k[1];
  t[2] = M_PI + k[2];
  vector_add(p, t, r);
  return 0;
}

double Tfi_den(const double * x){
  //x: q:0-2, pd:3-5, p:6-8, k:9-11
  double q[3] = {x[0], x[1], x[2]};
  double pd[3]= {x[3], x[4], x[5]};
  double p[3] = {x[6], x[7], x[8]};
  double k[3] = {x[9], x[10], x[11]};
  double p2[3];
  vector_sub(pd, k, p2);
  double tmp[3], p1[3];
  vector_sub(q, p, tmp);
  vector_sub(tmp, k, p1);
  double wf1 = wf_N_iso(p1[0], L1);
  double wf2 = wf_N_iso(p2[0], L2);
  double deno = EnergyDenominator(q[0], p1[0], p[0], k[0]);
  double qc = getrelativep(0.0, q, MN, p1);
  double tt = ampT(qc);
  double Q = getrelativep(Mmeson, k, MN, p2);//phi N relative
  double FF = ampF(Q);
  double result = wf1 * wf2 * FF / deno * tt;
  return result;
}

double integrand_Tfi(double * k, double * others){
  double x[12];
  for (int i = 0; i < 9; i++) x[i] = others[i];
  x[9] = tan(k[0]);
  x[10] = k[1]; x[11] = k[2];
  double TT = Tfi_den(x);
  return x[9] * x[9] * sin(x[10]) * TT / (cos(k[0]) * cos(k[0]));
}

double Tfi(double * par){
  //par: q:0-2, pd:3-5, p:6-8
  double q[3] = {par[0], par[1], par[2]};
  double pd[3]= {par[3], par[4], par[5]};
  double p[3] = {par[6], par[7], par[8]};
  double result;
  //Energy conservation check
  if (q[0] + MN - E0(Md, pd[0]) <= 0) result = 0;
  else {
    TF3 f("Integrand", integrand_Tfi, 0, M_PI/2.0, 0.0, M_PI, -M_PI, M_PI, 9);
    f.SetParameters(par);
    ROOT::Math::WrappedMultiTF1 wf1(f);
    ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetFunction(wf1);
    ig.SetRelTolerance(0.001);
    double xmin[] = {0, 0, -M_PI};
    double xmax[] = {M_PI/2.0, M_PI, M_PI};
    result = ig.Integral(xmin, xmax);
  }
  //std::cout << result <<  std::endl;
  return result;
}

double integrand_ds(double * x, size_t dim, void * params){
  double * getpar = (double *) params;
  double par[9];
  for (int i = 0; i < 6; i++) par[i] = getpar[i];
  double q[3] = {par[0], par[1], par[2]};
  double pd[3] = {par[3], par[4], par[5]};
  double result;
  if (q[0] + MN - E0(Md, pd[0]) <= 0) result = 0;
  else {
    par[6] = sqrt(pow(q[0] + 2.0*MN - E0(Md, pd[0]), 2) - MN*MN);
    par[7] = x[0];
    par[8] = x[1];
    result = pow(2.0*M_PI, 4) * pd[0] * pd[0] * E0(MN, par[6]) * par[6] * sin(par[7]) * pow(Tfi(par), 2);
  }
  return result;
}

double integrand_s(double * x, size_t dim, void * params){
  double * getpar = (double *) params;
  double par[9];
  for (int i = 0; i < 4; i++) par[i] = getpar[i];
  double pd[3] = {par[3], x[0], 2.0};
  double p[3] = {0.0, x[1], x[2]};
  double result;
  if (par[0] + MN - E0(Md, pd[0]) <= 0) result = 0;
  else {
    p[0] = sqrt(pow(par[0] + 2.0*MN - E0(Md, pd[0]), 2) - MN*MN);
    //x: pd, theta d, theta p, phi p
    for (int j = 0; j < 3; j++) par[j+3] = pd[j];
    for (int j = 0; j < 3; j++) par[j+6] = p[j];
    result = pow(2.0*M_PI, 5) * pd[0] * pd[0] * E0(MN, p[0]) * p[0] * sin(p[1]) * sin(pd[1]) * pow(Tfi(par), 2);
  }
  return result;
}
  
double dsigma(double * par){
  double res, err;
  double q[3] = {par[0], par[1], par[2]};
  double pd[3] = {par[3], par[4], par[5]};
  if (q[0] + MN - E0(Md, pd[0]) <= 0) res = 0;
  else {
    std::cout << "Monte carlo integrating..." << std::endl;
    double xl[2] = {0.0, -M_PI};
    double xu[2] = {M_PI, M_PI};
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_monte_function F;
    F.f = &integrand_ds;
    F.dim = 2;
    F.params = par;
    size_t calls = 1000;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_monte_vegas_state * s = gsl_monte_vegas_alloc(2);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, 500, r, s, &res, &err);
    //gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &res, &err);
    gsl_monte_vegas_free(s);
  }
  std::cout << res << "  " << err << std::endl;
  return res;
}

double sigma(double * q){
  double res, err;
  //double pdmax = sqrt( (q[0] + MN) * (q[0] * MN) - Md * Md );
  double xl[3] = {0.0, 0.0, -M_PI};
  double xu[3] = {M_PI, M_PI, M_PI};
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_monte_function F;
  F.f = &integrand_s;
  F.dim = 3;
  F.params = q;
  size_t calls = 500;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_monte_vegas_state * s = gsl_monte_vegas_alloc(3);
  std::cout << "Monte carlo integrating..." << std::endl;
  //gsl_monte_vegas_integrate(&F, xl, xu, 3, 100, r, s, &res, &err);
  gsl_monte_vegas_integrate(&F, xl, xu, 3, calls, r, s, &res, &err);
  gsl_monte_vegas_free(s);
  std::cout << res << "  " << err << std::endl;
  return res;
}

*/

#endif
