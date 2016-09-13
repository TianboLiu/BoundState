#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>

#include "TSystem.h"
#include "Math/Functor.h"
#include "Math/ParamFunctor.h"
#include "Math/Factory.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/GaussIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

using namespace std;

double func(const double * x, const double * par){
  return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) * par[0] * par[1];
}

int main(){
  double par[2] = {1, 2};
  ROOT::Math::WrappedParamFunction<> wf(&func, 3, 2);
  wf.SetParameters(par);
  double a[3] = {0, 0, 0};
  double b[3] = {1, 1, 1};
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
  ig.SetFunction(wf);
  double re = ig.Integral(a, b);

  cout << re << endl;
  
  return 0;
}
  
