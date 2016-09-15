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
#include "TFile.h"
#include "TTree.h"

using namespace std;

double func(const double * x, const double * par){
  return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) * par[0] * par[1];
}

int main(){
  TLorentzVector p(0.1, 0.2, 0.3, 0.4);
  TFile * fs = new TFile("t.root", "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  Ts->Branch("p", "TLorentzVector", &p);
  Ts->Fill();
  fs->Write();


  
  return 0;
}
  
