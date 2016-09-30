#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <omp.h>

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


int main(){
  double x;
#pragma omp parallel for
  for (int i = 0; i < 100; i++){
    x = i*i;
    cout << x << endl;
  }


  
  return 0;
}
  
