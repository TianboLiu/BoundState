#include "Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./quasi-elastic <Ebeam> <Nsim>" << endl;
    return 0;
  }

  double Ebeam = atof(argv[1]);
  Long64_t Nsim = (Long64_t) atof(argv[2]);

  Initialize();
  TLorentzVector ki(0, 0, Ebeam, Ebeam);
  TLorentzVector kf[2];

  double weight = 0;

  for (Long64_t i = 0; i < Nsim; i++){

    weight = GENERATE::Event_ep_QE(&ki, kf);
    
    cout << weight << endl;
  }

  return 0;
}
