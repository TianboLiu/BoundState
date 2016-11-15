#include "spectra.h"

using namespace std;

int main(int argc, char * argv[]){
  char * control = argv[1];
  if (strcmp("set", control) == 0){
    makesetgold();
  }
  if (strcmp("read", control) == 0){
    readset(argv[2]);
  }
  if (strcmp("fit", control) == 0){
    readset(argv[2]);
    FitMissingMomentum();
  }
  if (strcmp("energy", control) == 0){
    readset(argv[2]);
    FitMissingEnergy();
  }
  return 0;
}
