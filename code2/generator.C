#include "bound.h"

using namespace std;

int main(int argc, char * argv[]){

  LoadDS(argv[1]);
  genData(argv[2], 10000000);//ds grid file, create root file, Nsim

  return 0;
}


