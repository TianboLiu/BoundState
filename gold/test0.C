#include "Lgold.h"

using namespace std;

int main(){

  SetFunctions();

  cout << TF_BWd.Integral(TF_BWd.GetXmin(),TF_BWd.GetXmax()) << endl;

  return 0;
}
