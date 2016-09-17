#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "TCanvas.h"
#include "TGraph2D.h"

#include "bound.h"

using namespace std;

int main(){

  LoadDS("ds1450.dat");
  cout << sigmatotal()*3.89379e5*132.0 << endl;
  
  if (true){
    TCanvas * c1 = new TCanvas("", "", 800, 600);
    ds2D.Draw("surf1");
    
    c1->Print("c1.pdf");
  }
  return 0;
}
