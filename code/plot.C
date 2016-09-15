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
  cout << sigmatotal() << endl;
  
  if (false){
    double x, y, z;
    ifstream infile("ds1450.dat");
    
    TGraph2D * g2 = new TGraph2D(2800);
    
    for (int i = 0; i < 2800; i++){
    infile >> x >> y >> z;
    if (z < 1.0e-60) z = 1.0e-60;
    g2->SetPoint(i, x, y, z);
    }
    infile.close();
    
    TCanvas * c1 = new TCanvas("", "", 800, 600);
    g2->Draw("surf1");
    
    c1->Print("c1.pdf");
  }
  return 0;
}
