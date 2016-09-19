#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "TCanvas.h"
#include "TGraph2D.h"

#include "bound.h"

using namespace std;

int main(){
  if (false){
    LoadDS("DS/ds1450ss.dat");
    double ss = sigmatotal()*3.89379e5;
    cout << ss << endl;
    LoadDS("DS/ds1450sp.dat");
    double sp = sigmatotal()*3.89379e5;
    cout << sp << endl;
    LoadDS("DS/ds1450ps.dat");
    double ps = sigmatotal()*3.89379e5;
    cout << ps << endl;
    LoadDS("DS/ds1450pp.dat");
    double pp = sigmatotal()*3.89379e5;
    cout << pp << endl;
    cout << 12*ss + 32*sp + 32*ps + 56*pp << endl;
  }
  
  

  if (true){
    LoadDS("DS/ds1450t.dat");
    cout << sigmatotal()*3.89379e5 << endl;
    TCanvas * c1 = new TCanvas("", "", 800, 600);
    ds2D.Draw("surf1");
    
    c1->Print("c1.pdf");
  }
  return 0;
}
