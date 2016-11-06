#include "Lsidis0.h"

using namespace std;

int main(){

  const double rad = M_PI / 180.0;

  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);
  TLorentzVector lp, Ph;

  //lp.SetXYZM(p_ele * sin(theta_ele) * cos(phi_ele), p_ele * sin(theta_ele) * sin(phi_ele), p_ele * cos(theta_ele), 0.000511);
  //Ph.SetXYZM(p_pion * sin(theta_pion) * cos(phi_pion), p_pion * sin(theta_pion) * sin(phi_pion), p_pion * cos(theta_pion), 0.1395);

  Lsidis mysidis;

  TLorentzVector boo(0.1, 0.15, -0.3, 1.0);
  //l.Boost(boo.BoostVector());
  //P.Boost(boo.BoostVector());

  mysidis.SetNucleus(2.0, 1.0);
  mysidis.SetHadron("pi+");
  mysidis.SetPDFset("CT14lo");
  mysidis.SetInitialState(l, P);

  //cout << mysidis.FUUT() << endl;

  //mysidis.Test();

  mysidis.SetVariables(0.192311, 0.665824, 0.359942, 0.290156, 2.00175, 0.334765);
  mysidis.CalculateFinalState();
  lp = mysidis.GetLorentzVector("lp");
  Ph = mysidis.GetLorentzVector("Ph");

  //cout << lp.E() << " " << lp.Theta() * 180.0 / M_PI << " " << lp.Phi() * 180 / M_PI << endl;
  cout << Ph.E() << " " << Ph.Theta() * 180.0 / M_PI << " " << Ph.Phi() * 180 / M_PI << endl;

  mysidis.SetFinalState(lp, Ph);
  mysidis.CalculateVariables();
  mysidis.Test();
  


  return 0;
}
