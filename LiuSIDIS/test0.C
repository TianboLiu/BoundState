#include "Lsidis0.h"

using namespace std;

int main(){

  const double rad = M_PI / 180.0;

  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);
  TLorentzVector lp, Ph;
  double p_ele = 6.64219;
  double theta_ele = 10.3625 * rad;
  double phi_ele = -65.1553 * rad;
  double p_pion = 1.43295;
  double theta_pion = 10.2494 * rad;
  double phi_pion = -42.8162 * rad;

  lp.SetXYZM(p_ele * sin(theta_ele) * cos(phi_ele), p_ele * sin(theta_ele) * sin(phi_ele), p_ele * cos(theta_ele), 0.000511);
  Ph.SetXYZM(p_pion * sin(theta_pion) * cos(phi_pion), p_pion * sin(theta_pion) * sin(phi_pion), p_pion * cos(theta_pion), 0.1395);

  Lsidis mysidis;

  TLorentzVector boo(0.1, 0.15, -0.3, 1.0);
  //l.Boost(boo.BoostVector());
  //P.Boost(boo.BoostVector());

  mysidis.SetNucleus(2.0, 1.0);
  mysidis.SetHadron("pi+");
  mysidis.SetPDFset("CT14lo");
  mysidis.SetInitialState(l, P);
  //mysidis.SetFinalState(lp, Ph);
  //mysidis.CalculateVariables();

  //cout << mysidis.FUUT() << endl;

  //mysidis.Test();

  mysidis.SetVariables(0.291231, 0.396164, 0.330381, 0.600107, 0.1622, 1.13717);
  mysidis.CalculateFinalState();
  lp = mysidis.GetLorentzVector("Ph");
  //lp.Boost(-boo.BoostVector());

  cout << lp.E() << " " << lp.Theta() * 180.0 / M_PI << " " << lp.Phi() * 180 / M_PI << endl;


  return 0;
}
