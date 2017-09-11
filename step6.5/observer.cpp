//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include "systemparam.hpp"
#include "observer.hpp"
#include "ls_calculator.hpp"
//------------------------------------------------------------------------
double Observer::kinetic_energy(Variables *vars) {
  double k = 0;
  for (auto &a : vars->atoms) {
    k += a.px * a.px;
    k += a.py * a.py;
    k += a.pz * a.pz;
  }
  k /= static_cast<double>(vars->number_of_atoms());
  return k * 0.5;
};
//------------------------------------------------------------------------
double
Observer::potential_energy(Variables *vars, std::vector<Pair> &pairs) {
  double v = 0.0;
  const int pp = pairs.size();
  const int pn = vars->number_of_atoms();
  Atom *atoms = vars->atoms.data();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    v += 4.0 * (1.0 / r12 - 1.0 / r6) + C0;
  }
  v /= static_cast<double>(pn);
  return v;
}
//------------------------------------------------------------------------
double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs) {
  const double N = static_cast<double>(vars->number_of_atoms());
  const double V = Lx * Ly * Lz;
  const double T = temperature(vars);
  double phi = 0.0;
  const int pp = pairs.size();
  Atom *atoms = vars->atoms.data();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    phi += 48.0 / r12 - 24.0 / r6;
  }
  phi = phi / 3.0 / V;
  const double density = N / V;
  return density * T + phi;
}
//------------------------------------------------------------------------
void
Observer::local_pressure(Variables *vars, std::vector<Pair> &pairs) {
  static LS::LSCalculator<double> lscalculator({0.0, 0.0, 0.0}, {Lx, Ly, Lz},
                                               LS::BoundaryType::PERIODIC_XYZ,
                                               {1, 1, 160});
  Atom *atoms = vars->atoms.data();
  // kinetic term
  const int N = vars->number_of_atoms();
  for (int k = 0; k < N; k++) {
    lscalculator.calcLocalStressKin({atoms[k].qx, atoms[k].qy, atoms[k].qz},
                                    {atoms[k].px, atoms[k].py, atoms[k].pz},
                                    1.0);
  }

  // potential term
  const int pp = pairs.size();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2);
    lscalculator.calcLocalStressPot2({atoms[j].qx, atoms[j].qy, atoms[j].qz},
                                     {atoms[i].qx, atoms[i].qy, atoms[i].qz},
                                     {-df * dx, -df * dy, -df * dz},
                                     {df * dx, df * dy, df * dz});
  }

  lscalculator.nextStep();
}
//------------------------------------------------------------------------
