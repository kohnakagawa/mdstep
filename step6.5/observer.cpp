//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include "systemparam.hpp"
#include "observer.hpp"
#include <omp.h>
//------------------------------------------------------------------------
Observer::Observer() {
  lscalculators = LS::CalculatorFactory<double>::createOMP(omp_get_max_threads(),
                                                           {0.0, 0.0, 0.0}, {Lx, Ly, Lz},
                                                           LS::BoundaryType::PERIODIC_XYZ,
                                                           {1, 1, 120},
                                                           {"Kinetic", "LJ"});
}
//------------------------------------------------------------------------
Observer::~Observer() {
  LS::LSHelpers<double>::saveLocalStressDistOMP(lscalculators);
}
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
Observer::local_pressure(Variables *vars) {
  Atom *atoms = vars->atoms.data();
  const double mass = 1.0;
  // kinetic term
#pragma omp parallel
  {
    const auto tid = omp_get_thread_num();
    const int N = vars->number_of_atoms();
#pragma omp for nowait
    for (int k = 0; k < N; k++) {
      lscalculators[tid]->calcLocalStressKin({atoms[k].qx, atoms[k].qy, atoms[k].qz},
                                             {atoms[k].px, atoms[k].py, atoms[k].pz},
                                             mass,
                                             0);
    }
  }

  // potential term
#pragma omp parallel
  {
    const auto tid = omp_get_thread_num();
    const int N = vars->number_of_atoms();
    const int *neighbor_list = vars->neighbor_list.data();
    const int *i_position    = vars->i_position.data();
    const int *j_count       = vars->j_count.data();
#pragma omp for nowait
    for (int i = 0; i < N; i++) {
      const double qix = atoms[i].qx;
      const double qiy = atoms[i].qy;
      const double qiz = atoms[i].qz;
      const int ip = i_position[i];
      for (int k = 0; k < j_count[i]; k++) {
        const int j = neighbor_list[ip + k];
        const double qjx = atoms[j].qx;
        const double qjy = atoms[j].qy;
        const double qjz = atoms[j].qz;
        double dx = qjx - qix;
        double dy = qjy - qiy;
        double dz = qjz - qiz;
        adjust_periodic(dx, dy, dz);
        const double r2 = (dx * dx + dy * dy + dz * dz);
        if (r2 > CL2)continue;
        const double r6 = r2 * r2 * r2;
        const double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2);
        lscalculators[tid]->calcLocalStressPot2NoCheck({qix, qiy, qiz},
                                                       {qjx, qjy, qjz},
                                                       {df * dx, df * dy, df * dz},
                                                       {-df * dx, -df * dy, -df * dz},
                                                       1);
      }
    }
    lscalculators[tid]->nextStep();
  }

  LS::LSHelpers<double>::showPressureTotalOMP(lscalculators);
}
//------------------------------------------------------------------------
