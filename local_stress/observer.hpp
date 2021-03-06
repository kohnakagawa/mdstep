#pragma once
#include "variables.hpp"
#include "LSCMD/ls_calculator.hpp"
//------------------------------------------------------------------------
class Observer {
public:
  double kinetic_energy(Variables *vars);
  double potential_energy(Variables *vars, std::vector<Pair> & pairs);
  double temperature(Variables *vars) {return kinetic_energy(vars) / 1.5;}
  double pressure(Variables *vars, std::vector<Pair> &pairs);
  double total_energy(Variables *vars, std::vector<Pair> &pairs) {return kinetic_energy(vars) + potential_energy(vars, pairs);}
  void local_pressure(Variables *vars);
  void local_density(Variables *vars);
};
//------------------------------------------------------------------------
