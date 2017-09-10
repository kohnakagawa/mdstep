#pragma once
//------------------------------------------------------------------------
const double Lx = 40, Ly = 40, Lz = 80;
const double dt = 0.005;
const double CUTOFF = 3.0;
const double MARGIN = 0.5;
const double ML2 = (CUTOFF + MARGIN) * (CUTOFF + MARGIN);
const double CL2 = (CUTOFF*CUTOFF);
const double RC2 = 1.0 / CL2;
const double RC6 = RC2 * RC2 * RC2;
const double RC12 = RC6 * RC6;
const double C0 = - 4.0 * (RC12 - RC6);
inline void adjust_periodic(double &dx, double &dy, double &dz) {
  const double LHx = Lx * 0.5;
  const double LHy = Ly * 0.5;
  const double LHz = Lz * 0.5;
  if (dx < -LHx)dx += Lx;
  if (dx > LHx) dx -= Lx;
  if (dy < -LHy)dy += Ly;
  if (dy > LHy) dy -= Ly;
  if (dz < -LHz)dz += Lz;
  if (dz > LHz) dz -= Lz;
}
//------------------------------------------------------------------------
