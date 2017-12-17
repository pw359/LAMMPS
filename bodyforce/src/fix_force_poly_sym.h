/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(force/poly/sym,FixForcePolySym)

#else

#ifndef LMP_FIX_FORCE_POLY_SYM_H
#define LMP_FIX_FORCE_POLY_SYM_H

#include "fix.h"
#include <vector>


namespace LAMMPS_NS {

class FixForcePolySym : public Fix {
 public:
  FixForcePolySym(class LAMMPS *, int, char **);
  ~FixForcePolySym();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double memory_usage();

  double evaluate_polynomial(std::vector<double>&, const double);
  double evaluate_fit(const double);
  double dump_polynomials();

 private:

  double lb, ub, L, delta;

  int dimension;
  int degree;
  std::vector<double> coeffs, coeffs_der, coeffs_centre, coeffs_right;

  int symmetry;
  int ilevel_respa;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix addforce does not exist

Self-explanatory.

E: Variable name for fix addforce does not exist

Self-explanatory.

E: Variable for fix addforce is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix addforce

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix addforce

Must define an energy variable when applying a dynamic
force during minimization.

*/
