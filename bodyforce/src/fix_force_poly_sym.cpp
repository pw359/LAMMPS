/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_force_poly_sym.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include <vector>
//#define DEBUG


using namespace LAMMPS_NS;
using namespace FixConst;


// symmetric or asymmetric force profile

enum {NONE, SYM, ANTISYM};


/* ---------------------------------------------------------------------- */

FixForcePolySym::FixForcePolySym(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  if (narg < 8) error->all(FLERR,"Illegal fix force/poly/sym command");


  // pw359: not sure about dynamic group
  dynamic_group_allow = 1;
  //scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 3;
  global_freq = 1;
  //extscalar = 1;
  //extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // fix ID group force/poly/sym {"x"|"y"|"z"} {"sym"|"antisym"|"none"}  {N a_0 a_1 ... a_N} 

  // parse dimension

  if (strcmp(arg[3], "x") == 0) {
    dimension = 0;
  }
  else if (strcmp(arg[3], "y") == 0) {
    dimension = 1;
  }
  else if (strcmp(arg[3], "z") == 0) {
    dimension = 2;
  }
  else {
    error->all(FLERR, "Illegal fix force/poly/sym command: unknown dimension");
  }

  // parse symmetry argument

  if (strcmp(arg[4], "sym") == 0) {
    symmetry = SYM;
  }
  else if (strcmp(arg[4], "antisym") == 0) {
    symmetry = ANTISYM;
  }
  else if (strcmp(arg[4], "none") == 0) {
    symmetry = NONE;
  }
  else {
    error->all(FLERR, "Illegal fix force/poly/sym command: unknown symmetry option"); 
  }

  // parse delta 

  delta = force->numeric(FLERR,arg[5]);
 
  // parse degree of polynomial  

  degree  = atoi(arg[6]);

  if (degree < 0) {
    error->all(FLERR, "Illegal fix force/poly/sym command: degree of polynomial cannot be negative");
  } 

  if (narg < 8 + degree) {
    error->all(FLERR, "Illegal fix force/poly/sym command: require degree+1 coefficients");
  }

  // parse coefficients of polynomial 

  for (int i=0; i < degree+1; i++) {
    coeffs.push_back(force->numeric(FLERR,arg[7+i]));
  }  

  for (int i=1; i < degree+1; i++) {
    coeffs_der.push_back(i*coeffs[i]);
  }  
 
  // output settings


#ifdef DEBUG
  printf("SETTINGS:\n");
  printf("===================\n");
  printf("     dimension: %s\n", (dimension==0) ? "x" : (dimension==1) ? "y" : (dimension==2) ? "z" : "");
  printf("      symmetry: %s\n", (symmetry==SYM) ? "sym" : (symmetry == ANTISYM) ? "antisym" : (symmetry == NONE) ? "none" : "");
  printf("         delta: %lf\n", delta);
  printf("        degree: %d\n", degree);
  printf("  coefficients: ");

  for (std::vector<double>::iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
    printf("\t%lf", *it);
  }
  printf("\n");

  printf("  coefficients for derivative: ");
  for (std::vector<double>::iterator it = coeffs_der.begin(); it != coeffs_der.end(); ++it) {
    printf("\t%lf", *it);
  }
  printf("\n");
#endif

}

/* ---------------------------------------------------------------------- */

FixForcePolySym::~FixForcePolySym()
{
   //printf("FixForcePolySym::~FixForcePolySym");


}

/* ---------------------------------------------------------------------- */

int FixForcePolySym::setmask()
{
  // pw359: not entirely sure about datamask_read and datamask_modify
  datamask_read = datamask_modify = 0;
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixForcePolySym::init()
{

  //printf("FixForcePolySym::init()");


  lb = domain->boxlo[dimension];
  ub = domain->boxhi[dimension];
  L  = ub-lb;

  // check that domain extents are compatitble with symmetry arguments

  // check for periodicity, dimensionality and box type  

  if (domain->nonperiodic > 0 || domain->dimension != 3 || domain->triclinic == 1) {
     error->all(FLERR,"currently you can only use fix poly/force/sym for a 3-dimensional, orthogonal, periodic system ");
  }

  if ( (symmetry == SYM || symmetry == ANTISYM) && fabs(lb+ub) > 0) {
      error->all(FLERR, "the simulation box needs to be centred around the origin to be compatible with the imposed symmetry");
  }



  coeffs_centre.clear();
  coeffs_right.clear();



  if (symmetry == NONE) {

    double f_r, f_l, df_r, df_l, a0, a1, a2,a3;

    // right polynomial 

    df_l = evaluate_polynomial(this->coeffs_der, ub-delta/2.0);
    df_r = evaluate_polynomial(this->coeffs_der, lb+delta/2.0);
    f_l  = evaluate_polynomial(this->coeffs,     ub-delta/2.0);
    f_r  = evaluate_polynomial(this->coeffs,     lb+delta/2.0);

    a0 = 1./8. * (delta * df_l - delta * df_r + 4.0 * f_l + 4.0 * f_r);
    a1 = - (delta * df_l + delta * df_r + 6.0 * f_l - 6.0 * f_r)/(4.0 * delta);
    a2 = - (df_l - df_r)/(2.0 * delta);
    a3 = - (-delta * df_l - delta * df_r - 2.0 * f_l + 2.0 * f_r)/(delta*delta*delta);

    coeffs_right.push_back(a0);
    coeffs_right.push_back(a1);
    coeffs_right.push_back(a2);
    coeffs_right.push_back(a3);
  }
  else if (symmetry == SYM) {

    double f_r, f_l, df_r, df_l, a0, a1, a2,a3;

    // central polynomial 

    df_r = evaluate_polynomial(this->coeffs_der, delta/2.0);
    df_l = -df_r;
    f_r  = evaluate_polynomial(this->coeffs,     delta/2.0);
    f_l  = f_r;

    a0 = 1./8. * (delta * df_l - delta * df_r + 4.0 * f_l + 4.0 * f_r);
    a1 = - (delta * df_l + delta * df_r + 6.0 * f_l - 6.0 * f_r)/(4.0 * delta);
    a2 = - (df_l - df_r)/(2.0 * delta);
    a3 = - (-delta * df_l - delta * df_r - 2.0 * f_l + 2.0 * f_r)/(delta*delta*delta);

    coeffs_centre.push_back(a0);
    coeffs_centre.push_back(a1);
    coeffs_centre.push_back(a2);
    coeffs_centre.push_back(a3);

    // right polynomial 

    df_l = evaluate_polynomial(this->coeffs_der, ub-delta/2.0);
    df_r = -df_l;
    f_l  = evaluate_polynomial(this->coeffs,     ub-delta/2.0);
    f_r  = f_l;

    a0 = 1./8. * (delta * df_l - delta * df_r + 4.0 * f_l + 4.0 * f_r);
    a1 = - (delta * df_l + delta * df_r + 6.0 * f_l - 6.0 * f_r)/(4.0 * delta);
    a2 = - (df_l - df_r)/(2.0 * delta);
    a3 = - (-delta * df_l - delta * df_r - 2.0 * f_l + 2.0 * f_r)/(delta*delta*delta);

    coeffs_right.push_back(a0);
    coeffs_right.push_back(a1);
    coeffs_right.push_back(a2);
    coeffs_right.push_back(a3);

  }


  else if (symmetry == ANTISYM) {

    double f_r, f_l, df_r, df_l, a0, a1, a2,a3;

    // central polynomial 

    df_r = evaluate_polynomial(this->coeffs_der, delta/2.0);
    df_l = df_r;
    f_r  = evaluate_polynomial(this->coeffs,     delta/2.0);
    f_l  = -f_r;

    a0 = 1./8. * (delta * df_l - delta * df_r + 4.0 * f_l + 4.0 * f_r);
    a1 = - (delta * df_l + delta * df_r + 6.0 * f_l - 6.0 * f_r)/(4.0 * delta);
    a2 = - (df_l - df_r)/(2.0 * delta);
    a3 = - (-delta * df_l - delta * df_r - 2.0 * f_l + 2.0 * f_r)/(delta*delta*delta);

    coeffs_centre.push_back(a0);
    coeffs_centre.push_back(a1);
    coeffs_centre.push_back(a2);
    coeffs_centre.push_back(a3);

    // right polynomial 

    df_l = evaluate_polynomial(this->coeffs_der, ub-delta/2.0);
    df_r = df_l;
    f_l  = evaluate_polynomial(this->coeffs,     ub-delta/2.0);
    f_r  = -f_l;

    a0 = 1./8. * (delta * df_l - delta * df_r + 4.0 * f_l + 4.0 * f_r);
    a1 = - (delta * df_l + delta * df_r + 6.0 * f_l - 6.0 * f_r)/(4.0 * delta);
    a2 = - (df_l - df_r)/(2.0 * delta);
    a3 = - (-delta * df_l - delta * df_r - 2.0 * f_l + 2.0 * f_r)/(delta*delta*delta);

    coeffs_right.push_back(a0);
    coeffs_right.push_back(a1);
    coeffs_right.push_back(a2);
    coeffs_right.push_back(a3);

  }

#ifdef DEBUG
  dump_polynomials();
#endif

}



/* ---------------------------------------------------------------------- */

void FixForcePolySym::setup(int vflag)
{
//  printf("FixForcePolySym::setup()");

  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixForcePolySym::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double xmapped[3];
 
  if (update->ntimestep % nevery) return;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xmapped[0] = x[i][0];
      xmapped[1] = x[i][1];
      xmapped[2] = x[i][2];
      domain->remap(xmapped);
      f[i][dimension] += evaluate_fit(xmapped[dimension]);
    }
  } 
}

/* ---------------------------------------------------------------------- */

void FixForcePolySym::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixForcePolySym::memory_usage()
{
  return 0.0;
}


/* ----------------------------------------------------------------------
   evaluate polynomial force
------------------------------------------------------------------------- */

double FixForcePolySym::evaluate_polynomial(std::vector<double>& coeffs, const double x) 
{
    double x_to_i = 1.0;
    double sum    = 0.0; 
    for (std::vector<double>::iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
      sum += *it * x_to_i;
      x_to_i *= x; 
    }
    return sum;
}


/* ----------------------------------------------------------------------
   evaluate fitted force
------------------------------------------------------------------------- */

double FixForcePolySym::evaluate_fit(const double x) 
{
    double value = 0;
    double y = fabs(x);

    if (y > ub) {
       // TODO: map x back into box and recompute y
       error->all(FLERR, "polynomial evaluation out of range");
    }

    if (symmetry == NONE) {
      if (y <= ub-delta/2.0) {
         value = evaluate_polynomial(coeffs, x);
      }
      else if (x >= ub-delta/2.0){
         value = evaluate_polynomial(coeffs_right, x-ub);
      }
      else {
         value = evaluate_polynomial(coeffs_right, x-lb);
      }
    }

    else if (symmetry == SYM) {
      if (y <= delta/2.0) {
         value = evaluate_polynomial(coeffs_centre, y);
      }
      else if (delta/2.0 <= y && y <= ub-delta/2.0) {
         value = evaluate_polynomial(coeffs, y);
      }
      else {
         value = evaluate_polynomial(coeffs_right, ub-y);
      }
    }

    else if (symmetry == ANTISYM) {
      if (y <= delta/2.0) {
         value = evaluate_polynomial(coeffs_centre, x);
      }
      else if (delta/2.0 <= x && x <= ub-delta/2.0) {
         value = evaluate_polynomial(coeffs, x);
      }
      else if (ub-delta/2.0 <= x) {
         value = -evaluate_polynomial(coeffs_right, ub-y);
      }
      else if (lb + delta/2.0 <= x && x <= -delta/2.0) {
         value = -evaluate_polynomial(coeffs, y);
      }
      else {
         value = evaluate_polynomial(coeffs_right, ub-y);
      }
    }

    return value;
}


/* ----------------------------------------------------------------------
    write force profile to file for debugging
------------------------------------------------------------------------ */


#ifdef DEBUG
double FixForcePolySym::dump_polynomials() 
{
   const int res = 1000;
   const double dx = L/res;
   double x;
   FILE * fp = fopen("force.log", "w");
   if (fp) {
     for (int i=0; i<res; i++) {
       x   = lb + (i+0.5) * dx;
       fprintf(fp, "%lf   %lf\n", x, evaluate_fit(x)); 
     }
     fclose(fp);
   }
}
#endif








