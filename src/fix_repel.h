/* -*- c++ -*- ----------------------------------------------------------
  July 2015 Henrik Andersen Sveinsson
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

FixStyle(repel,FixRepel)

#else

#ifndef LMP_FIX_REPEL_H
#define LMP_FIX_REPEL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRepel : public Fix {
 public:
  FixRepel(class LAMMPS *, int, char **);
  ~FixRepel();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int istyle,scaleflag,side;
  char *xstr,*ystr,*zstr,*astr,*bstr,*cstr,*kstr;
  int xvar,yvar,zvar,avar,bvar,cvar,kvar;
  double xvalue,yvalue,zvalue,avalue,bvalue,cvalue,kvalue;
  int repeller_flag,planeside;
  double repeller[4],repeller_all[4];
  int cdim,varflag;
  int nlevels_respa;

  void options(int, char **);
};

}

#endif
#endif
