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

#ifdef DIHEDRAL_CLASS

DihedralStyle(opls/omp,DihedralOPLSOMP)

#else

#ifndef LMP_DIHEDRAL_OPLS_OMP_H
#define LMP_DIHEDRAL_OPLS_OMP_H

#include "stdio.h"
#include "dihedral_omp.h"

namespace LAMMPS_NS {

class DihedralOPLSOMP : public DihedralOMP {
 public:
  DihedralOPLSOMP(class LAMMPS *);
  ~DihedralOPLSOMP();
  virtual void compute(int, int);
  virtual void coeff(int, int, char **);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);

  virtual double memory_usage();

 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval();

private:
  double *k1,*k2,*k3,*k4;

  void allocate();
};

}

#endif
#endif
