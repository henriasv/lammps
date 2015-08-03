/*
July 2015 Henrik Andersen Sveinsson
*/
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

/* ----------------------------------------------------------------------
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_repel.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,ELLIPSOID};
enum{INSIDE,OUTSIDE};

/* ---------------------------------------------------------------------- */

FixRepel::FixRepel(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix repel command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    kstr = new char[n];
    strcpy(kstr,&arg[3][2]);
  } else kvalue = force->numeric(FLERR,arg[3]);

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors to geometry

  if (istyle == ELLIPSOID) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;
    if (!astr) avalue *= xscale;
    if (!bstr) bvalue *= yscale;
    if (!cstr) cvalue *= zscale;
    } else error->all(FLERR,"Illegal fix indent command");

  varflag = 0;
  if (xstr || ystr || zstr || astr || bstr || cstr) varflag = 1;
  repeller_flag = 0;
  repeller[0] = repeller[1] = repeller[2] = repeller[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixRepel::~FixRepel()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] astr;
  delete [] bstr;
  delete [] cstr;
}

/* ---------------------------------------------------------------------- */

int FixRepel::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRepel::init()
{
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix repel does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix repel is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix repel does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix repel is not equal style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix repel does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix repel is not equal style");
  }
  if (astr) {
    avar = input->variable->find(astr);
    if (avar < 0)
      error->all(FLERR, "Variable name for fix repel does nor exist");
    if (!input->variable->equalstyle(avar))
      error->all(FLERR, "Variable for fix repel is not equal style");
  }
  if (bstr) {
    bvar = input->variable->find(bstr);
    if (avar < 0)
      error->all(FLERR, "Variable name for fix repel does nor exist");
    if (!input->variable->equalstyle(avar))
      error->all(FLERR, "Variable for fix repel is not equal style");
  }
  if (cstr) {
    cvar = input->variable->find(cstr);
    if (avar < 0)
      error->all(FLERR, "Variable name for fix repel does nor exist");
    if (!input->variable->equalstyle(avar))
      error->all(FLERR, "Variable for fix repel is not equal style");
  }
  if (kstr) {
    kvar = input->variable->find(kstr);
    if (kvar < 0)
      error->all(FLERR, "Variable name for fix repel does nor exist");
    if (!input->variable->equalstyle(kvar))
      error->all(FLERR, "Variable for fix repel is not equal style");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRepel::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixRepel::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRepel::post_force(int vflag)
{
  // repeller values, 0 = energy, 1-3 = force components
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  repeller_flag = 0;
  repeller[0] = repeller[1] = repeller[2] = repeller[3] = 0.0;

  // Ellipsoidal repeller

  if (istyle == ELLIPSOID) {

    // ctr = current repeller center
    // remap into periodic box

    double ctr[3];
    double a, b, c,k;
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    if (astr) a = input->variable->compute_equal(avar);
    else a = avalue;
    if (bstr) b = input->variable->compute_equal(bvar);
    else b = bvalue;
    if (cstr) c = input->variable->compute_equal(cvar);
    else c = cvalue;
    if (kstr) k = input->variable->compute_equal(kvar);
    else k = kvalue;
    domain->remap(ctr);

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag,fx,fy,fz;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx,dely,delz);
        fmag = 2*k*exp(-delx*delx/a/a-dely*dely/b/b-delz*delz/c/c);
        dr = sqrt(delx*delx/a/a+dely*dely/b/b+delz*delz/c/c);
        if (dr < 1)
        {
          fx = fmag*delx/a/a;
          fy = fmag*dely/b/b;
          fz = fmag*delz/c/c;
        }
        else
        {
          fx = fy = fz = 0;
        }
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          repeller[0] -= k * dr*dr*dr;
          repeller[1] -= fx;
          repeller[2] -= fy;
          repeller[3] -= fz;
      }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixRepel::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRepel::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of repeller interaction
------------------------------------------------------------------------- */

double FixRepel::compute_scalar()
{
  // only sum across procs one time

  if (repeller_flag == 0) {
    MPI_Allreduce(repeller,repeller_all,4,MPI_DOUBLE,MPI_SUM,world);
    repeller_flag = 1;
  }
  return repeller_all[0];
}

/* ----------------------------------------------------------------------
   components of force on repeller
------------------------------------------------------------------------- */

double FixRepel::compute_vector(int n)
{
  // only sum across procs one time

  if (repeller_flag == 0) {
    MPI_Allreduce(repeller,repeller_all,4,MPI_DOUBLE,MPI_SUM,world);
    repeller_flag = 1;
  }
  return repeller_all[n+1];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixRepel::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix repel command");

  istyle = NONE;
  xstr = ystr = zstr = astr = bstr = cstr = NULL;
  xvalue = yvalue = zvalue = avalue = bvalue = cvalue = kvalue = 0.0;
  scaleflag = 0;
  side = INSIDE;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ellipsoid") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal fix repel command");

      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        astr = new char[n];
        strcpy(astr,&arg[iarg+4][2]);
      } else avalue = force->numeric(FLERR,arg[iarg+4]);
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        int n = strlen(&arg[iarg+5][2]) + 1;
        bstr = new char[n];
        strcpy(bstr,&arg[iarg+5][2]);
      } else bvalue = force->numeric(FLERR,arg[iarg+5]);
      if (strstr(arg[iarg+6],"v_") == arg[iarg+6]) {
        int n = strlen(&arg[iarg+6][2]) + 1;
        cstr = new char[n];
        strcpy(cstr,&arg[iarg+6][2]);
      } else cvalue = force->numeric(FLERR,arg[iarg+6]);

      istyle = ELLIPSOID;
      iarg += 7;
    }
      else if (strcmp(arg[iarg],"units") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix repel command");
        if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
        else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
        else error->all(FLERR,"Illegal fix repel command");
        iarg += 2;

      } else if (strcmp(arg[iarg],"side") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix repel command");
        if (strcmp(arg[iarg+1],"in") == 0) side = INSIDE;
        else if (strcmp(arg[iarg+1],"out") == 0) side = OUTSIDE;
        else error->all(FLERR,"Illegal fix repel command");
        iarg += 2;
      } else error->all(FLERR,"Illegal fix repel command");
    }
    std::cout << "xvalue " << xvalue << std::endl;
    std::cout << "yvalue " << yvalue << std::endl;
    std::cout << "zvalue " << zvalue << std::endl;
    std::cout << "avalue " << avalue << std::endl;
    std::cout << "bvalue " << bvalue << std::endl;
    std::cout << "cvalue " << cvalue << std::endl;
}
