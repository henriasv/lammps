#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "region_ellipsoid.h"
#include "domain.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

RegEllipsoid::RegEllipsoid(class LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
	options(narg-8, &arg[8]);
	xmid = force->numeric(FLERR, arg[2]);
	ymid = force->numeric(FLERR, arg[3]);
	zmid = force->numeric(FLERR, arg[4]);
	a = force->numeric(FLERR, arg[5]);
	b = force->numeric(FLERR, arg[6]);
	c = force->numeric(FLERR, arg[7]);
}

RegEllipsoid::~RegEllipsoid()
{
	;
}

int RegEllipsoid::inside(double x, double y, double z)
{
	double dx = x-xmid;
	double dy = y-ymid;
	double dz = z-zmid;
	return ( (dx/a)*(dx/a)+(dy/b)*(dy/b)+(dz/c)*(dz/c) ) < 1.0;
}

int RegEllipsoid::surface_interior(double *x, double cutoff)
{
	return 0;
}

int RegEllipsoid::surface_exterior(double *x, double cutoff)
{
	return 0;
}
