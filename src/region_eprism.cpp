#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "region_eprism.h"
#include "domain.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

RegEprism::RegEprism(class LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
	options(narg-6, &arg[6]);
	xmid = force->numeric(FLERR, arg[2]);
	ymid = force->numeric(FLERR, arg[3]);
	a = force->numeric(FLERR, arg[4]);
	b = force->numeric(FLERR, arg[5]);
}

RegEprism::~RegEprism() 
{
	;
}

int RegEprism::inside(double x, double y, double z) 
{
	double dx = x-xmid;
	double dy = y-ymid;
	return ( (dx/a)*(dx/a)+(dy/b)*(dy/b) ) < 1.0;
}

int RegEprism::surface_interior(double *x, double cutoff)
{
	return 0;
}

int RegEprism::surface_exterior(double *x, double cutoff)
{
	return 0;
}
