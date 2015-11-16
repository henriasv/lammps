#ifdef REGION_CLASS
RegionStyle(ellipsoid, RegEllipsoid)
#else

#ifndef	LMP_REGION_ELLIPSOID_H
#define	LMP_REGION_ELLIPSOID_H

#include "region.h"

namespace LAMMPS_NS
{

class RegEllipsoid : public Region
{
 public:
 	RegEllipsoid(class LAMMPS *, int, char **);
 	~RegEllipsoid();
 	int inside(double, double, double);
 	int surface_interior(double *, double);
 	int surface_exterior(double *, double);

 private:
 	double a, b, c, xmid, ymid, zmid;
};

}

#endif
#endif
