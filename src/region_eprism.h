#ifdef REGION_CLASS
RegionStyle(eprism, RegEprism)
#else

#ifndef	LMP_REGION_EPRISM_H
#define	LMP_REGION_EPRISM_H

#include "region.h"

namespace LAMMPS_NS 
{

class RegEprism : public Region
{
 public:
 	RegEprism(class LAMMPS *, int, char **);
 	~RegEprism();
 	int inside(double, double, double);
 	int surface_interior(double *, double);
 	int surface_exterior(double *, double);

 private:
 	double a, b, xmid, ymid;
};

}

#endif
#endif