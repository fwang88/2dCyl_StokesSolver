/* RK2DPeriod solving the level set equation with a dt time interval */

#include <petscksp.h>
#include <petscdmda.h>
#include <time.h>
#include <math.h>
#include "myheader.h"
//#include <mpi.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void RK2DPeriod(DM da, DM da3, Vec G, Vec uwp, PetscScalar dz, PetscScalar dr, PetscInt nz, PetscInt nr, PetscScalar dt) /* dz dr nz nr */
{
	Vec G1, G2;
	VecDuplicate(G,&G1);
	VecDuplicate(G,&G2);

	Euler_mex2D(da,da3,G1, G, uwp, dz, dr, nz, nr, dt);
	Euler_mex2D(da,da3,G2, G1, uwp, dz, dr, nz, nr, dt);
	VecCopy(G,G1);
	VecScale(G1,0.75);
	VecAXPY(G1,0.25,G2);
	VecSet(G2,0);
	Euler_mex2D(da,da3,G2, G1, uwp, dz, dr, nz, nr, dt);
	VecScale(G,1.0/3.0);
	VecAXPY(G,2.0/3.0,G2);
	VecDestroy(&G1);
	VecDestroy(&G2);
}

