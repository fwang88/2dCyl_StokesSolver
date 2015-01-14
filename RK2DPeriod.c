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
	//VecSet(G1,0);
	//VecSet(G2,0);

	//MPI_Barrier(MPI_COMM_WORLD);
	//*dt = timestep(da3,uwp,dr,dz);
	//printf("dt=%f\n",*dt);
	//*dt = 1.0;
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

/****************************************************************************************
double timestep(DM da3, Vec uwp, PetscScalar dr, PetscScalar dz) {
	StokesDOF **array;
	Vec uwp2;
	//DMCreateGlobalVector(da3,&uwp2);
	VecDuplicate(uwp,&uwp2);
	VecCopy(uwp,uwp2);
	DMDAVecGetArray(da3,uwp2,&array);
	PetscInt llr,llz,lsizer,lsizez,i,j;
	double maxv;
	DMDAGetCorners(da3,&llr,&llz,0,&lsizer,&lsizez,0);
	for(i=llz;i<llz+lsizez;i++) 
		for(j=llr;j<llr+lsizer;j++) 
			array[i][j].p = 0;
	VecMax(uwp2,NULL,&maxv);
	printf("maxv=%f\n",maxv);
	return dr/maxv;
}
****************************************************************************************/
