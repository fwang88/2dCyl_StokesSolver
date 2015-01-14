/* RK2DReinit solving the reinitialization equation with a dt time interval */

#include <petscksp.h>
#include <petscdmda.h>
#include <time.h>
#include <math.h>
#include "myheader.h"
#include <stdio.h>

#define MIN(a,b) (((a)<(b))?(a):(b))

void RK2DReinit(Vec G, PetscInt ReinitStep, PetscScalar dz, PetscScalar dr, PetscInt nz, PetscInt nr, DM da) /* dz dr nz nr */
{
  PetscScalar epsilon, dtau;
  PetscInt i;
  PetscInt llr, llz, lsizer, lsizez, rank;
  Vec G1, G2;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  epsilon = 2.0 * MIN(dz,dr);
  dtau = 0.1 * MIN(dz,dr);
  
  DMDAGetCorners(da,&llr,&llz,0,&lsizer,&lsizez,0);
  
  VecDuplicate(G, &G1);
  VecDuplicate(G, &G2);
  
  for(i=0; i<ReinitStep; i++) {
    Euler_Reinit_2D(G1, G, dtau, epsilon, dz, dr, nz, nr, da);
    Euler_Reinit_2D(G2, G1, dtau, epsilon, dz, dr, nz, nr, da);
    VecCopy(G,G1);
    VecScale(G1,0.75);
    VecAXPY(G1,0.25,G2);
    Euler_Reinit_2D(G2, G1, dtau, epsilon, dz, dr, nz, nr, da);
    VecScale(G,1.0/3.0);
    VecAXPY(G,2.0/3.0,G2);
  }
  VecDestroy(&G1);
  VecDestroy(&G2);
}
