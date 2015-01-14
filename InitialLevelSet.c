#include <petscksp.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
void InitialLevelSet(PetscScalar dr, PetscScalar dz, PetscInt nr, PetscInt nz, PetscScalar interfacei, PetscScalar r0, PetscScalar pertb, PetscScalar **array, PetscInt llr, PetscInt llz, PetscInt lsizer, PetscInt lsizez)
{
  PetscScalar r, number, gij, interfaceR, z;
  PetscInt i, j, rank;
  PetscScalar lambdaz;
  lambdaz = nz*dz - dz/2;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  srandom ( rank );
  for (j=llz;j<llz+lsizez;j++) {
    z = j * dz;
    interfaceR = interfacei + (pertb/dr) * cos(z/lambdaz * 2*M_PI);
    for(i=llr;i<llr+lsizer;i++) {
      array[j][i] = (i - interfaceR)*dr;
    } 
  }  
}
