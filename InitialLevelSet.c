#include <petscksp.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
void InitialLevelSet(PetscScalar dr, PetscScalar dz, PetscInt nr, PetscInt nz, PetscScalar interfacei, PetscScalar r0, PetscScalar pertb, PetscScalar **array, PetscInt llr, PetscInt llz, PetscInt lsizer, PetscInt lsizez, char *mode)
{
  PetscScalar r, number, gij, interfaceR, z;
  PetscInt i, j, j0, rank;
  PetscScalar lambdaz;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if(strcmp(mode,"periodic")==0) {
    lambdaz = nz*dz - dz/2;
    
    srandom ( rank );
    for (j=llz;j<llz+lsizez;j++) {
      z = j * dz;
      interfaceR = interfacei + (pertb/dr) * cos(z/lambdaz * 2*M_PI);
      for(i=llr;i<llr+lsizer;i++) {
        array[j][i] = (i - interfaceR)*dr;
      } 
    }
    return;  
  }
  else if(strcmp(mode,"end")==0) {
    j0 = nz - 50;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    srandom ( rank );
    for (j=llz;j<llz+lsizez;j++) {
      if(j<j0) {
        for(i=llr;i<llr+lsizer;i++) {
          array[j][i] = sqrt( (i)*(i)*dr*dr + (j-j0)*(j-j0)*dz*dz ) - r0;
        }
      }
      else {
        number = ((double) random() / (RAND_MAX)) - 0.5;
        for(i=llr;i<llr+lsizer;i++) {
          interfaceR = interfacei + number * pertb;
          gij = i - interfaceR;
          array[j][i] = gij * dr;
        }
      }
    }
    return;
  }
  else {
    if(rank==0) printf("Please check the initial condition configuration: periodic or end \n");
    return;
  }

}
