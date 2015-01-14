
#include <petscksp.h>
#include <petscdmda.h>
#include <time.h>
#include <math.h>
#include "myheader.h"

double timestep(DM da3, Vec uwp, PetscScalar dr, PetscScalar dz, PetscScalar vf) {
  StokesDOF **array;
  DMDAVecGetArray(da3,uwp,&array);
  PetscInt llr,llz,lsizer,lsizez,i,j;
  
  DMDAGetCorners(da3,&llr,&llz,0,&lsizer,&lsizez,0);
  double maxu_proc=0, maxw_proc=0;
  for(i=llz;i<llz+lsizez;i++) {
    for(j=llr;j<llr+lsizer;j++) {
      if(fabs(array[i][j].u)>maxu_proc) maxu_proc=fabs(array[i][j].u);
      if(fabs(array[i][j].w)>maxw_proc) maxw_proc=fabs(array[i][j].w);
    }
  }
  
  double dtu = dr/maxu_proc;
  double dtw = dz/maxw_proc;
  double dtmin = MIN(dtu,dtw);
  double dtresult;
  (void) MPI_Allreduce(&dtmin, &dtresult, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
  return ((dtresult*1.1) > (dr/vf) ? (dr/vf) : (dtresult*1.1));
}
