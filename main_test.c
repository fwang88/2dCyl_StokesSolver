#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <math.h>
#include <time.h>
#include "stokes.h"
#include <stdio.h>
#include <petscviewerhdf5.h>
#include <string.h>

#include <getopt.h>
#include <sys/stat.h>

int main(int argc, char **args) {
  levelset_vec G;


  char mode[10] = "end";
  PetscScalar pertb = 0.02, r0 = 0.97;
  
  //  PetscScalar *array;
  PetscViewer viewer;

  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  
  PetscScalar maxr, maxz, dr, dz, r0, pertb;
  PetscInt period_or_end;
  PetscScalar mui, muo;
  PetscScalar vf;
  PetscInt temp_profile;
  PetscScalar tlow, thigh, twidth, lowtwidth;
  PetscScalar restart, trestart;
  PetscScalar outputdt;

  get_input(argc, args, 
	    &maxr, &maxz, &dr, &dz, &r0, &pertb,
	    &period_or_end, &mui, &muo, &vf,
	    &temp_profile, &tlow, &thigh, &twidth,
	    &lowtwidth,
	    &restart, &trestart,
	    &outputdt
	    );

  PetscInt nghostlayer = 3;
  G = create_levelset(maxR, maxZ, dr, dz, nghostlayer);
  PetscObjectSetName((PetscObject)G.data, "levelset");

  initial_levelset(&G, r0, pertb, mode);
  





















  /** check initial level set **/
  /*
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"./test_G",
                          FILE_MODE_WRITE,&viewer);
    VecView(G.data,viewer);
    PetscViewerDestroy(&viewer);
  */

  

//  VecGetArray(G.data,&array);

  destroy_levelset(G);
  //free(array);
  PetscFinalize();
  return 0;
}
