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
  PetscViewer viewer;

  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  
  PetscScalar maxr, maxz, dr, dz, r0;
  PetscScalar tension;
  PetscScalar pertb;
  PetscInt period_or_end;
  PetscScalar mui, muo;
  PetscScalar vf;
  PetscInt temp_profile;
  PetscScalar tlow, thigh, twidth, lowtwidth;
  PetscScalar restart, trestart;
  PetscScalar outputdt;
  char mode[10];
  PetscScalar *mu1=NULL, *mu2=NULL;

  get_input(argc, args, 
	    &maxr, &maxz, &dr, &dz, &r0,
	    &tension,
	    &pertb,
	    &period_or_end, &mui, &muo, &vf,
	    &temp_profile, &tlow, &thigh, &twidth,
	    &lowtwidth,
	    &restart, &trestart,
	    &outputdt,
	    mode,
	    mu1, mu2
	    );

  PetscInt nghostlayer = 3;
  G = create_levelset(maxr, maxz, dr, dz, nghostlayer);
  PetscObjectSetName((PetscObject)G.data, "levelset");

  PetscInt reinitstep =2;
  
  if(restart == 0) {
    initial_levelset(&G, r0, pertb, mode);
    output(&G, 0);
    reinitiate(&G, reinitstep);
    trestart = 0;
  }
  else {
    load_levelset(&G,trestart);
    G.t = trestart;
  }
  
  



















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
