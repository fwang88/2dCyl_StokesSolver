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

  viscosity *mu;
  parameter *para;
  mu = malloc(sizeof(viscosity));
  para = malloc(sizeof(parameter));
  
  get_input(argc, args, mu, para);
  
  G = create_levelset(para);
  PetscObjectSetName((PetscObject)G.data, "levelset");

  mk_dir(para);

  if(para->restart == 0) {
    initial_levelset(&G, para);
    output(&G, 0, para);
    //reinitiate(&G, para->reinitstep);
    para->trestart = 0;
    printf("%g\n",G.t);
  }
  else {
    load_levelset(&G, para);
    G.t = para->trestart;
  }
  
  //  VecView(G.data, PETSC_VIEWER_DEFAULT);
  



















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
