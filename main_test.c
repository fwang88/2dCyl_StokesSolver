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
    reinit(&G, para);
    para->trestart = 0;
    // printf("%g\n",G.t);
  }
  else {
    load_levelset(&G, para);
    G.t = para->trestart;
  }

  PetscViewer viewer;
  //VecView(G.data, PETSC_VIEWER_DEFAULT);
  
  stokes_matrix B;
  B = create_stokes_matrix(para);
  assembly_stokes_matrix(&B);
  
  stokes_force force, uwp;
  force = create_stokes_force(para, B.da);  
  uwp = create_stokes_force(para, B.da);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  //VecView(force.data, PETSC_VIEWER_DEFAULT);

  /** time evolution of stokes equation **/
  
  //  MatView(B.data, PETSC_VIEWER_DEFAULT);
  //  get_force_from_G(&force, &G, para);

  while (0) {
    
    get_B_from_G(&B, &G, para, mu);  
    
    get_force_from_G(&force, &G, para);
    /*
    if(itstep >= 11) {
      stokes_solver_direct(&ksp, &B, &force, &uwp);
    }
    else {
      stokes_solver_precond(&ksp, &B, &force, &uwp);
    }
    
    */
    break;
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
