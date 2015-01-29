#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <math.h>
#include <time.h>
#include "stokes.h"
#include <stdio.h>
#include <petscviewerhdf5.h>
#include <string.h>
#include <stdlib.h>
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
  /*  
  int i;
  for(i=0; i<para->nz; i++) {
    mu->mu2[i] = (double) (10.0 - 2.0)*1.0/para->nz * i + 2;
    //mu->mu2[i] = 2.0;
    mu->mu1[i] = 1.0;
  }
  */

  G = create_levelset(para);
  PetscObjectSetName((PetscObject)G.data, "levelset");

  mk_dir(para);

  if(para->restart == 0) {
    initial_levelset(&G, para);
    output_levelset(&G, para);
    reinit(&G, para);
    para->trestart = 0;
  }
  else {
    load_levelset(&G, para);
    G.t = para->trestart;
  }

  
  stokes_matrix B;
  B = create_stokes_matrix(para);  
  assembly_stokes_matrix(&B);
 
  stokes_force force, uwp;
  force = create_stokes_force(para, B.da);  
  uwp = create_stokes_force(para, B.da);

  MYKSP solver;
  initialize_MYKSP(&solver);

  PetscScalar dt, t1, t2, t_end;
  PetscInt output_count = 1, evolve_count = 0;
  t1 = para->trestart;
  t2 = t1;
  t_end = para->maxz / para->vf * 5;


  while (t2 < t_end) {
    evolve_count++;

    get_B_from_G(&B, &G, para, mu);
    get_force_from_G(&force, &G, para);

    stokes_solve(&uwp, 
		 &solver, 
		 &B, &force, 
		 &G, 
		 para, mu);
    
    
    lab_frame_shift_w(&uwp, para->vf);
    
    dt = advection_timestep(&uwp, para);

    advection_evolve(&G, &uwp, para, dt);

    if( evolve_count % 5 == 0) 
      reinit(&G, para);
    
    t2 = t1 + dt;
    G.t = t2;
    uwp.t = t2;
    force.t = t2;

    //    printf("t1 = %g, t2 = %g, dt = %g\n", t1, t2, dt);

    if( t1 < output_count * para->outputdt + para->trestart
	&& t2 >= output_count * para->outputdt + para->trestart ) {
      G.t = output_count * para->outputdt + para->trestart;
      uwp.t = G.t;
      force.t = G.t;
      
      output_levelset(&G, para);
      output_uwp(&uwp, para);
      output_force(&force, para);

      output_count++;
    }

    t1 = t2;
  }




  /** free memory **/
  free(mu->mu1);
  free(mu->mu2);
  free(mu);
  free(para);

  VecDestroy(&G.data);
  VecDestroy(&force.data);
  VecDestroy(&uwp.data);
  MatDestroy(&B.data);
  KSPDestroy(&solver.ksp);

  PetscFinalize();
  return 0;
}
