#include <stdio.h>
#include <petscdmda.h>
#include <math.h>
#include "stokes.h"
#include <getopt.h>

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif


levelset_vec create_levelset(parameter *para) {
  levelset_vec X;
  PetscInt nr, nz;
  DM da;
  
  nr = round( para->maxr/para->dr );
  nz = round( para->maxz/para->dz );
  DMDACreate2d(PETSC_COMM_WORLD,
               DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
               DMDA_STENCIL_BOX,
               nr, nz,
               1, PETSC_DECIDE, 1,
               para->nghostlayer, NULL, NULL, &da);
  
  X.nr = nr;
  X.nz = nz;
  X.da = da;
  X.t = 0;
  DMCreateGlobalVector(X.da,&X.data);
  VecSet(X.data,0.0);
  
  return X;
}

void destroy_levelset(levelset_vec X) {
  free(X.data);
}


stokes_matrix create_stokes_matrix(parameter *para) {
  
  stokes_matrix B;
  DM da;
  int dof = 3;
  
  DMDACreate2d(PETSC_COMM_WORLD, 
	       DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
	       DMDA_STENCIL_BOX,
	       para->nr, para->nz,
	       1, PETSC_DECIDE, dof,
	       1, NULL, NULL, &da);
  
  DMSetMatrixPreallocateOnly(da, PETSC_TRUE);
  DMSetUp(da);

  B.nz = para->nz;
  B.nr = para->nr;
  B.da = da;
  
  return B;
}

void assembly_stokes_matrix(stokes_matrix *B) {

  DMCreateMatrix(B->da, MATAIJ, &B->data);
  
  MatAssemblyBegin(B->data, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B->data, MAT_FINAL_ASSEMBLY);

}


stokes_force create_stokes_force(parameter *para, DM da) {
  stokes_force X;

  X.nr = para->nr;
  X.nz = para->nz;
  X.da = da;

  DMCreateGlobalVector(X.da, &X.data);
  VecSet(X.data, 0.1);
  
  return X;
  
}


/*******************************************************/

void initial_levelset(levelset_vec *G, parameter *para) {
  
  PetscInt llr, llz, lsizer, lsizez, rank;

  DMDAGetCorners(G->da, &llr, &llz, 0, &lsizer, &lsizez, 0);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscScalar interfacei, interfaceR;
  PetscScalar z, lambdaz;
  
  interfacei = para->r0/para->dr;
  
  PetscInt i, j, j0;
  PetscScalar number;
  PetscScalar **array;
  DMDAVecGetArray(G->da, G->data, &array);
  
  if( strcmp(para->mode,"periodic") == 0 ) {
    lambdaz = para->maxz - para->dz/2;
    srandom(rank);
    
    for(j=llz; j<llz+lsizez; j++) {
      z = j * para->dz;
      interfaceR = interfacei + (para->pertb/para->dr) * cos(z/lambdaz*2*M_PI);
      for(i=llr; i<llr+lsizer; i++) {
	array[j][i] = (i-interfaceR) * para->dr;
      }
    }
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  else if( strcmp(para->mode,"end") == 0 ) {
    j0 = G->nz - 50;
    srandom(rank);
    printf("j0 = %d\n",j0);
    for( j=llz; j<llz+lsizez; j++ ) {
      if(j<j0) {
	for( i=llr; i<llr+lsizer; i++) {
	  array[j][i] = sqrt( (i)*(i)*para->dr*para->dr + (j-j0)*(j-j0)*para->dz*para->dz ) - para->r0;
	}
      }
      else {
	number = ((double) random() / (RAND_MAX)) - 0.5;
	for( i=llr; i<llr+lsizer; i++) {
	  interfaceR = interfacei + number * para->pertb;
	  array[j][i] = (i - interfaceR) * para->dr;
	}
      }
    }
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  else {
    if(rank==0) printf("Please check the initial condition configuration: periodic or end \n");
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  
  DMDAVecRestoreArray(G->da, G->data, &array);
  
  
  return;
}


/****************************************************************************/

void get_input(int argc, char **args, viscosity *mu, parameter *para) {
  
  int opt;
  
  struct option opts[] = {
    { "maxr",          1, NULL, 1 },
    { "maxz",          1, NULL, 2 },
    { "dr",            1, NULL, 3 },
    { "dz",            1, NULL, 4 },
    { "r0",            1, NULL, 5 },
    { "tension",       1, NULL, 6 },
    { "pertb",         1, NULL, 7 },
    { "period_or_end", 1, NULL, 8 },
    { "mui",           1, NULL, 9 },
    { "muo",           1, NULL, 10 },
    { "vf",            1, NULL, 11 },
    { "temp_profile",  1, NULL, 12 },
    { "tlow",          1, NULL, 13 },
    { "thigh",         1, NULL, 14 },
    { "twidth",        1, NULL, 15 },
    { "lowtwidth",     1, NULL, 16 },
    { "restart",       1, NULL, 17 },
    { "trestart",      1, NULL, 18 },
    { "outputdt",      1, NULL, 19 },
    { "nghostlayer",   1, NULL, 20 },
    { "epsilon",       1, NULL, 21 },
    { "reinitstep",    1, NULL, 22 },
     { NULL,   0, NULL, 0 }
  };

  while ((opt = getopt_long(argc, args, "h", opts, NULL)) != -1) {
    switch (opt) {
    case 1:
      para->maxr = atof(optarg);
      break;
    case 2:
      para->maxz = atof(optarg);
      break;
    case 3:
      para->dr = atof(optarg);
      break;
    case 4:
      para->dz = atof(optarg);
      break;
    case 5:
      para->r0 = atof(optarg);
      break;
    case 6:
      para->tension = atof(optarg);
      break;
    case 7:
      para->pertb = atof(optarg);
      break;
    case 8:
      para->period_or_end = atof(optarg);
      break;
    case 9:
      para->mui = atof(optarg);
      break;
    case 10:
      para->muo = atof(optarg);
      break;
    case 11:
      para->vf = atof(optarg);
      break;
    case 12:
      para->temp_profile = atof(optarg);
      break;
    case 13:
      para->tlow = atof(optarg);
      break;
    case 14:
      para->thigh = atof(optarg);
      break;
    case 15:
      para->twidth = atof(optarg);
      break;
    case 16:
      para->lowtwidth  = atof(optarg);
      break;
    case 17:
      para->restart = atof(optarg);
      break;
    case 18:
      para->trestart = atof(optarg);
      break;
    case 19:
      para->outputdt = atof(optarg);
      break;
    case 20:
      para->nghostlayer = atof(optarg);
      break;
    case 21:
      para->epsilon = atof(optarg);
      break;
    case 22:
      para->reinitstep = atof(optarg);
      break;
    }
  }

  para->nr = (int)(para->maxr/para->dr);
  para->nz = (int)(para->maxz/para->dz);

  if(para->period_or_end) strcpy(para->mode, "periodic");
  else strcpy(para->mode, "end");

  int i, nz;
  PetscScalar t_z, z_width;
  nz = para->nz;
  z_width = para->twidth / para->dz;
  mu->mu1 = malloc( nz * sizeof(PetscScalar) );
  mu->mu2 = malloc( nz * sizeof(PetscScalar) );

  for(i=0; i<nz; i++) {
    if(para->temp_profile == 1) {
      t_z = (para->thigh - para->tlow) * tanh( (nz + 1.0*para->lowtwidth/para->dz - i) / z_width ) + para->tlow;
      mu->mu2[i] = pow(10, 26909.0/(t_z+273)-7.2348)/1.0e3;
      mu->mu1[i] = pow(10, 819.0/(t_z+273)-3.727)/1.0e3;
    }
    else {
      mu->mu2[i] = para->muo;
      mu->mu1[i] = para->mui;
    }
  }
  
  if(para->temp_profile == 1) para->tension = para->tension * 1000;
  
  printf("running simulation with maxr = %g, maxz = %g, nr = %d, nz = %d, r0 = %g, tension = %g, pertb = %g, mode = %s, mu1 = %g, mu2 = %g\n", para->maxr, para->maxz, para->nr, para->nz, para->r0, para->tension, para->pertb, para->mode, mu->mu1[0], mu->mu2[0]);
  
  return;
}
