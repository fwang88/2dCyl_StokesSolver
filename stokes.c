#include <petscdmda.h>
#include <math.h>
#include "stokes.h"
#include <getopt.h>


levelset_vec create_levelset(PetscScalar maxR, 
			     PetscScalar maxZ, 
			     PetscScalar dr, 
			     PetscScalar dz,
			     PetscInt nghostlayer) {
  levelset_vec X;
  PetscInt nr, nz;
  DM da;
  
  nr = round( maxR/dr );
  nz = round( maxZ/dz );
  DMDACreate2d(PETSC_COMM_WORLD,
               DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
               DMDA_STENCIL_BOX,
               nr, nz,
               1, PETSC_DECIDE, 1,
               nghostlayer, NULL, NULL, &da);
  
  X.maxR = maxR;
  X.maxZ = maxZ;
  X.nr = nr;
  X.nz = nz;
  X.dr = dr;
  X.dz = dz;
  X.da = da;
  DMCreateGlobalVector(X.da,&X.data);
  VecSet(X.data,0.0);

  return X;
}

void destroy_levelset(levelset_vec X) {
  free(X.data);
}

/*******************************************************/

void initial_levelset(levelset_vec *G,
		      PetscScalar r0, PetscScalar pertb,
                      char *mode) {

  PetscInt llr, llz, lsizer, lsizez, rank;

  DMDAGetCorners(G->da, &llr, &llz, 0, &lsizer, &lsizez, 0);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscScalar interfacei, interfaceR;
  PetscScalar dr, dz, z, lambdaz;
  dr = G->dr;
  dz = G->dz;
  interfacei = r0/dr;

  PetscInt i, j, j0;
  PetscScalar number;
  PetscScalar **array;
  DMDAVecGetArray(G->da, G->data, &array);
  
  if( strcmp(mode,"periodic") == 0 ) {
    lambdaz = G->maxZ - dz/2;
    srandom(rank);
    
    for(j=llz; j<llz+lsizez; j++) {
      z = j * dz;
      interfaceR = interfacei + (pertb/dr) * cos(z/lambdaz*2*M_PI);
      for(i=llr; i<llr+lsizer; i++) {
	array[j][i] = (i-interfaceR) * dr;
      }
    }
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  else if( strcmp(mode,"end") == 0 ) {
    j0 = G->nz - 50;
    srandom(rank);
    printf("j0 = %d\n",j0);
    for( j=llz; j<llz+lsizez; j++ ) {
      if(j<j0) {
	for( i=llr; i<llr+lsizer; i++) {
	  array[j][i] = sqrt( (i)*(i)*dr*dr + (j-j0)*(j-j0)*dz*dz ) - r0;
	}
      }
      else {
	number = ((double) random() / (RAND_MAX)) - 0.5;
	for( i=llr; i<llr+lsizer; i++) {
	  interfaceR = interfacei + number * pertb;
	  array[j][i] = (i - interfaceR) * dr;
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


/**********************************************************************/

void get_input(int argc, char **args, 	       
	       PetscScalar *maxr, PetscScalar *maxz,
	       PetscScalar *dr, PetscScalar *dz, 
	       PetscScalar *r0, PetscScalar *pertb,
	       PetscInt *period_or_end,
	       PetscScalar *mui, PetscScalar *muo,
	       PetscScalar *vf,
	       PetscInt *temp_profile,
	       PetscScalar *tlow, PetscScalar *thigh, 
	       PetscScalar *twidth, PetscScalar *lowtwidth,
	       PetscScalar *restart, PetscScalar *trestart,
	       PetscScalar *outputdt
	       ) {
  
  int opt;
  
  struct option opts[] = {
    { "maxR", 1, NULL, 1 },
    { "maxZ", 1, NULL, 2 },
    { "dr",   1, NULL, 3 },
    { "dz",   1, NULL, 4 },
    { "r0",   1, NULL, 5 },
    { "pertb",1, NULL, 6 },
    { "period_or_end", 1, NULL, 7 },
    { "mui",  1, NULL, 8 },
    { "muo",  1, NULL, 9 },
    { "vf",   1, NULL, 10 },
    { "temp_profile",  1, NULL, 11 },
    { "tlow", 1, NULL, 12 },
    { "thigh",1, NULL, 13 },
    { "twidth",        1, NULL, 14 },
    { "lowtwidth",     1, NULL, 15 },
    { "restart",       1, NULL, 16 },
    { "trestart",      1, NULL, 17 },
    { "outputdt",      1, NULL, 18 },
    { NULL,   0, NULL, 0 }
  };

  while ((opt = getopt_long(argc, args, "h", opts, NULL)) != -1) {
    switch (opt) {
    case 1:
      *maxR = atof(optarg);
      break;
    case 2:
      *maxZ = atof(optarg);
      break;
    case 3:
      *dr   = atof(optarg);
      break;
    case 4:
      *dz   = atof(optarg);
      break;
    }
  }

}
