#include <stdio.h>
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
  X.t = 0;
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
	       PetscScalar *r0, 
	       PetscScalar *tension, 
	       PetscScalar *pertb,
	       PetscInt *period_or_end,
	       PetscScalar *mui, PetscScalar *muo,
	       PetscScalar *vf,
	       PetscInt *temp_profile,
	       PetscScalar *tlow, PetscScalar *thigh, 
	       PetscScalar *twidth, PetscScalar *lowtwidth,
	       PetscScalar *restart, PetscScalar *trestart,
	       PetscScalar *outputdt,
	       char *mode,
	       PetscScalar *mu1, PetscScalar *mu2
	       ) {
  
  int opt;
  
  struct option opts[] = {
    { "maxr", 1, NULL, 1 },
    { "maxz", 1, NULL, 2 },
    { "dr",   1, NULL, 3 },
    { "dz",   1, NULL, 4 },
    { "r0",   1, NULL, 5 },
    { "tension", 1, NULL, 6 },
    { "pertb",1, NULL, 7 },
    { "period_or_end", 1, NULL, 8 },
    { "mui",  1, NULL, 9 },
    { "muo",  1, NULL, 10 },
    { "vf",   1, NULL, 11 },
    { "temp_profile",  1, NULL, 12 },
    { "tlow", 1, NULL, 13 },
    { "thigh",1, NULL, 14 },
    { "twidth",        1, NULL, 15 },
    { "lowtwidth",     1, NULL, 16 },
    { "restart",       1, NULL, 17 },
    { "trestart",      1, NULL, 18 },
    { "outputdt",      1, NULL, 19 },
    { NULL,   0, NULL, 0 }
  };

  while ((opt = getopt_long(argc, args, "h", opts, NULL)) != -1) {
    switch (opt) {
    case 1:
      *maxr = atof(optarg);
      break;
    case 2:
      *maxz = atof(optarg);
      break;
    case 3:
      *dr   = atof(optarg);
      break;
    case 4:
      *dz   = atof(optarg);
      break;
    case 5:
      *r0   = atof(optarg);
      break;
    case 6:
      *tension = atof(optarg);
      break;
    case 7:
      *pertb = atof(optarg);
      break;
    case 8:
      *period_or_end = atof(optarg);
      break;
    case 9:
      *mui = atof(optarg);
      break;
    case 10:
      *muo = atof(optarg);
      break;
    case 11:
      *vf = atof(optarg);
      break;
    case 12:
      *temp_profile = atof(optarg);
      break;
    case 13:
      *tlow = atof(optarg);
      break;
    case 14:
      *thigh = atof(optarg);
      break;
    case 15:
      *twidth = atof(optarg);
      break;
    case 16:
      *lowtwidth = atof(optarg);
      break;
    case 17:
      *restart = atof(optarg);
      break;
    case 18:
      *trestart = atof(optarg);
      break;
    case 19:
      *outputdt = atof(optarg);
      break;
    }
  }
  if(*period_or_end) strcpy(mode, "periodic");
  else strcpy(mode, "end");
  
  int i, nz;
  PetscScalar t_z, z_width;
  nz = round((*maxz) / (*dz));
  z_width = (*twidth)/(*dz);
  mu1 = malloc( nz * sizeof(PetscScalar) );
  mu2 = malloc( nz * sizeof(PetscScalar) );
  
  for(i=0; i<nz; i++) {
    if((*temp_profile) == 1) {
      t_z = ((*thigh)-(*tlow)) * tanh( (nz + 1.0*(*lowtwidth)/(*dz) - i) / (z_width) ) + (*tlow);
      mu2[i] = pow(10, 26909.0/(t_z+273)-7.2348)/1.0e3;
      mu1[i] = pow(10, 819.0/(t_z+273)-3.727)/1.0e3;
    }
    else {
      mu2[i] = *muo;
      mu1[i] = *mui;
    }
  }
  
  if((*temp_profile) == 1) *tension = *tension * 1000;

  printf("running simulation with \n\
maxr = %g, maxz = %g, dr = %g, dz = %g, r0 = %g\n\
tension = %g, pertb = %g\n", *maxr, *maxz, *dr, *dz, *r0, *tension, *pertb);
  
  return;
}
