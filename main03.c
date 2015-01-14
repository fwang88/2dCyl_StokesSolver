#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <math.h>
#include <time.h>
#include "myheader.h"
#include <stdio.h>
#include <petscviewerhdf5.h>

#include <getopt.h>
#include <sys/stat.h>

int main(int argc, char **args) {
  double begin, end;
  Vec G;
  PetscInt nr, nz, i,j;
  PetscScalar muRat, vf, r0, tension, maxR, maxZ, dr, dz, interfacei, pertb, epsilon;
  PetscMPIInt size, rank, llr, llz, lsizer, lsizez, nghostlayer=3;
  DM da, da3;
  DMDABoundaryType br = DMDA_BOUNDARY_NONE, bz = DMDA_BOUNDARY_NONE;
  DMDAStencilType stype = DMDA_STENCIL_BOX;
  PetscViewer    viewer;
  KSP ksp;
  PC pc;
  
  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  r0 = 0.97;
  tension = 1.5; 
  muRat = 0.91;
  maxR = 10;
  maxZ = 10;
  dr = 0.1;
  dz = 0.1;
  pertb = 1.0;
  int hdf5 = 1;
  vf = 0.0;
  PetscScalar dt=0.001;
  PetscScalar outputstep=dt;
  epsilon = 0.3;

  struct option opts[] = {
    { "r0", 1, NULL, 1 },
    { "muRat", 1, NULL, 2 },
    { "tension", 1, NULL, 3 },
    { "maxR", 1, NULL, 4 },
    { "maxZ", 1, NULL, 5 },
    { "dr", 1, NULL, 6 },
    { "dz", 1, NULL, 7 },
    { "pertb", 1, NULL, 8 },
    { "dt", 1, NULL, 9 },
    { "outputstep", 1, NULL, 10},
    { "epsilon", 1, NULL, 11},
    { "hdf5", 0, NULL, 'h' },
    { NULL, 0, NULL, 0 }
  };
  int opt;
  while ((opt = getopt_long(argc, args, "h", opts, NULL)) != -1) {
    switch (opt) {
    case 1:
      r0 = atof(optarg);
      break;
    case 2:
      muRat = atof(optarg);
      break;
    case 3:
      tension = atof(optarg);
      break;
    case 4:
      maxR = atof(optarg);
      break;
    case 5:
      maxZ = atof(optarg);
      break;
    case 6:
      dr = atof(optarg);
      break;
    case 7:
      dz = atof(optarg);
      break;
    case 8:
      pertb = atof(optarg);
      break;
    case 9:
      dt = atof(optarg);
      break;
    case 10:
      outputstep = atof(optarg);
      break;
    case 11:
      epsilon = atof(optarg);
      break;
    case 'h':
      hdf5 = 1;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "unrecognized option");
    }
  }
  
  interfacei = r0/dr;
  nr = round(maxR/dr);
  nz = round(maxZ/dz);
  
  dt = dt * (dr/0.1);
  outputstep = outputstep * (dr/0.1);
  
  DMDACreate2d(PETSC_COMM_WORLD, br, bz, stype, nr, nz, 1, PETSC_DECIDE, 1,nghostlayer,NULL,NULL,&da);
  DMCreateGlobalVector(da,&G);
  VecSet(G,0);
  DMDAGetCorners(da,&llr,&llz,0,&lsizer,&lsizez,0);
  PetscObjectSetName((PetscObject)G, "levelset");
  PetscScalar **array;
  DMDAVecGetArray(da,G,&array);
  MPI_Barrier(MPI_COMM_WORLD);
  InitialLevelSet(dr, dz, nr, nz, interfacei, r0, pertb, array, llr,llz,lsizer,lsizez); 
  DMDAVecRestoreArray(da,G,&array);
  
  PetscScalar t1, t2, z_turnh, z_turnl, z_width, Tatz, mu1[nz], mu2[nz];
  char buffer[200];
  PetscInt count, itsteps, its;
  
  t1 = 0; t2 = t1+dt; 
  count = 1;
  
  for(i=0;i<nz;i++) {
    mu2[i] = 1.0;
    mu1[i] = mu2[i]*muRat;
  }
  
  Mat B, F; 
  Vec Force, uwp;
  PetscInt dof=3;
  nghostlayer = 1;
  DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, nr, nz, 1, PETSC_DECIDE, dof, nghostlayer,NULL,NULL,&da3);
  DMCreateGlobalVector(da3,&Force);
  VecSet(Force,0);
  VecDuplicate(Force,&uwp);
  PetscObjectSetName((PetscObject)uwp,"velocity");
  PetscObjectSetName((PetscObject)Force,"force");
  DMSetMatrixPreallocateOnly(da3, PETSC_TRUE);
  DMSetUp(da3);
  DMCreateMatrix(da3,MATAIJ,&B);
  DMCreateMatrix(da3,MATAIJ,&F);

  
  itsteps = 0;
  its = 51;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  DMDAGetCorners(da3,&llr,&llz,0,&lsizer,&lsizez,0);
  sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f",maxZ,r0,dr,dz,pertb,muRat);
  int status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  /**** output initial G ***/
  sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f/initialG.h5",maxZ,r0,dr,dz,pertb,muRat);
  PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
  VecView(G,viewer);
  PetscViewerDestroy(&viewer);

  /*** Reinitiate G and output ***/
  PetscInt ReinitStep = 2;
  MPI_Barrier(MPI_COMM_WORLD);
  RK2DReinit(G, ReinitStep, dz, dr, nz, nr, da);
  sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f/intial_after_reinit_G.h5",maxZ,r0,dr,dz,pertb,muRat);
  PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
  VecView(G,viewer);
  PetscViewerDestroy(&viewer);
  /*
  if (rank == 0) {
    printf("running maxR= %g x maxZ= %g simulation   \
    //with muRat = %g, epsilon= %g, pertb = %g, \n", maxR, maxZ, muRat, epsilon, pertb);
  }
  */

  /***************** time evolution of stokes equation ********************/
  while(count<=5) {
    
    if (hdf5) {
      if( t1<count*outputstep && t2>=count*outputstep) {
#ifdef PETSC_HAVE_HDF5
        sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f/outputG_t=%.6f.h5",maxZ,r0,dr,dz,pertb,muRat,t1);
        PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
        VecView(G,viewer);
        PetscViewerDestroy(&viewer);
        /*
          sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_pertb.%.3f_muRat.%.3f/uwp_t=%.6f.h5",maxZ,r0,dr,pertb,muRat,t2);
        PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
        VecView(uwp,viewer);
        PetscViewerDestroy(&viewer);        
        */
        count++;
#else
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "not compiled with HDF5");
#endif
      }
    }

    /********************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    itsteps++;
    AssembleB(B,da,da3,epsilon,nr,nz,dr,dz,G,mu1,mu2);
    AssembleForce(Force,G,da,da3,nr,nz,dr,dz,tension,epsilon);

    /* PetscObjectSetName((PetscObject)B,"matrix"); */
    /* 
    sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f/matrixB_t=%.6f.h5",maxZ,r0,dr,dz,pertb,muRat,t1);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,buffer,&viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(B,viewer);
    PetscViewerDestroy(&viewer);
    
    sprintf(buffer,"./wavelength.%.0f_r0.%.3f_dr.%.3f_dz.%.3f_pertb.%.3f_muRat.%.3f/Force_t=%.6f.h5",maxZ,r0,dr,dz,pertb,muRat,t1);
    PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
    VecView(Force,viewer);
    PetscViewerDestroy(&viewer);
    */
    
    if(its>=0) {
      KSPSetOperators(ksp,B,B,DIFFERENT_NONZERO_PATTERN);
      KSPSetDM(ksp,da3);
      KSPSetDMActive(ksp,PETSC_FALSE);
      KSPSetFromOptions(ksp);
      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
      KSPSetUp(ksp);

      // PCFactorGetMatrix(pc,&F);
      // MatMumpsSetIcntl(F,14,1000);

      KSPSolve(ksp,Force,uwp);
      its = 0;
    }
    else {
      KSPSetOperators(ksp,B,B,SAME_PRECONDITIONER);
      KSPSetUp(ksp);
      KSPSolve(ksp,Force,uwp);
      KSPGetIterationNumber(ksp,&its);
    }
   
    RK2DPeriod(da,da3,G,uwp,dz,dr,nz,nr,dt);
    
    if( itsteps % 1 == 0 ) {
      MPI_Barrier(MPI_COMM_WORLD);
      RK2DReinit(G,2,dz,dr,nz,nr,da);
    }
    t1 = t2;
    t2 = t1+dt;

  end = MPI_Wtime();

  if(rank==0)
    printf("its = %d, t = %15.10f, count = %d \n", its, t2, count);
}
  /************************************************************************************/
  PetscFinalize();
  return 0;
}
