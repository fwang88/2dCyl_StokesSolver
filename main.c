#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <math.h>
#include <time.h>
#include "myheader.h"
#include <stdio.h>
#include <petscviewerhdf5.h>
#include <string.h>

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
  int initflag=1;
  char* mode;
  int temp_profile = 0;

  mode = malloc(11 * sizeof(char));
  strncpy(mode,"end",10);
  mode[10] = '\0';

  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  r0 = 0.97;
  muRat = 0.91;
  tension = 1.5; 
  maxR = 10;
  maxZ = 100;
  dr = 0.1;
  dz = 0.1;
  pertb = 0.02;
  PetscScalar Tlow = 1400;
  PetscScalar Thigh = 1800;
  PetscScalar Twidth = 295;
  PetscScalar outputstep = 1;
  vf = 2.0;
  PetscScalar LowTWidth = 350;
  PetscScalar trestart = 0.0;
  PetscInt restart=0;
  initflag = 1; /* 1=p; 0=end */
  temp_profile = 0;
  int hdf5 = 1;

  struct option opts[] = {
    { "r0", 1, NULL, 1 },
    { "muRat", 1, NULL, 2 },
    { "tension", 1, NULL, 3 },
    { "maxR", 1, NULL, 4 },
    { "maxZ", 1, NULL, 5 },
    { "dr", 1, NULL, 6 },
    { "dz", 1, NULL, 7 },
    { "pertb", 1, NULL, 8 },
    { "Tlow", 1, NULL, 9 },
    { "Thigh", 1, NULL, 10 },
    { "Twidth", 1, NULL, 11},
    { "outputstep", 1, NULL, 12},
    { "vf", 1, NULL, 13},
    { "LowTWidth", 1, NULL, 14},
    { "trestart", 1, NULL, 15},
    { "restart", 1, NULL, 16},
    { "initflag", 1, NULL, 17},    
    { "temp_profile", 1, NULL, 18},
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
      Tlow = atof(optarg);
      break;
    case 10:
      Thigh = atof(optarg);
      break;
    case 11:
      Twidth = atof(optarg);
      break;
    case 12:
      outputstep = atof(optarg);
      break;
    case 13:
      vf = atof(optarg);
      break;
    case 14:
      LowTWidth = atof(optarg);
      break;
    case 15:
      trestart = atof(optarg);
      break;
    case 16:
      restart = atof(optarg);
      break;
    case 17:
      initflag = atof(optarg);
      break;
    case 18:
      temp_profile = atof(optarg);
      break;
    case 'h':
      hdf5 = 1;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "unrecognized option");
    }
  }

  switch (initflag) {
  case 0:
    mode = "end";
    break;
  case 1:
    mode = "periodic";
    break;
  }  

  printf("%s\n",mode);

  interfacei = r0/dr;
  nr = round(maxR/dr);
  nz = round(maxZ/dz);
  
  if (rank == 0) {
    printf("running nr= %d x nz= %d simulation (maxR = %g, maxZ = %g) with restart time = %d, \
    //vf = %g, outputstep = %g \n ", nr, nz, maxR, maxZ, restart, vf, outputstep);
  }
  
  DMDACreate2d(PETSC_COMM_WORLD, br, bz, stype, nr, nz, 1, PETSC_DECIDE, 1,nghostlayer,NULL,NULL,&da);
  DMCreateGlobalVector(da,&G);
  VecSet(G,0);
  DMDAGetCorners(da,&llr,&llz,0,&lsizer,&lsizez,0);
  PetscObjectSetName((PetscObject)G, "levelset");
  
  epsilon = 3.0 * dr;
  PetscScalar t1, t2, z_width, Tatz, mu1[nz], mu2[nz];
  PetscScalar dt;
  PetscInt count, itsteps, its;
  char buffer[200];
  count = 1;
  z_width = Twidth/dz;
  
  for(i=0;i<nz;i++) {
    if(temp_profile == 1) {
      Tatz = (Thigh-Tlow)*(tanh((nz+1.0*LowTWidth/dz-i)/z_width)) + Tlow;
      mu2[i] = pow(10, 26909.0/(Tatz+273)-7.2348);
      mu1[i] = pow(10, 819.0/(Tatz+273)-3.727);
    }
    else {
        mu2[i] = 1.0;
        mu1[i] = mu2[i] * muRat;
      }
  }
  /** create matrix B, Force, uwp and DMDA structure **/
  Mat B; 
  Vec Force, uwp;
  StokesDOF **p_uwp;
  PetscInt dof=3;
  nghostlayer = 1;
  DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, nr, nz, 1, PETSC_DECIDE, dof, nghostlayer,NULL,NULL,&da3);
  DMCreateGlobalVector(da3,&Force);
  VecSet(Force,0);
  VecDuplicate(Force,&uwp);
  PetscObjectSetName((PetscObject)uwp,"velocity");
  DMSetMatrixPreallocateOnly(da3, PETSC_TRUE);
  DMSetUp(da3);
  DMCreateMatrix(da3,MATAIJ,&B);
  /** set up some parameters for ksp solver **/
  itsteps = 0;
  its = 51;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  DMDAGetCorners(da3,&llr,&llz,0,&lsizer,&lsizez,0);

  /** create the output directories and output initial level set functions **/
  sprintf(buffer,"./master_branch_tanh_restart");
  int status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  sprintf(buffer,"./master_branch_tanh_restart/mass_g_maxZ_%.2f_maxR_%.2f_dz_%.2f_twidth_%.2f_vf_%.2f_Tlow_%.1f_Thigh_%.1f_lowtwidth_%.1f",maxZ,maxR,dz,Twidth,vf,Tlow,Thigh,LowTWidth);
  status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 

  if (restart) {
    sprintf(buffer,"./master_branch_tanh_restart/mass_g_maxZ_%.2f_maxR_%.2f_dz_%.2f_twidth_%.2f_vf_%.2f_Tlow_%.1f_Thigh_%.1f_lowtwidth_%.1f/outputG_t_%.6f.h5",maxZ,maxR,dz,Twidth,vf,Tlow,Thigh,LowTWidth,trestart);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,buffer,FILE_MODE_READ,&viewer);
    VecLoad(G,viewer);
    PetscViewerDestroy(&viewer);
    sprintf(buffer,"./master_branch_tanh_restart/mass_g_maxZ_%.2f_maxR_%.2f_dz_%.2f_twidth_%.2f_vf_%.2f_Tlow_%.1f_Thigh_%.1f_lowtwidth_%.1f/begin_outputG_t_%.6f.h5",maxZ,maxR,dz,Twidth,vf,Tlow,Thigh,LowTWidth,trestart);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
    VecView(G,viewer);
    PetscViewerDestroy(&viewer);
  }
  else {
    trestart = 0.0;
    PetscScalar **array;
    DMDAVecGetArray(da,G,&array);
    InitialLevelSet(dr, dz, nr, nz, interfacei, r0, pertb, array, llr,llz,lsizer,lsizez,mode);
    DMDAVecRestoreArray(da,G,&array);
    PetscInt ReinitStep = 2;
    //    MPI_Barrier(MPI_COMM_WORLD);
    RK2DReinit(G, ReinitStep, dz, dr, nz, nr, da);
    sprintf(buffer,"./master_branch_tanh_restart/mass_g_maxZ_%.2f_maxR_%.2f_dz_%.2f_twidth_%.2f_vf_%.2f_Tlow_%.1f_Thigh_%.1f_lowtwidth_%.1f/outputG_t_%.6f.h5",maxZ,maxR,dz,Twidth,vf,Tlow,Thigh,LowTWidth,0.0);

    PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
    VecView(G,viewer);
    PetscViewerDestroy(&viewer);
  }
  
  t1 = trestart; t2 = t1; 
  
  /***************** time evolution of stokes equation *******************************/
  while(t2<(maxZ/vf)*3+trestart) {
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    itsteps++;
    AssembleB(B,da,da3,epsilon,nr,nz,dr,dz,G,mu1,mu2);
    AssembleForce(Force,G,da,da3,nr,nz,dr,dz,tension,epsilon);
    /** sparse direct solver and preconditioning iterative solver to solve stokes equation **/
    if(its>=11) {
      KSPSetOperators(ksp,B,B,DIFFERENT_NONZERO_PATTERN);
      KSPSetDM(ksp,da3);
      KSPSetDMActive(ksp,PETSC_FALSE);
      KSPSetFromOptions(ksp);
      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
      //PCFactorSetMatSolverPackage(pc,MATSOLVERPASTIX);
      KSPSetUp(ksp);
      KSPSolve(ksp,Force,uwp);
      its = 0;
    }
    else {
      KSPSetOperators(ksp,B,B,SAME_PRECONDITIONER);
      KSPSetUp(ksp);
      KSPSolve(ksp,Force,uwp);
      KSPGetIterationNumber(ksp,&its);
    }
    /** include the feed speed; get the time step by CFL condition **/
    DMDAVecGetArray(da3,uwp,&p_uwp);
    for(j=llz;j<llz+lsizez;j++) {
      for(i=llr;i<llr+lsizer;i++) {
        p_uwp[j][i].w -= vf;
      }
    }
    DMDAVecRestoreArray(da3,uwp,&p_uwp);
   
    dt = timestep(da3,uwp,dr,dz,vf);    
    /** iterate one step; reinitialize 5 steps of the level set every 6 time steps; update time sequence **/
    RK2DPeriod(da,da3,G,uwp,dz,dr,nz,nr,dt);

    if( itsteps % 6 == 0 ) {
      MPI_Barrier(MPI_COMM_WORLD);
      RK2DReinit(G,5,dz,dr,nz,nr,da);
    }
    t2 = t1+dt;
    /** output level set field and velocity field data in binary format **/
    if (1) {
      if( t1<count*outputstep+trestart && t2>=count*outputstep+trestart) {
        sprintf(buffer,"./master_branch_tanh_restart/mass_g_maxZ_%.2f_maxR_%.2f_dz_%.2f_twidth_%.2f_vf_%.2f_Tlow_%.1f_Thigh_%.1f_lowtwidth_%.1f/outputG_t_%.6f.h5",maxZ,maxR,dz,Twidth,vf,Tlow,Thigh,LowTWidth,count*outputstep+trestart);

        if( rank == 0) printf("output successfully at t2 = %.4f", t2);
        MPI_Barrier(MPI_COMM_WORLD);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
        VecView(G,viewer);
        PetscViewerDestroy(&viewer);
        /*
          sprintf(buffer,"./maxZ_%.2f/twidth_%.2f_vf_%.2f/uwp_t_%.6f.h5",maxZ,Twidth,vf,count*outputstep);
          PetscViewerHDF5Open(PETSC_COMM_WORLD,buffer,FILE_MODE_WRITE,&viewer);
          VecView(uwp,viewer);
          PetscViewerDestroy(&viewer);
        */
        count++;
      }
    }
    end = MPI_Wtime();
    if(rank==0 && itsteps<=10)
      printf("t1 = %.4f, t2 = %.4f, solve time = %.4f \n", t1, t2, end-begin);
    
    t1=t2;
}
  free(mode);
  PetscFinalize();
  return 0;
}
