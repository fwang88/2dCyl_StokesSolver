#include <petscksp.h>
#include <sys/stat.h>
#include <string.h>
#include "stokes.h"

void mk_dir(parameter *para){
  char buffer[1000];
  
  sprintf(buffer, "./data");
  mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth, para->lowtwidth);
  mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  return;
}


void output_levelset(levelset_vec *G, parameter *para) {
  
  PetscViewer viewer;
  char buffer[1000];

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth, para->lowtwidth, G->t);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(G->data, viewer);
  PetscViewerDestroy(&viewer);
  
  return;
}

void output_uwp(stokes_force *uwp, parameter *para){
  
  PetscViewer viewer;
  char buffer[1000];

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f/outputUWP_t_%.6f", para->maxz, para->maxr, para->dz, para->vf,para->tension, para->r0, para->tlow, para->thigh, para->twidth, para->lowtwidth, uwp->t);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(uwp->data, viewer);
  PetscViewerDestroy(&viewer);
  //  printf("output uwp \n");  
  return;
}

void output_force(stokes_force *force, parameter *para){
  
  PetscViewer viewer;
  char buffer[1000];

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f/outputforce_t_%.6f", para->maxz, para->maxr, para->dz, para->vf,para->tension, para->r0, para->tlow, para->thigh, para->twidth, para->lowtwidth, force->t);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(force->data, viewer);
  PetscViewerDestroy(&viewer);
  //  printf("output uwp \n");  
  return;
}


void load_levelset(levelset_vec *G, parameter *para) {
  char buffer[1000];
  PetscViewer viewer;

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth, para->lowtwidth, para->trestart);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_READ, &viewer);
  VecLoad(G->data, viewer);
  PetscViewerDestroy(&viewer);
}
