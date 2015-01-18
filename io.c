#include <petscksp.h>
#include <sys/stat.h>
#include <string.h>
#include "stokes.h"

void mk_dir(parameter *para){
  char buffer[1000];
  
  sprintf(buffer, "./data");
  int status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f", para->maxz, para->maxr, para->dz, para->vf, para->tlow, para->thigh, para->twidth, para->lowtwidth);
  status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  return;
}


void output(levelset_vec *G, PetscScalar time, parameter *para) {
  
  PetscViewer viewer;
  char buffer[1000];

  sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.2f_tlow_%.1f_thigh_%.1f_twidth_%.1f_lowtwidth_%.1f/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tlow, para->thigh, para->twidth, para->lowtwidth, time);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(G->data, viewer);
  PetscViewerDestroy(&viewer);
  
  return;
}
