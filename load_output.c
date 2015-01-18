#include <petscksp.h>
#include <sys/stat.h>
#include <string.h>

void output(levelset_vec *G, PetscScalar time) {

  PetscViewer viewer;
  char buffer[300];
  sprintf(buffer,"./data");
  int status = mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  sprintf(buffer, "./data/outputG_t_%.6f", time);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(G,viewer);
  PetscViewerDestroy(&viewer);
}
