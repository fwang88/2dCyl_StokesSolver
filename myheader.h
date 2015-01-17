
typedef struct {
  PetscScalar u, w, p;
} StokesDOF;


typedef struct {
  Vec data;
  PetscInt nr, nz;
  DM da;
} stokes_vec;


stokes_vec create_stokesvec(PetscInt, PetscInt, DM);


PETSC_EXTERN void InitialLevelSet(PetscScalar, PetscScalar, PetscInt, PetscInt, PetscScalar, PetscScalar, PetscScalar, PetscScalar **, PetscInt, PetscInt, PetscInt, PetscInt, char *);

PETSC_EXTERN void RK2DReinit(Vec, PetscInt, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_Reinit_2D(Vec, Vec, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_mex2D(DM, DM, Vec, Vec, Vec, PetscScalar,  PetscScalar, PetscInt, PetscInt, PetscScalar);

PETSC_EXTERN void RK2DPeriod(DM, DM, Vec, Vec, PetscScalar, PetscScalar, PetscInt, PetscInt, PetscScalar );

PETSC_EXTERN void AssembleB( Mat, DM, DM , PetscScalar , PetscInt , PetscInt , PetscScalar , PetscScalar , Vec , PetscScalar *, PetscScalar *);

void AssembleForce(Vec , Vec , DM , DM , PetscInt , PetscInt , PetscScalar , PetscScalar , PetscScalar, PetscScalar );

double timestep(DM, Vec, PetscScalar, PetscScalar, PetscScalar);
