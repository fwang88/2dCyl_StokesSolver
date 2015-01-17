
typedef struct {
  PetscScalar u, w, p;
} StokesDOF;


typedef struct {
  Vec data;
  PetscScalar maxR, maxZ, dr, dz;
  PetscInt nr, nz;
  DM da;
} levelset_vec;


levelset_vec create_levelset(PetscScalar maxR, PetscScalar maxZ,
			    PetscScalar dr,   PetscScalar dz,
			    PetscInt nghostlayer);

void destroy_levelset(levelset_vec);

void initial_levelset(levelset_vec *G, PetscScalar r0, PetscScalar pertb, char *mode);

void get_input(int argc, char **args,
	       PetscScalar *maxR, PetscScalar *maxZ,
	       PetscScalar *dr,   PetscScalar *dz);

PETSC_EXTERN void InitialLevelSet(PetscScalar, PetscScalar, PetscInt, PetscInt, PetscScalar, PetscScalar, PetscScalar, PetscScalar **, PetscInt, PetscInt, PetscInt, PetscInt, char *);

PETSC_EXTERN void RK2DReinit(Vec, PetscInt, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_Reinit_2D(Vec, Vec, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_mex2D(DM, DM, Vec, Vec, Vec, PetscScalar,  PetscScalar, PetscInt, PetscInt, PetscScalar);

PETSC_EXTERN void RK2DPeriod(DM, DM, Vec, Vec, PetscScalar, PetscScalar, PetscInt, PetscInt, PetscScalar );

PETSC_EXTERN void AssembleB( Mat, DM, DM , PetscScalar , PetscInt , PetscInt , PetscScalar , PetscScalar , Vec , PetscScalar *, PetscScalar *);

void AssembleForce(Vec , Vec , DM , DM , PetscInt , PetscInt , PetscScalar , PetscScalar , PetscScalar, PetscScalar );

double timestep(DM, Vec, PetscScalar, PetscScalar, PetscScalar);
