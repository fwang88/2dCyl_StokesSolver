
typedef struct {
  PetscScalar u, w, p;
} StokesDOF;

typedef struct {
  Vec data;
  PetscInt nr, nz;
  DM da;
  PetscScalar t;
} levelset_vec;

typedef struct {
  PetscScalar maxr, maxz, dr, dz, r0;
  PetscScalar tension;
  PetscScalar pertb;
  PetscInt period_or_end;
  PetscScalar mui, muo;
  PetscScalar vf;
  PetscInt temp_profile;
  PetscScalar tlow, thigh, twidth, lowtwidth;
  PetscScalar restart, trestart;
  PetscScalar outputdt;
  char mode[10];
  PetscInt nghostlayer;
  PetscScalar epsilon;
  PetscInt reinitstep;
  PetscInt nr, nz;
} parameter;

typedef struct {
  PetscScalar *mu1;
  PetscScalar *mu2;
} viscosity;

/****************************************************************************/

levelset_vec create_levelset(parameter *para);

void destroy_levelset(levelset_vec);

void initial_levelset(levelset_vec *G, parameter *para);

void get_input(int argc, char **args, viscosity *mu, parameter *para);

void mk_dir(parameter *para);

void output(levelset_vec *G, PetscScalar time, parameter *para);

void load_levelset(levelset_vec *G, parameter *para);

void reinit(levelset_vec *G, parameter *para);

void reinit_kernal(Vec G1, Vec G,
                   PetscScalar dtau, PetscScalar epsilon,
                   PetscScalar dz, PetscScalar dr,
                   PetscInt nz, PetscInt nr,
                   DM da
                   );

  double upwind_WENO_reinit(double v1,double v2,double v3,double v4,double v5);



PETSC_EXTERN void RK2DReinit(Vec, PetscInt, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_Reinit_2D(Vec, Vec, PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscInt, PetscInt, DM);

PETSC_EXTERN void Euler_mex2D(DM, DM, Vec, Vec, Vec, PetscScalar,  PetscScalar, PetscInt, PetscInt, PetscScalar);

PETSC_EXTERN void RK2DPeriod(DM, DM, Vec, Vec, PetscScalar, PetscScalar, PetscInt, PetscInt, PetscScalar );

PETSC_EXTERN void AssembleB( Mat, DM, DM , PetscScalar , PetscInt , PetscInt , PetscScalar , PetscScalar , Vec , PetscScalar *, PetscScalar *);

void AssembleForce(Vec , Vec , DM , DM , PetscInt , PetscInt , PetscScalar , PetscScalar , PetscScalar, PetscScalar );

double timestep(DM, Vec, PetscScalar, PetscScalar, PetscScalar);
