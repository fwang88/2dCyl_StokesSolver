#ifndef HEADER_STOKES
#define HEADER_STOKES

typedef struct {
  Vec data;
  PetscInt nr, nz;
  DM da;
  PetscScalar t;
} levelset_vec;

typedef struct {
  Vec data;
  PetscInt nr, nz;
  DM da;
  PetscScalar t;
} stokes_force;

typedef struct {
  Mat data;
  PetscInt nr, nz;
  DM da;
} stokes_matrix;


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

typedef struct {
  PetscScalar u, w, p;
} StokesDOF;

typedef struct {
  KSP ksp;
  PetscInt niter;
} MYKSP;

/****************************************************************************/

levelset_vec create_levelset(parameter *para);

void destroy_levelset(levelset_vec);

stokes_matrix create_stokes_matrix(parameter *para);
//void create_stokes_matrix(stokes_matrix **B, parameter *para);

void assembly_stokes_matrix(stokes_matrix *B);

void initial_levelset(levelset_vec *G, parameter *para);

void get_input(int argc, char **args, viscosity *mu, parameter *para);

void mk_dir(parameter *para);

void output_levelset(levelset_vec *G, parameter *para);

void load_levelset(levelset_vec *G, parameter *para);

void output_uwp(stokes_force *uwp, parameter *para);

void output_force(stokes_force *force, parameter *para);

void reinit(levelset_vec *G, parameter *para);

void reinit_kernal(Vec G1, Vec G,
                   PetscScalar dtau, PetscScalar epsilon,
                   PetscScalar dz, PetscScalar dr,
                   PetscInt nz, PetscInt nr,
                   DM da
                   );

double upwind_WENO_reinit(double v1,double v2,double v3,double v4,double v5);


stokes_force create_stokes_force(parameter *para, DM da);

void get_B_from_G(stokes_matrix *B, levelset_vec *G, parameter *para, viscosity *mu);

void get_force_from_G(stokes_force *force, levelset_vec *G, parameter *para);

void initialize_MYKSP(MYKSP *solver);

void stokes_solver_direct(MYKSP *solver, 
			  stokes_matrix *B, stokes_force *force, 
			  stokes_force *uwp);

void stokes_solver_precond(MYKSP *solver,
			   stokes_matrix *B, stokes_force *force,
			   stokes_force *uwp);

void lab_frame_shift_w(stokes_force *uwp, PetscScalar vf);

PetscScalar advection_timestep(stokes_force *uwp, 
			       parameter *para);

void advection_evolve(levelset_vec *G, stokes_force *uwp,
		      parameter *para, PetscScalar dt);

void advection_evolve_kernal(DM da, DM da3,
                             Vec G1, Vec G,
                             Vec uwp,
                             PetscScalar dz, PetscScalar dr,
                             PetscInt nz, PetscInt nr,
                             PetscScalar dt);


void stokes_solve(stokes_force *uwp,
                  MYKSP *solver,
                  stokes_matrix *B,
                  stokes_force *force,
                  levelset_vec *G,
                  parameter *para,
                  viscosity *mu);

#endif
