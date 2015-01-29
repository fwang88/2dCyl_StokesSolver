#include <stdio.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <math.h>
#include "stokes.h"


#define PI 3.14159265358979323846
#define uIL(i) ((i)<0 ? (-1-(i)) : ((i)>=nr ? (2*nr-2-(i)) : (i)))
#define uJL(i) ((i)<0 ? (-(i)) : ((i)>=nz ? (2*nz-1-(i)) : (i)))
#define wIL(i) ((i)<0 ? (-(i)) : ((i)>=nr ? (2*nr-1-(i)) : (i)))
#define wJL(i) ((i)<0 ? (-1-(i)) : ((i)>=nz ? (2*nz-2-(i)) : (i)))
#define pIL(i) ((i)<0 ? (-(i)) : ((i)>=nr ? (2*nr-1-(i)) : (i)))
#define pJL(i) ((i)<0 ? (-(i)) : ((i)>=nz ? (2*nz-1-(i)) : (i)))

#define drinv (1.0/dr)
#define dzinv (1.0/dz)
#define drinv2 (1.0/(dr * dr))
#define dzinv2 (1.0 / (dz * dz))
#define drzinv (1.0 / (dr * dz))

/*#define rPI(i) ((i) < 0 ? -(i)*dr : ((i) >= nr ? ((2*nr-1)-(i))*dr :\
  (i)*dr))*/
#define rPI(i) ((i)*dr)

#define suIJ(i,j) ((i) < 0 || (i) >= nr ? -1 : +1)

/*#define ruI(i) ((i) < 0 ? (-0.5-(i))*dr : ((i) >= nr ? ((2*nr-1.5)-(\
  i))*dr : ((i)+0.5)*dr))*/
#define ruI(i) (((i)+0.5)*dr)

#define swIJ(i,j) ((j) < 0 || (j) >= nz ? -1 : +1)
#define rwI(i) rPI(i)


void getElementEqnums(MatStencil eqnu[], MatStencil eqnw[], MatStencil eqnp[],
		      PetscInt i, PetscInt j,
		      PetscInt nr, PetscInt nz) {
  
  eqnu[0].i = uIL(i)  ; eqnu[0].j = uJL(j)  ; eqnu[0].c = 0;
  eqnu[1].i = uIL(i+1); eqnu[1].j = uJL(j)  ; eqnu[1].c = 0;
  eqnu[2].i = uIL(i-1); eqnu[2].j = uJL(j)  ; eqnu[2].c = 0;
  eqnu[3].i = uIL(i)  ; eqnu[3].j = uJL(j+1); eqnu[3].c = 0;
  eqnu[4].i = uIL(i)  ; eqnu[4].j = uJL(j-1); eqnu[4].c = 0;

  eqnu[5].i = wIL(i)  ; eqnu[5].j = wJL(j)  ; eqnu[5].c = 1;
  eqnu[6].i = wIL(i+1); eqnu[6].j = wJL(j)  ; eqnu[6].c = 1;
  eqnu[7].i = wIL(i)  ; eqnu[7].j = wJL(j-1); eqnu[7].c = 1;
  eqnu[8].i = wIL(i+1); eqnu[8].j = wJL(j-1); eqnu[8].c = 1;

  eqnu[9].i = pIL(i)  ; eqnu[9].j = pJL(j)  ; eqnu[9].c = 2;
  eqnu[10].i = pIL(i+1);eqnu[10].j = pJL(j) ; eqnu[10].c = 2;

  /***/
  eqnw[0].i = uIL(i)  ; eqnw[0].j = uJL(j)  ; eqnw[0].c = 0;
  eqnw[1].i = uIL(i-1); eqnw[1].j = uJL(j)  ; eqnw[1].c = 0;
  eqnw[2].i = uIL(i)  ; eqnw[2].j = uJL(j+1); eqnw[2].c = 0;
  eqnw[3].i = uIL(i-1); eqnw[3].j = uJL(j+1); eqnw[3].c = 0;

  eqnw[4].i = wIL(i)  ; eqnw[4].j = wJL(j)  ; eqnw[4].c = 1;
  eqnw[5].i = wIL(i+1); eqnw[5].j = wJL(j)  ; eqnw[5].c = 1;
  eqnw[6].i = wIL(i-1); eqnw[6].j = wJL(j)  ; eqnw[6].c = 1;
  eqnw[7].i = wIL(i)  ; eqnw[7].j = wJL(j+1); eqnw[7].c = 1;
  eqnw[8].i = wIL(i)  ; eqnw[8].j = wJL(j-1); eqnw[8].c = 1;

  eqnw[9].i = pIL(i)  ; eqnw[9].j = pJL(j)  ; eqnw[9].c = 2;
  eqnw[10].i = pIL(i) ; eqnw[10].j = pJL(j+1);eqnw[10].c = 2;

  /***/
  eqnp[0].i = uIL(i)  ; eqnp[0].j = uJL(j)  ; eqnp[0].c = 0;
  eqnp[1].i = uIL(i-1); eqnp[1].j = uJL(j)  ; eqnp[1].c = 0;
  eqnp[2].i = wIL(i)  ; eqnp[2].j = wJL(j)  ; eqnp[2].c = 1;
  eqnp[3].i = wIL(i)  ; eqnp[3].j = wJL(j-1); eqnp[3].c = 1;

}

/** function heaviside function **/

static double heavi(double G, double epsilon, double* val1, double* val2, int j, int nz)
{ 
  if (j >= nz) j = nz-1;
  if (j < 0) j = 1;
  if (G <= -epsilon) return val1[j];
  if (G >= epsilon) return val2[j];
  return val1[j] + (val2[j]-val1[j])*(1+G/epsilon+sin(PI*G/epsilon)/PI)/2;
}

/** get matrix coefficients **/
void eqnval(DM da, 
	    PetscScalar *valu, PetscScalar *valw, PetscScalar *valp, 
	    PetscInt i, PetscInt j,
	    PetscScalar *mu1, PetscScalar *mu2,
	    Vec lG, 
	    PetscScalar dr, PetscScalar dz, 
	    PetscInt nr, PetscInt nz,
	    PetscScalar epsilon
	    ) {

  PetscScalar **array;
  DMDAVecGetArray(da,lG,&array);
  
  PetscInt rank, k;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  for (k=0; k<11; k++) {
    valu[k]=0; 
    valw[k]=0;
  }

  for (k=0; k<4; k++) {
    valp[k]=0;
  }
  
#define mu(gij,j) heavi(gij,epsilon,mu1,mu2,j,nz)
#define array(i,j) array[pJL(j)][pIL(i)] //i --- r; j --- z;                                       
  {
    double avgmuij = 0.25 * (mu(array(i,j),j) + mu(array(i+1,j),j) + mu(array(i,j+1),j+1) + mu(array(i+1,j+1),j+1));
    double avgmuijm1 = 0.25 * (mu(array(i,j),j) + mu(array(i+1,j),j) + mu(array(i,j-1),j-1) + mu(array(i+1,j-1),j-1));
    double avg2muij = 0.5 * (mu(array(i,j),j) + mu(array(i+1,j),j));

    double cu = 2 * rPI(i) * mu(array(i,j),j) * drinv2;
    double cup = 2 * rPI(i+1) * mu(array(i+1,j),j) * drinv2;
    double ruinv = 1.0 / ruI(i);
    double cu2 = 2.0 * avg2muij * ruinv * ruinv;
    double mup = avgmuij * dzinv2;
    double mum = avgmuijm1 * dzinv2;

    valu[0] = suIJ(i,j) * (-(cu+cup)*ruinv - (mup+mum) - cu2);

    valu[1] = suIJ(i+1,j) * cup * ruinv;
    valu[2] = suIJ(i-1,j) * cu * ruinv;
    valu[3] = suIJ(i,j+1) * mup;
    valu[4] = suIJ(i,j-1) * mum;

    mup = avgmuij * drzinv;
    mum = avgmuijm1 * drzinv;
    valu[5] = swIJ(i,j) * -mup;
    valu[6] = swIJ(i+1,j) * mup;
    valu[7] = swIJ(i,j-1) * mum;
    valu[8] = swIJ(i+1,j-1) * -mum;

    valu[9] = drinv; //i,j                                                                         
    valu[10] = -drinv; //i+1,j                                                 
  }

  {
    double avgmuij = 0.25 * (mu(array(i,j),j) + mu(array(i+1,j),j) + mu(array(i,j+1),j+1) + mu(array(i+1,j+1),j+1));
    double avgmuim1j = 0.25 * (mu(array(i,j),j) + mu(array(i-1,j),j) + mu(array(i,j+1),j+1) + mu(array(i-1,j+1),j+1));

    double mup = (i ? ruI(i)/rwI(i) : 2) * avgmuij * drzinv;
    double mum = (i ? ruI(i-1)/rwI(i) : 2) * avgmuim1j * drzinv;

    valw[0] = suIJ(i,j) * -mup;
    valw[1] = suIJ(i-1,j) * mum;
    valw[2] = suIJ(i,j+1) * mup;
    valw[3] = suIJ(i-1,j+1) * -mum;

    double cu = 2 * mu(array(i,j),j) * dzinv2;
    double cup = 2 * mu(array(i,j+1),j+1) * dzinv2;

    mup = (i ? ruI(i)/rwI(i) : 2) * (avgmuij*drinv2);
    mum = (i ? ruI(i-1)/rwI(i) : 2) * (avgmuim1j*drinv2);
    valw[4] = swIJ(i,j) * (-(cu+cup) - (mup+mum));
    valw[5] = swIJ(i+1,j) * mup;
    valw[6] = swIJ(i-1,j) * mum;
    valw[7] = swIJ(i,j+1) * cup;
    valw[8] = swIJ(i,j-1) * cu;

    valw[9] = +dzinv;
    valw[10] = -dzinv;
  }

  {
    double cu = (i ? ruI(i)/rPI(i) : 2) * drinv;
    double cm = (i ? ruI(i-1)/rPI(i) : 2) * drinv;
    valp[0] = suIJ(i,j) * cu;
    valp[1] = -suIJ(i-1,j) * cm;
    valp[2] = swIJ(i,j) * dzinv;
    valp[3] = -swIJ(i,j-1) * dzinv;
  }

  DMDAVecRestoreArray(da,lG,&array);
}


/**  function get_B_from_G **/

void get_B_from_G(stokes_matrix *B, levelset_vec *G, parameter *para, viscosity *mu) {
  
  PetscInt rank, bgnr, bgnz, szr, szz;
  DM da, da3;

  da = G->da;
  da3 = B->da;
  
  Vec lG;
  DMGetLocalVector(da, &lG);
  VecSet(lG, 0);
  DMGlobalToLocalBegin(da, G->data, INSERT_VALUES, lG);
  DMGlobalToLocalEnd  (da, G->data, INSERT_VALUES, lG);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  DMDAGetCorners(da3, &bgnr, &bgnz, 0, &szr, &szz, 0);
  
  MatZeroEntries(B->data);
  
  PetscInt i, j;
  MatStencil ueqn, weqn, peqn,
    eqnu[11], eqnw[11], eqnp[4];
  PetscScalar valu[11]={0}, valw[11]={0}, valp[4]={0};
  
  for(j = bgnz; j < bgnz + szz; j++) {
    for(i = bgnr; i < bgnr + szr; i++) {
      
      ueqn.i = i; ueqn.j = j; ueqn.c = 0;
      weqn.i = i; weqn.j = j; weqn.c = 1;
      peqn.i = i; peqn.j = j; peqn.c = 2;

      getElementEqnums(eqnu, eqnw, eqnp,
		       i, j, 
		       para->nr, para->nz
		       );
      
      eqnval(da, 
	     valu, valw, valp,
	     i, j,
	     mu->mu1, mu->mu2,
	     lG,
	     para->dr, para->dz,
	     para->nr, para->nz,
	     para->epsilon
	     );

      MatSetValuesStencil(B->data, 1, &ueqn, 11, eqnu, valu, ADD_VALUES);
      /*
      int kk;
      for(kk=0;kk<11;kk++)
	if(i==0 && j==0) printf("%g\n",valu[kk]);
      */
      MatSetValuesStencil(B->data, 1, &weqn, 11, eqnw, valw, ADD_VALUES);

      if(i==(para->nr-1) && j==(para->nz-1)) {
        MatSetValue(B->data, 
		    3*para->nr*para->nz-1, 3*para->nr*para->nz-1, 
		    1.0, ADD_VALUES);
      }
      else {
        MatSetValuesStencil(B->data, 1, &peqn, 4, eqnp, valp, ADD_VALUES);
      }

    }
  }
  
  MatAssemblyBegin(B->data, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B->data, MAT_FINAL_ASSEMBLY);

  DMRestoreLocalVector(da, &lG);

}



/** function get_force_from_G **/

#define IL(i) pIL(i)
#define JL(i) pJL(i)

#define sr(j,i) ((i) < 0 || (i) >= nr ? -1 : +1)
#define sz(j,i) ((j) < 0 || (j) >= nz ? -1 : +1)

#if !defined(MIN)
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#endif

static double delta(double G, double epsilon ) {
  if (G <= -epsilon) return 0;
  if (G >= epsilon) return 0;
  return  (1+ cos(PI*G/epsilon))/epsilon/2;
}


double minmag(const double a,const double b) {
  return fabs(a)<fabs(b) ? a:b;
}

double fsign(const double phi) {
  if(phi<=0.0) return -1.0;
  else return 1.0;
}

/** curvature function in 2d case **/

double curvature2D(double **garray, 
		   int i, int j, 
		   int nz, int nr, 
		   double dz, double dr) {
  
  double denominator, 
    llz, llr, llzz, llrr, llzr, 
    kappa, mindx;
  
  mindx = MIN(dz, dr);

  llz = ( garray[JL(j+1)][IL(i)] \
	  - garray[JL(j-1)][IL(i)] ) / (2.0*dz);
  llr = ( garray[JL(j)][IL(i+1)] \
	  - garray[JL(j)][IL(i-1)] ) / (2.0*dr);

  llzz = ( garray[JL(j+1)][IL(i)] \
	   - 2.0*garray[JL(j)][IL(i)] \
	   + garray[JL(j-1)][IL(i)] )  / (dz*dz);

  llrr = ( garray[JL(j)][IL(i+1)] \
	   - 2.0*garray[JL(j)][IL(i)] \
	   + garray[JL(j)][IL(i-1)] ) / (dr*dr);

  llzr = ( garray[JL(j+1)][IL(i+1)] \
	   - garray[JL(j-1)][IL(i+1)] \
	   - garray[JL(j+1)][IL(i-1)] \
	   + garray[JL(j-1)][IL(i-1)] ) / (4.0*dr*dz);

  denominator = sqrt( (llz*llz) + (llr*llr) );

  if(denominator >= 1.0e-14 && i>0) {
    kappa = ((llz*llz)*llrr - 2.0*llz*llr*llzr + (llr*llr)*llzz) \
      / (denominator*denominator*denominator) \
      + 1/(i*dr) * llr / denominator;
  }
  else if (denominator >= 1.0e-14 && i==0) {
    kappa = +((llz*llz)*llrr - 2.0*llz*llr*llzr + (llr*llr)*llzz) \
      / (denominator*denominator*denominator) \
      + llrr / denominator;
  }
  else { /* max resolution = 1/dx or -1/dx */
    kappa=fsign(garray[JL(j)][IL(i)])/mindx ;
  }

  kappa=minmag(kappa,fsign(kappa)/mindx);

  return kappa;
}


void get_force_from_G(stokes_force *force, levelset_vec *G, parameter *para)
{
  
  Vec l_G, l_force;
  PetscScalar **ptr_g;
  StokesDOF **ptr_f;

  PetscInt llr, llz, lsizer, lsizez,
    i, j,
    nr, nz;

  double gijz, gijr,
    length_gij, 
    deltaij, 
    curveij, 
    coeff_gij,
    dr, dz;

  DM da, da3;

  nr = para->nr;
  nz = para->nz;

  dr = para->dr;
  dz = para->dz;

  da = G->da;
  da3 = force->da;

  VecSet(force->data, 0);

  DMGetLocalVector(da, &l_G);
  DMGetLocalVector(da3, &l_force);

  VecSet(l_G,0);
  VecSet(l_force,0);

  DMGlobalToLocalBegin(da, G->data, INSERT_VALUES, l_G);
  DMGlobalToLocalEnd  (da, G->data, INSERT_VALUES, l_G);

  DMGlobalToLocalBegin(da3, force->data, INSERT_VALUES, l_force);
  DMGlobalToLocalEnd  (da3, force->data, INSERT_VALUES, l_force);

  DMDAVecGetArray(da, l_G, &ptr_g);
  DMDAVecGetArray(da3, l_force, &ptr_f);

  DMDAGetCorners(da3, &llr, &llz,0, &lsizer, &lsizez, 0);


  for(j=llz; j<llz+lsizez; j++) {
    for(i=llr; i<llr+lsizer; i++) {

      deltaij = delta(ptr_g[JL(j)][IL(i)], para->epsilon);

      curveij = curvature2D(ptr_g, 
			    i, j, 
			    para->nz, para->nr, 
			    para->dz, para->dr
			    );

      gijz = ( ptr_g[JL(j+1)][IL(i)] - ptr_g[JL(j-1)][IL(i)] ) / (2.0*dz);
      gijr = ( ptr_g[JL(j)][IL(i+1)] - ptr_g[JL(j)][IL(i-1)] ) / (2.0*dr);

      length_gij = sqrt((gijr*gijr) + (gijz*gijz));

      if(length_gij>1e-14) {
        coeff_gij = deltaij * curveij / length_gij;
      } 
      else {
        coeff_gij = 0;
      }

      ptr_f[j][i].u = para->tension * coeff_gij * gijr;
      ptr_f[j][i].w = para->tension * coeff_gij * gijz;
    }
  }

  DMLocalToGlobalBegin(da3, l_force, INSERT_VALUES, force->data);
  DMLocalToGlobalEnd  (da3, l_force, INSERT_VALUES, force->data);

  VecAssemblyBegin(force->data);
  VecAssemblyEnd  (force->data);

  DMGlobalToLocalBegin(da3, force->data, INSERT_VALUES, l_force);
  DMGlobalToLocalEnd  (da3, force->data, INSERT_VALUES, l_force);

  DMDAVecRestoreArray(da3, l_force, &ptr_f);
  DMDAVecGetArray    (da3, l_force, &ptr_f);

  for(j=llz; j<llz+lsizez; j++) {
    for(i=llr; i<llr+lsizer; i++) {
      
      ptr_f[j][i].u = 0.5 * (ptr_f[j][i].u + ptr_f[j][pIL(i+1)].u * sr(j, i+1));
      ptr_f[j][i].w = 0.5 * (ptr_f[j][i].w + ptr_f[pJL(j+1)][i].w * sz(j+1, i));
    }
  }

  DMLocalToGlobalBegin(da3, l_force, INSERT_VALUES, force->data);
  DMLocalToGlobalEnd  (da3, l_force, INSERT_VALUES, force->data);

  VecAssemblyBegin(force->data);
  VecAssemblyEnd  (force->data);

  DMDAVecRestoreArray(da3, l_force, &ptr_f);
  DMDAVecRestoreArray(da, l_G, &ptr_g);

  DMRestoreLocalVector(da, &l_G);
  DMRestoreLocalVector(da3, &l_force);

}


/** solver **/
void initialize_MYKSP(MYKSP *solver) {
  KSPCreate(PETSC_COMM_WORLD, &solver->ksp);
  solver->niter = 100;
}

void stokes_solver_direct(MYKSP *solver, 
			  stokes_matrix *B, stokes_force *force,
			  stokes_force *uwp) {
  PC pc;

  KSPSetOperators(solver->ksp, B->data, B->data, 
		  DIFFERENT_NONZERO_PATTERN);

  KSPSetDM(solver->ksp, B->da);

  KSPSetDMActive(solver->ksp, PETSC_FALSE);
  KSPSetFromOptions(solver->ksp);

  KSPGetPC(solver->ksp, &pc);
  PCSetType(pc, PCLU);
  PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS);

  KSPSetUp(solver->ksp);

  KSPSolve(solver->ksp, force->data, uwp->data);

  solver->niter = 0;
}

void stokes_solver_precond(MYKSP *solver,
			   stokes_matrix *B,
			   stokes_force *force,
			   stokes_force *uwp) {
  
  KSPSetOperators(solver->ksp,
		  B->data, B->data,
		  SAME_PRECONDITIONER);

  KSPSetUp(solver->ksp);

  KSPSolve(solver->ksp, force->data, uwp->data);

  KSPGetIterationNumber(solver->ksp, &solver->niter);
  
}


void lab_frame_shift_w(stokes_force *uwp, PetscScalar vf) {
  
  PetscInt i, j;
  PetscInt llz, lsizez, llr, lsizer;
  StokesDOF **ptr_uwp;
  
  DMDAVecGetArray(uwp->da, uwp->data, &ptr_uwp);
  
  DMDAGetCorners(uwp->da, 
		 &llr, &llz,0, 
		 &lsizer, &lsizez, 0);

  for(j=llz; j<llz+lsizez; j++) {
    for(i=llr; i<llr+lsizer; i++) {

      ptr_uwp[j][i].w -= vf;
    
    }
  }
 
  DMDAVecRestoreArray(uwp->da, uwp->data, &ptr_uwp);
  
}


void stokes_solve(stokes_force *uwp, 
		  MYKSP *solver, 
		  stokes_matrix *B,
		  stokes_force *force,
		  levelset_vec *G,
		  parameter *para,
		  viscosity *mu) {

  //  get_B_from_G(B, G, para, mu);

  //  get_force_from_G(force, G, para);

  if(solver->niter >= 11) {
    stokes_solver_direct(solver, B, force, uwp);
  }
  else {
    stokes_solver_precond(solver, B, force, uwp);
  }
  

}
