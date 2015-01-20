#include <petscksp.h>
#include <petscdmda.h>
#include <stdio.h>
#include "stokes.h"

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#if !defined(IDL)
#define IDL(i, j) (((i)<0 ? (-(i)) : ((i)>=nz) ? (2*nz-1-(i)) : (i))*nr + ((j)<0 ? (-(j)) : ((j)>=nr) ? (2*nr-1-(j)) : (j)) - llz*nr)
#endif

/*****************************************************************************/

double sqr(double a) {
  return pow(a,2);
}

int Sign(const double phi) {
  if(phi <= 0.0) return -1;
  else return 1;
}


void reinit(levelset_vec *G, parameter *para) {
  
  int i, rank, llr, llz, lsizer, lsizez;
  PetscScalar epsilon, dtau;
  Vec G1, G2;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  epsilon = 2.0 * MIN(para->dz, para->dr);
  dtau = 0.1 * MIN(para->dz,para->dr);
  
  DMDAGetCorners(G->da,&llr,&llz,0,&lsizer,&lsizez,0);
  
  VecDuplicate(G->data, &G1);
  VecDuplicate(G->data, &G2);
  
  for(i=0; i<para->reinitstep; i++) {
    reinit_kernal(G1, G->data,
		  dtau, epsilon, 
		  para->dz, para->dr, 
		  para->nz, para->nr, 
		  G->da);
    
    reinit_kernal(G2, G1, 
		  dtau, epsilon, 
		  para->dz, para->dr, 
		  para->nz, para->nr, 
		  G->da);
    
    VecCopy(G->data, G1);
    VecScale(G1, 0.75);
    VecAXPY(G1, 0.25, G2);
    
    reinit_kernal(G2, G1, 
		  dtau, epsilon, 
		  para->dz, para->dr, 
		  para->nz, para->nr, 
		  G->da);
    
    VecScale(G->data, 1.0/3.0);
    VecAXPY(G->data, 2.0/3.0, G2);
    
  }
  
  VecDestroy(&G1);
  VecDestroy(&G2);
  
  return;
}


void reinit_kernal(Vec G1, Vec G, 
		   PetscScalar dtau, PetscScalar epsilon,
		   PetscScalar dz, PetscScalar dr,
		   PetscInt nz, PetscInt nr,
		   DM da
		   ) {
  
  PetscInt i, j, llr, llz, lsizer, lsizez, rank;
  Vec localG;
  PetscScalar **plg, **plg1;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  DMCreateLocalVector(da, &localG);
  VecSet(localG, 0);
  DMGlobalToLocalBegin(da, G, INSERT_VALUES, localG);
  DMGlobalToLocalEnd(da, G, INSERT_VALUES, localG);
  DMDAVecGetArray(da, localG, &plg);
  DMDAVecGetArray(da, G1, &plg1);
  
  DMDAGetCorners(da,
		 &llr, &llz, 0,
		 &lsizer, &lsizez, 0);
  
  static double v1, v2, v3, v4, v5, vtmp;
  double der_z_minus, der_z_plus, 
    der_r_minus, der_r_plus,
    G_z_squared, G_r_squared,
    smoothing_factor;
  
  for (i=llz; i<llz+lsizez; i++)
    for (j=llr; j<llr+lsizer; j++) {

      v1 = (plg[llz][IDL(i-2,j)] - plg[llz][IDL(i-3,j)])/dz;
      v2 = (plg[llz][IDL(i-1,j)] - plg[llz][IDL(i-2,j)])/dz;
      v3 = (plg[llz][IDL(i,j)] - plg[llz][IDL(i-1,j)])/dz;
      v4 = (plg[llz][IDL(i+1,j)] - plg[llz][IDL(i,j)])/dz;
      v5 = (plg[llz][IDL(i+2,j)] - plg[llz][IDL(i+1,j)])/dz;
      der_z_minus = upwind_WENO_reinit(v1, v2, v3, v4, v5);
      
      v1 = (plg[llz][IDL(i+3,j)] - plg[llz][IDL(i+2,j)])/dz;
      vtmp = v2; v2 = v5; v5 = vtmp;
      vtmp = v3; v3 = v4; v4 = vtmp;
      der_z_plus = upwind_WENO_reinit(v1,v2,v3,v4,v5);
      
      v1 = (plg[llz][IDL(i,j-2)] - plg[llz][IDL(i,j-3)])/dr;
      v2 = (plg[llz][IDL(i,j-1)] - plg[llz][IDL(i,j-2)])/dr;
      v3 = (plg[llz][IDL(i,j)] - plg[llz][IDL(i,j-1)])/dr;
      v4 = (plg[llz][IDL(i,j+1)] - plg[llz][IDL(i,j)])/dr;
      v5 = (plg[llz][IDL(i,j+2)] - plg[llz][IDL(i,j+1)])/dr;
      der_r_minus = upwind_WENO_reinit(v1,v2,v3,v4,v5);
      
      v1 = (plg[llz][IDL(i,j+3)] - plg[llz][IDL(i,j+2)])/dr;
      vtmp = v2; v2 = v5; v5 = vtmp;
      vtmp = v3; v3 = v4; v4 = vtmp;
      der_r_plus = upwind_WENO_reinit(v1,v2,v3,v4,v5);
      
      if (Sign(plg[llz][IDL(i,j)]) > 0) {
	G_z_squared = MAX( sqr(MAX(der_z_minus,0)) , sqr( MIN(der_z_plus,0)) );
	G_r_squared = MAX( sqr(MAX(der_r_minus,0)) , sqr( MIN(der_r_plus,0)) );
      }
      else {
	G_z_squared = MAX( sqr(MIN(der_z_minus,0)) , sqr( MAX(der_z_plus,0)) );
	G_r_squared = MAX( sqr(MIN(der_r_minus,0)) , sqr( MAX(der_r_plus,0)) );
      }
      
      smoothing_factor = plg[llz][IDL(i,j)] / sqrt(sqr(plg[llz][IDL(i,j)]) + sqr(epsilon));
      plg1[i][j] = plg[llz][IDL(i,j)] - dtau * smoothing_factor * (sqrt(G_z_squared + G_r_squared)-1.0);
    }
  
  DMDAVecRestoreArray(da,localG,&plg);
  DMDAVecRestoreArray(da,G1,&plg1);
  
  VecDestroy(&localG);
  
  return;
}


double upwind_WENO_reinit(double v1,double v2,double v3,double v4,double v5)
{
  double data1 = v1/3.0 - 7.0*v2/6.0 + 11.0*v3/6.0;
  double data2 = -v2/6.0 + 5.0*v3/6.0 + v4/3.0;
  double data3 = v3/3.0 + 5.0*v4/6.0 - v5/6.0;

  double S1 = (13.0/12.0)*pow(v1-2.0*v2+v3, 2) + 0.25*pow(v1-4.0*v2+3.0*v3, 2);
  double S2 = (13.0/12.0)*pow(v2-2.0*v3+v4, 2) + 0.25*pow(v2-v4, 2);
  double S3 = (13.0/12.0)*pow(v3-2.0*v4+v5, 2) + 0.25*pow(3.0*v3-4.0*v4+v5, 2);

  double epsilon = 1e-6*(MAX(v1*v1,MAX(v2*v2,MAX(v3*v3,MAX(v4*v4,v5*v5))))) + 1e-99;

  double alpha1 = 0.1/(pow(S1+epsilon,2));
  double alpha2 = 0.6/(pow(S2+epsilon,2));
  double alpha3 = 0.3/(pow(S3+epsilon,2));
  double a_total = alpha1+alpha2+alpha3;

  return (alpha1/a_total)*data1 + (alpha2/a_total)*data2 + (alpha3/a_total)*data3;

}
