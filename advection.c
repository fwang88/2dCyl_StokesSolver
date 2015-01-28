#include <stdio.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <time.h>
#include <math.h>
#include "stokes.h"


#define IDL(i,j) (((i)<0 ? (-(i)) : ((i)>=nz) ? (2*nz-1-(i)) : (i))\
		  *nr + ((j)<0 ? (-(j)) : ((j)>=nr) ? (2*nr-1-(j)) : (j)) - llz*nr)

#if !defined(MAX)
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#endif

double upwind_WENO(double v1,double v2,double v3,double v4,double v5);


PetscScalar advection_timestep(stokes_force *uwp, 
			parameter *para){

  if (para->period_or_end == 1) {
    return 0.01;
  }
  
  StokesDOF **ptr_uwp;
  DMDAVecGetArray(uwp->da, uwp->data, &ptr_uwp);

  PetscInt llr, llz, lsizer, lsizez, i, j;
  DMDAGetCorners(uwp->da,
		 &llr, &llz, 0, 
		 &lsizer, &lsizez, 0);

  PetscScalar maxu_proc = 0, maxw_proc = 0;

  for (i=llz; i<llz+lsizez; i++) {
      for (j=llr; j<llr+lsizer; j++) {

	if( fabs(ptr_uwp[i][j].u) > maxu_proc ) 
	  maxu_proc = fabs(ptr_uwp[i][j].u);

	if( fabs(ptr_uwp[i][j].w) > maxw_proc ) 
	  maxw_proc = fabs(ptr_uwp[i][j].w);

      }
  }
  
  PetscScalar dtu = para->dr/maxu_proc;
  PetscScalar dtw = para->dz/maxw_proc;
  PetscScalar dtmin = MIN(dtu, dtw);

  double dtresult;
  (void) MPI_Allreduce(&dtmin, &dtresult, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
  
  dtresult = 1.1 * dtresult;
  return (dtresult);
}


void advection_evolve(levelset_vec *G, stokes_force *uwp, 
		      parameter *para, PetscScalar dt) {
  
  Vec G1, G2;
  VecDuplicate(G->data, &G1);
  VecDuplicate(G->data, &G2);
  
  DM da, da3;
  da = G->da;
  da3 = uwp->da;

  advection_evolve_kernal(da, da3,
			  G1, G->data,
			  uwp->data,
			  para->dz, para->dr, 
			  para->nz, para->nr, 
			  dt);

  advection_evolve_kernal(da, da3,
			  G2, G1, 
			  uwp->data,
			  para->dz, para->dr,
			  para->nz, para->nr,
			  dt);

  VecCopy(G->data, G1);
  VecScale(G1, 0.75);
  VecAXPY(G1, 0.25, G2);
  VecSet(G2, 0);
  
  advection_evolve_kernal(da, da3,
			  G2, G1, 
			  uwp->data, 
			  para->dz, para->dr, 
			  para->nz, para->nr, 
			  dt);
  
  VecScale(G->data, 1.0/3.0);
  VecAXPY(G->data, 2.0/3.0, G2);
  
  VecDestroy(&G1);
  VecDestroy(&G2);

  return;
}


void advection_evolve_kernal(DM da, DM da3,
			     Vec G1, Vec G,
			     Vec uwp, 
			     PetscScalar dz, PetscScalar dr, 
			     PetscInt nz, PetscInt nr, 
			     PetscScalar dt) {

  PetscInt llr, llz, lsizer, lsizez, rank, i, j;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  Vec localG, localuwp;

  DMGetLocalVector(da, &localG);
  DMGetLocalVector(da3, &localuwp);

  VecSet(localG, 0);
  VecSet(localuwp, 0);

  DMGlobalToLocalBegin(da, G, INSERT_VALUES, localG);
  DMGlobalToLocalEnd  (da, G, INSERT_VALUES, localG);

  DMGlobalToLocalBegin(da3, uwp, INSERT_VALUES, localuwp);
  DMGlobalToLocalEnd  (da3, uwp, INSERT_VALUES, localuwp);

  PetscScalar **ptr_localG, **ptr_G1;
  StokesDOF **vel;

  DMDAVecGetArray(da, localG, &ptr_localG);
  DMDAVecGetArray(da, G1, &ptr_G1);
  DMDAVecGetArray(da3, localuwp, &vel);

  DMDAGetCorners(da, 
		 &llr, &llz, 0,
		 &lsizer, &lsizez, 0);

  static double v1,v2,v3,v4,v5;
  static double phi_conv_z, phi_conv_r;
  static double u_center, w_center;

  for (i=llz; i<llz+lsizez; i++)
    for (j=llr; j<llr+lsizer; j++) {

      if (i==0) {
	/*w_center = vel[llz][IDL(i,j)].w; */
	w_center = 0;
      }
      else {
	w_center = 0.5*(vel[llz][IDL(i,j)].w + vel[llz][IDL(i-1,j)].w);
      }
      
      if (j==0) {
	/*u_center = vel[llz][IDL(i,j)].u;*/
	u_center = 0;
      }
      else {
	u_center = 0.5*(vel[llz][IDL(i,j)].u + vel[llz][IDL(i,j-1)].u);
      }

      if(w_center > 0.0) {
	v1 = (ptr_localG[llz][IDL(i-2,j)] - ptr_localG[llz][IDL(i-3,j)])/dz;
	v2 = (ptr_localG[llz][IDL(i-1,j)] - ptr_localG[llz][IDL(i-2,j)])/dz;
	v3 = (ptr_localG[llz][IDL(i,j)]   - ptr_localG[llz][IDL(i-1,j)])/dz;
	v4 = (ptr_localG[llz][IDL(i+1,j)] - ptr_localG[llz][IDL(i,j)])  /dz;
	v5 = (ptr_localG[llz][IDL(i+2,j)] - ptr_localG[llz][IDL(i+1,j)])/dz;
	phi_conv_z = upwind_WENO(v1, v2, v3, v4, v5);
      }
      else {
	v1 = (ptr_localG[llz][IDL(i+3,j)] - ptr_localG[llz][IDL(i+2,j)])/dz;
	v2 = (ptr_localG[llz][IDL(i+2,j)] - ptr_localG[llz][IDL(i+1,j)])/dz;
	v3 = (ptr_localG[llz][IDL(i+1,j)] - ptr_localG[llz][IDL(i,j)])  /dz;
	v4 = (ptr_localG[llz][IDL(i,j)]   - ptr_localG[llz][IDL(i-1,j)])/dz;;
	v5 = (ptr_localG[llz][IDL(i-1,j)] - ptr_localG[llz][IDL(i-2,j)])/dz;
	phi_conv_z = upwind_WENO(v1, v2, v3, v4, v5);
      }

      if(u_center>0.0) {
	v1 = (ptr_localG[llz][IDL(i,j-2)] - ptr_localG[llz][IDL(i,j-3)])/dr;
	v2 = (ptr_localG[llz][IDL(i,j-1)] - ptr_localG[llz][IDL(i,j-2)])/dr;
	v3 = (ptr_localG[llz][IDL(i,j)]   - ptr_localG[llz][IDL(i,j-1)])/dr;
	v4 = (ptr_localG[llz][IDL(i,j+1)] - ptr_localG[llz][IDL(i,j)])  /dr;
	v5 = (ptr_localG[llz][IDL(i,j+2)] - ptr_localG[llz][IDL(i,j+1)])/dr;
	phi_conv_r = upwind_WENO(v1, v2, v3, v4, v5);
      }
      else {
	v1 = (ptr_localG[llz][IDL(i,j+3)] - ptr_localG[llz][IDL(i,j+2)])/dr;
	v2 = (ptr_localG[llz][IDL(i,j+2)] - ptr_localG[llz][IDL(i,j+1)])/dr;
	v3 = (ptr_localG[llz][IDL(i,j+1)] - ptr_localG[llz][IDL(i,j)])  /dr;
	v4 = (ptr_localG[llz][IDL(i,j)]   - ptr_localG[llz][IDL(i,j-1)])/dr;
	v5 = (ptr_localG[llz][IDL(i,j-1)] - ptr_localG[llz][IDL(i,j-2)])/dr;
	phi_conv_r = upwind_WENO(v1, v2, v3, v4, v5);
      }

      ptr_G1[llz][IDL(i,j)] = ptr_localG[llz][IDL(i,j)] - (dt) \
	* ( w_center*phi_conv_z + u_center*phi_conv_r);
    }
  
  DMDAVecRestoreArray(da, localG, &ptr_localG);
  DMDAVecRestoreArray(da, G1, &ptr_G1);
  DMDAVecRestoreArray(da3, localuwp, &vel);

  DMRestoreLocalVector(da, &localG);
  DMRestoreLocalVector(da3, &localuwp);

}


double upwind_WENO(double v1,double v2,double v3,double v4,double v5) 
{
  double data1 = v1/3.0 - 7.0*v2/6.0 + 11.0*v3/6.0;
  double data2 = -v2/6.0 + 5.0*v3/6.0 + v4/3.0;
  double data3 = v3/3.0 + 5.0*v4/6.0 - v5/6.0;

  double S1 = (13.0/12.0)*pow(v1-2.0*v2+v3, 2) + 0.25*pow(v1-4.0*\
							  v2+3.0*v3, 2);
  double S2 = (13.0/12.0)*pow(v2-2.0*v3+v4, 2) + 0.25*pow(v2-v4, \
							  2);
  double S3 = (13.0/12.0)*pow(v3-2.0*v4+v5, 2) + 0.25*pow(3.0*v3-\
							  4.0*v4+v5, 2);

  double epsilon = 1e-6*(MAX(v1*v1,MAX(v2*v2,MAX(v3*v3,MAX(v4*v4,\
							   v5*v5))))) + 1e-99;

  double alpha1 = 0.1/(pow(S1+epsilon,2));
  double alpha2 = 0.6/(pow(S2+epsilon,2));
  double alpha3 = 0.3/(pow(S3+epsilon,2));
  double a_total = alpha1+alpha2+alpha3;

  return (alpha1/a_total)*data1 + (alpha2/a_total)*data2 + (alpha3/a_total)*data3;
}


