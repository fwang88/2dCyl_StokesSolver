#include <petscksp.h>
#include <petscdmda.h>
#include <time.h>
#include <math.h>
#include "myheader.h"

#define IDL(i,j) (((i)<0 ? (-(i)) : ((i)>=nz) ? (2*nz-1-(i)) : (i))*nr + ((j)<0 ? (-(j)) : ((j)>=nr) ? (2*nr-1-(j)) : (j)) - llz*nr)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

double upwind_WENO(double v1,double v2,double v3,double v4,double v5);


/****************************************************************************************/
void Euler_mex2D(DM da, DM da3, Vec G1, Vec G, Vec uwp,  PetscScalar dz, PetscScalar dr, PetscInt nz, PetscInt nr, PetscScalar dt) {
	PetscInt i,j;
	Vec lG, luwp;
	PetscInt llr, llz, lsizer, lsizez, rank;
	PetscScalar **lgarray, **g1array;
	StokesDOF **vel;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
        DMGetLocalVector(da,&lG);
        DMGetLocalVector(da3,&luwp);
        VecSet(lG,0);
        VecSet(luwp,0);
        DMGlobalToLocalBegin(da,G,INSERT_VALUES,lG);
        DMGlobalToLocalEnd(da,G,INSERT_VALUES,lG);
        DMGlobalToLocalBegin(da3,uwp,INSERT_VALUES,luwp);
        DMGlobalToLocalEnd(da3,uwp,INSERT_VALUES,luwp);
    
        DMDAVecGetArray(da,lG,&lgarray);
        DMDAVecGetArray(da,G1,&g1array);
        DMDAVecGetArray(da3,luwp,&vel);
    
        DMDAGetCorners(da,&llr,&llz,0,&lsizer,&lsizez,0);
    
	static double v1,v2,v3,v4,v5;
	static double phi_conv_z, phi_conv_r;
	static double u_center, w_center;

        //DMDAVecGetArray(da3,uwp,&vel);
    
	for (i=llz;i<llz+lsizez;i++) 
          for (j=llr;j<llr+lsizer;j++) {       
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

            if(w_center>0.0) {
              v1 = (lgarray[llz][IDL(i-2,j)] - lgarray[llz][IDL(i-3,j)])/dz;
              v2 = (lgarray[llz][IDL(i-1,j)] - lgarray[llz][IDL(i-2,j)])/dz;
              v3 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i-1,j)])/dz;
              v4 = (lgarray[llz][IDL(i+1,j)] - lgarray[llz][IDL(i,j)])/dz;
              v5 = (lgarray[llz][IDL(i+2,j)] - lgarray[llz][IDL(i+1,j)])/dz;
              phi_conv_z = upwind_WENO(v1, v2, v3, v4, v5);
            }
            else {
              v1 = (lgarray[llz][IDL(i+3,j)] - lgarray[llz][IDL(i+2,j)])/dz;
              v2 = (lgarray[llz][IDL(i+2,j)] - lgarray[llz][IDL(i+1,j)])/dz;
              v3 = (lgarray[llz][IDL(i+1,j)] - lgarray[llz][IDL(i,j)])/dz;
              v4 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i-1,j)])/dz;;
              v5 = (lgarray[llz][IDL(i-1,j)] - lgarray[llz][IDL(i-2,j)])/dz;
              phi_conv_z = upwind_WENO(v1,v2,v3,v4,v5);
            }
                
            if(u_center>0.0) {
              v1 = (lgarray[llz][IDL(i,j-2)] - lgarray[llz][IDL(i,j-3)])/dr;
              v2 = (lgarray[llz][IDL(i,j-1)] - lgarray[llz][IDL(i,j-2)])/dr;
              v3 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i,j-1)])/dr;
              v4 = (lgarray[llz][IDL(i,j+1)] - lgarray[llz][IDL(i,j)])/dr;
              v5 = (lgarray[llz][IDL(i,j+2)] - lgarray[llz][IDL(i,j+1)])/dr;
              phi_conv_r = upwind_WENO(v1,v2,v3,v4,v5);   
            }
            else {
              v1 = (lgarray[llz][IDL(i,j+3)] - lgarray[llz][IDL(i,j+2)])/dr;
              v2 = (lgarray[llz][IDL(i,j+2)] - lgarray[llz][IDL(i,j+1)])/dr;
              v3 = (lgarray[llz][IDL(i,j+1)] - lgarray[llz][IDL(i,j)])/dr;
              v4 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i,j-1)])/dr;
              v5 = (lgarray[llz][IDL(i,j-1)] - lgarray[llz][IDL(i,j-2)])/dr;
              phi_conv_r = upwind_WENO(v1,v2,v3,v4,v5);
            }               

            g1array[llz][IDL(i,j)] = lgarray[llz][IDL(i,j)] - (dt) * ( w_center*phi_conv_z + u_center*phi_conv_r);
          }
        //printf("dt=%f\n",dt);
        DMDAVecRestoreArray(da,lG,&lgarray);     
        DMDAVecRestoreArray(da,G1,&g1array);
        DMDAVecRestoreArray(da3,luwp,&vel);
        DMRestoreLocalVector(da,&lG);
        DMRestoreLocalVector(da3,&luwp);
}

/****************************************************************************************/
double upwind_WENO(double v1,double v2,double v3,double v4,double v5) {
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
