#include <petscksp.h>
#include <petscdmda.h>
#include <time.h>
#include <math.h>

#define IDL(i,j) (((i)<0 ? (-(i)) : ((i)>=nz) ? (2*nz-1-(i)) : (i))*nr + ((j)<0 ? (-(j)) : ((j)>=nr) ? (2*nr-1-(j)) : (j)) - llz*nr)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

double upwind_WENO_Reinit(double v1,double v2,double v3,double v4,double v5);

double sqr(double a) {
	return pow(a,2);
}

int Sign(const double phi) {
	if(phi <= 0.0) return -1;
	else return 1;
}



void Euler_Reinit_2D(Vec G1, Vec G, PetscScalar dtau, PetscScalar epsilon, PetscScalar dz, PetscScalar dr, PetscInt nz, PetscInt nr, DM da) /* dz=dx dr=dy nz nr */
{
	PetscInt i, j;
	Vec lG;
	PetscInt llr, llz, lsizer, lsizez, rank;
	PetscScalar **lgarray, **g1array;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
        DMCreateLocalVector(da,&lG);
	VecSet(lG,0);
	DMGlobalToLocalBegin(da,G,INSERT_VALUES,lG);
        DMGlobalToLocalEnd(da,G,INSERT_VALUES,lG);
	DMDAVecGetArray(da,lG,&lgarray);
	DMDAVecGetArray(da,G1,&g1array);
	
	DMDAGetCorners(da,&llr,&llz,0,&lsizer,&lsizez,0);
	
	static double v1, v2, v3, v4, v5, vtmp;
	double der_z_minus, der_z_plus, der_r_minus, der_r_plus;
	double G_z_squared, G_r_squared; 
        double smoothing_factor;
	
        for (i=llz; i<llz+lsizez; i++) 
          for (j=llr; j<llr+lsizer; j++) {   
            v1 = (lgarray[llz][IDL(i-2,j)] - lgarray[llz][IDL(i-3,j)])/dz;
            v2 = (lgarray[llz][IDL(i-1,j)] - lgarray[llz][IDL(i-2,j)])/dz;
            v3 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i-1,j)])/dz;
            v4 = (lgarray[llz][IDL(i+1,j)] - lgarray[llz][IDL(i,j)])/dz;
            v5 = (lgarray[llz][IDL(i+2,j)] - lgarray[llz][IDL(i+1,j)])/dz;
            der_z_minus = upwind_WENO_Reinit(v1, v2, v3, v4, v5);
            
            v1 = (lgarray[llz][IDL(i+3,j)] - lgarray[llz][IDL(i+2,j)])/dz;
            vtmp = v2; v2 = v5; v5 = vtmp;
            vtmp = v3; v3 = v4; v4 = vtmp;
            der_z_plus = upwind_WENO_Reinit(v1,v2,v3,v4,v5);
            
            v1 = (lgarray[llz][IDL(i,j-2)] - lgarray[llz][IDL(i,j-3)])/dr;
            v2 = (lgarray[llz][IDL(i,j-1)] - lgarray[llz][IDL(i,j-2)])/dr;
            v3 = (lgarray[llz][IDL(i,j)] - lgarray[llz][IDL(i,j-1)])/dr;
            v4 = (lgarray[llz][IDL(i,j+1)] - lgarray[llz][IDL(i,j)])/dr;
            v5 = (lgarray[llz][IDL(i,j+2)] - lgarray[llz][IDL(i,j+1)])/dr;
            der_r_minus = upwind_WENO_Reinit(v1,v2,v3,v4,v5);   				    
            
            v1 = (lgarray[llz][IDL(i,j+3)] - lgarray[llz][IDL(i,j+2)])/dr;
            vtmp = v2; v2 = v5; v5 = vtmp;
            vtmp = v3; v3 = v4; v4 = vtmp;
            der_r_plus = upwind_WENO_Reinit(v1,v2,v3,v4,v5);
            
            if (Sign(lgarray[llz][IDL(i,j)]) > 0) {
              G_z_squared = MAX( sqr(MAX(der_z_minus,0)) , sqr( MIN(der_z_plus,0)) );
              G_r_squared = MAX( sqr(MAX(der_r_minus,0)) , sqr( MIN(der_r_plus,0)) ); 
            }
            else { 
              G_z_squared = MAX( sqr(MIN(der_z_minus,0)) , sqr( MAX(der_z_plus,0)) );
              G_r_squared = MAX( sqr(MIN(der_r_minus,0)) , sqr( MAX(der_r_plus,0)) );        
            }
            
            smoothing_factor = lgarray[llz][IDL(i,j)]/sqrt(sqr(lgarray[llz][IDL(i,j)])+sqr(epsilon));
            g1array[i][j] = lgarray[llz][IDL(i,j)] - (dtau) * smoothing_factor * (sqrt(G_z_squared + G_r_squared)-1.0);
          }
     	
        DMDAVecRestoreArray(da,lG,&lgarray);
     	DMDAVecRestoreArray(da,G1,&g1array);
	
        VecDestroy(&lG);
}



double upwind_WENO_Reinit(double v1,double v2,double v3,double v4,double v5)
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
