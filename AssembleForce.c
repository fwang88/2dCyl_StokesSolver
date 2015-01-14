#include <petscdmda.h>
#include <math.h>
#include "myheader.h"

#define PI 3.14159265358979323846

#define pIL(i) ((i)<0 ? (-(i)) : ((i)>=nr ? (2*nr-1-(i)) : (i)))
#define pJL(i) ((i)<0 ? (-(i)) : ((i)>=nz ? (2*nz-1-(i)) : (i)))
#define IL(i) pIL(i)
#define JL(i) pJL(i)

#define sr(j,i) ((i) < 0 || (i) >= nr ? -1 : +1)
#define sz(j,i) ((j) < 0 || (j) >= nz ? -1 : +1)

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
/****************************************************************************************/
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

/****************************************************************************************/
double curvature2D(double **garray, int i, int j, int nz, int nr, double dz, double dr) {
       /* x--z; y--r */
  double denominator,llz,llr, llzz,llrr, llzr, kappa, mindx;
  mindx = MIN(dz, dr);
  
  llz = (garray[JL(j+1)][IL(i)] - garray[JL(j-1)][IL(i)])/(2.0*dz);
  llr = (garray[JL(j)][IL(i+1)] - garray[JL(j)][IL(i-1)])/(2.0*dr);
  
  llzz = (garray[JL(j+1)][IL(i)]-2.0*garray[JL(j)][IL(i)]+garray[JL(j-1)][IL(i)])/(dz*dz);
  llrr = (garray[JL(j)][IL(i+1)]-2.0*garray[JL(j)][IL(i)]+garray[JL(j)][IL(i-1)])/(dr*dr);
  llzr = (garray[JL(j+1)][IL(i+1)] - garray[JL(j-1)][IL(i+1)] - garray[JL(j+1)][IL(i-1)] + garray[JL(j-1)][IL(i-1)])/(4.0*dr*dz);
  
  denominator = sqrt( (llz*llz) + (llr*llr) );
            
  if(denominator >= 1.0e-14 && i>0) {
    kappa = +((llz*llz)*llrr - 2.0*llz*llr*llzr + (llr*llr)*llzz) / (denominator*denominator*denominator) + 1/(i*dr) * llr / denominator;
  }
  else if (denominator >= 1.0e-14 && i==0) {
    kappa = +((llz*llz)*llrr - 2.0*llz*llr*llzr + (llr*llr)*llzz) / (denominator*denominator*denominator) + llrr / denominator;
  }
  else { /* max resolution = 1/dx or -1/dx */
    kappa=fsign(garray[JL(j)][IL(i)])/mindx ;
  }
  
  kappa=minmag(kappa,fsign(kappa)/mindx);
  
  return kappa;
}
/****************************************************************************************/

void AssembleForce(Vec Force, Vec G, DM da, DM da3, PetscInt nr, PetscInt nz, PetscScalar dr, PetscScalar dz, PetscScalar tension, PetscScalar epsilon) {

  Vec l_G, l_Force, l_Force2;
  PetscScalar **garray;
  StokesDOF **farray, **farray2;
  PetscInt llr,llz,lsizer,lsizez, i, j;
  double gijz, gijr, length_gij, deltaij, curveij, coeff_gij;
  VecSet(Force,0);
  DMGetLocalVector(da,&l_G);
  DMGetLocalVector(da3,&l_Force);
  DMGetLocalVector(da3,&l_Force2);
  VecSet(l_G,0); 
  VecSet(l_Force,0);
  DMGlobalToLocalBegin(da,G,INSERT_VALUES,l_G);
  DMGlobalToLocalEnd(da,G,INSERT_VALUES,l_G);
  DMGlobalToLocalBegin(da3,Force,INSERT_VALUES,l_Force);
  DMGlobalToLocalEnd(da3,Force,INSERT_VALUES,l_Force);
  //  VecCopy(l_Force,l_Force2);
  DMDAVecGetArray(da,l_G,&garray);
  DMDAVecGetArray(da3,l_Force,&farray);
  //  DMDAVecGetArray(da3,l_Force2,&farray2);
  DMDAGetCorners(da3,&llr,&llz,0,&lsizer,&lsizez,0);

  for(j=llz;j<llz+lsizez;j++) {
    for(i=llr;i<llr+lsizer;i++) {
      deltaij = delta(garray[JL(j)][IL(i)],epsilon);
      curveij = curvature2D(garray,i,j,nz,nr,dz,dr);
      gijz = (garray[JL(j+1)][IL(i)]-garray[JL(j-1)][IL(i)])/(2.0*dz);
      gijr = (garray[JL(j)][IL(i+1)]-garray[JL(j)][IL(i-1)])/(2.0*dr);
      length_gij = sqrt((gijr*gijr) + (gijz*gijz));
      if(length_gij>1e-14) {
        coeff_gij = deltaij * curveij / length_gij;
      } else {
        coeff_gij = 0;
      }
      farray[j][i].u = tension * coeff_gij * gijr;
      farray[j][i].w = tension * coeff_gij * gijz;
    }
  }
    
  DMLocalToGlobalBegin(da3,l_Force,INSERT_VALUES,Force);
  DMLocalToGlobalEnd(da3,l_Force,INSERT_VALUES,Force);
  VecAssemblyBegin(Force);
  VecAssemblyEnd(Force);
  DMGlobalToLocalBegin(da3,Force,INSERT_VALUES,l_Force);
  DMGlobalToLocalEnd(da3,Force,INSERT_VALUES,l_Force);
  
  DMDAVecGetArray(da3,l_Force,&farray);

  for(j=llz;j<llz+lsizez;j++) {
    for(i=llr;i<llr+lsizer;i++) {
      farray[j][i].u = 0.5 * (farray[j][i].u + farray[j][pIL(i+1)].u * sr(j,i+1));
      farray[j][i].w = 0.5 * (farray[j][i].w + farray[pJL(j+1)][i].w * sz(j+1,i));
    }
  }
  DMLocalToGlobalBegin(da3,l_Force,INSERT_VALUES,Force);
  DMLocalToGlobalEnd(da3,l_Force,INSERT_VALUES,Force);
  VecAssemblyBegin(Force);
  VecAssemblyEnd(Force);
  DMDAVecRestoreArray(da3,l_Force,&farray);
  //  DMDAVecRestoreArray(da3,l_Force2,&farray2);
  DMDAVecRestoreArray(da,l_G,&garray);
  DMRestoreLocalVector(da,&l_G);
  DMRestoreLocalVector(da3,&l_Force);
  //  DMRestoreLocalVector(da3,&l_Force2);

}
