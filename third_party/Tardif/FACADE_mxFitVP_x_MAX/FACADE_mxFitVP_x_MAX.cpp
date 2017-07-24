#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#define MATLAB

#ifdef MATLAB
#include<mex.h>
#endif 

#include <math.h>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <stdlib.h>



#include "../CMINPACK/f2c.h"

//#include "lapack.h"
#ifdef __cplusplus
extern "C" 
{
#endif
	/* Subroutine */  int DSYEV(char *jobz, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);
//#include "clapack.h"
#ifdef __cplusplus
}
#endif


#include "../CMINPACK/minpack.h"


//#define max(a,b) ((a) >= (b) ? (a) : (b))
//#define min(a,b) ((a) <= (b) ? (a) : (b))
#define pow2(a) ((a)*(a))

void printMatrix(const double* mat, int r,int c)
{ 
  int i,j;

  for (j=0;j<r;j++)
    {
      for (i=0;i<c;i++)
        printf("\t%f ", (double) *(mat+j+r*i) );
      printf("\n");
    }  
  printf("\n");
}


int nbEdges = 0;

#define MAX_LINE 1000

double* pvL[MAX_LINE];
double* pvC[MAX_LINE];
double* vW = NULL;
int nbW = 0;;

double threshold =  0.003;



void fcn_noW(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
  //printf("%f %f\n", x[0], x[1]);
  double vVP_0 = cos(x[0])*sin(x[1]);
  double vVP_1 = sin(x[0])*sin(x[1]);
  double vVP_2 = cos(x[1]);
  int i =0;
  for (i=0;i<nbEdges;i++)
    {
      double* vL = pvL[i];
      double* vC = pvC[i];
      fvec[i] = vL[0]*vVP_0 + vL[1]*vVP_1 + vL[2]*vVP_2;
      fvec[i] /= sqrt(pow2(-vVP_2*vC[1] + vVP_1) + pow2(vVP_2*vC[0] - vVP_0));

    }

  return;
}

void fcn_W(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
  /* if (nbW <= 0){ */
  /*   printf("wrong W\n"); */
  /*   return;  */
  /* } */
    
  //assert(nbW>0);
  double vVP_0 = cos(x[0])*sin(x[1]);
  double vVP_1 = sin(x[0])*sin(x[1]);
  double vVP_2 = cos(x[1]);
  int i =0;
  for (i=0;i<nbEdges;i++)
    {
      double* vL = pvL[i];
      double* vC = pvC[i];
      fvec[i] = vL[0]*vVP_0 + vL[1]*vVP_1 + vL[2]*vVP_2;
      fvec[i] /= sqrt(pow2(-vVP_2*vC[1] + vVP_1) + pow2(vVP_2*vC[0] - vVP_0));

      /* //robustity */
      /* if (fvec[i] > threshold) */
      /* 	fvec[i] = threshold; */
      
	fvec[i] *= vW[i];
      
    }

  return;
}


void
mexFunction(int nout, mxArray *out[], 
            int nin, const mxArray *in[])
{
  enum {vsEdges_i=0,   vW_i} ;
  enum {vVP_i=0};


  // int i =cvUseOptimized(1);
  // printf("optimized -> %d\n", i);
  /* ------------------------------------------------------------------
  **                                                Check the arguments
  ** --------------------------------------------------------------- */ 
  if (nin >2 || nin < 1) {
    mexErrMsgTxt(" 1 or 2 arguments required");
    return;
  } 
  if (nout > 1) {
    mexErrMsgTxt("Too many output arguments");
    return;
  }

//structure
  const mxArray* vsEdges = in[vsEdges_i];  
  nbEdges  = mxGetN(vsEdges);
  
  int nbIdx = nbEdges;
 

  out[vVP_i] = (mxArray*) mxCreateDoubleMatrix(3,1,mxREAL) ; 
  double*  vVP    = (double*) mxGetPr(out[vVP_i]);


  if (nin>=2)
    {
      vW  = (double*) mxGetData(in[vW_i]) ;
      nbW = mxGetN(in[vW_i]);
 
      if (nbEdges != nbW) {
	mexErrMsgTxt("vsEdges and vW not same size");
	return;
      } 
   }

 
  if (nbIdx==2)
    {
      mxArray* vL_mx;
      vL_mx = mxGetField(vsEdges, 0,"vL");
      double* vL1     = mxGetPr(vL_mx);
      vL_mx = mxGetField(vsEdges, 1,"vL");
      double* vL2     = mxGetPr(vL_mx);

      vVP[0] = -vL1[2]*vL2[1] + vL1[1]*vL2[2];
      vVP[1] =  vL1[2]*vL2[0] - vL1[0]*vL2[2];
      vVP[2] = -vL1[1]*vL2[0] + vL1[0]*vL2[1];

      double nrm = sqrt(vVP[0]*vVP[0]+ vVP[1]*vVP[1] + vVP[2]*vVP[2]) ;
      vVP[0] /= nrm;
      vVP[1] /= nrm;
      vVP[2] /= nrm;

      return;
    }



  double mAtA_all[9] = {0.};

  int i =0;
  for (i = 0; i<nbIdx; i++)
    {
      mxArray* mAtA_mx;
      //mAtA_mx = mxGetField(vsEdges, i,"mAtA");
      mAtA_mx = mxGetField(vsEdges, i,"vLtL");
      
      double* mAtA     = mxGetPr(mAtA_mx);
      double w = 1;
      if (nin>=2)
      	{
      	  w = vW[i ];
      	  w*=w;
      	}
      
      mAtA_all[0]+= w * mAtA[0];
      mAtA_all[1]+= w * mAtA[1];
      mAtA_all[2]+= w * mAtA[2];
      mAtA_all[3]+= w * mAtA[3];
      mAtA_all[4]+= w * mAtA[4];
      mAtA_all[5]+= w * mAtA[5];
      mAtA_all[6]+= w * mAtA[6];
      mAtA_all[7]+= w * mAtA[7];
      mAtA_all[8]+= w * mAtA[8];
    }

  //printMatrix(mAtA_all,3,3);
  doublereal w[3];

  {
    char jobz = 'V';
    char UPLO = 'L';
    integer n = 3;
    integer lda = 3;
    doublereal work[50];
    integer lwork = 50;
    integer info;
    DSYEV(&jobz, &UPLO, &n, mAtA_all,
	   &lda, w, work, &lwork, 
	   &info);
  }
  w[0] =fabs(w[0]);
  w[1] =fabs(w[1]);
  w[2] =fabs(w[2]);

  //printMatrix(mAtA_all,3,3);
	      

  if (w[0] < w[1])
    {
      if (w[0] < w[2])
	{
	  vVP[0] = mAtA_all[0];
	  vVP[1] = mAtA_all[1];
	  vVP[2] = mAtA_all[2];
	}
      else
	{
	  vVP[0] = mAtA_all[6];
	  vVP[1] = mAtA_all[7];
	  vVP[2] = mAtA_all[8];
	}
    }
  else
    {
      if (w[1] < w[2])
	{
	  vVP[0] = mAtA_all[3];
	  vVP[1] = mAtA_all[4];
	  vVP[2] = mAtA_all[5];
	}
      else
	{
	  vVP[0] = mAtA_all[6];
	  vVP[1] = mAtA_all[7];
	  vVP[2] = mAtA_all[8];
	}

    }


  /*----------------------------------------------*/
  {    
    for (i = 0; i<nbIdx; i++)
      {
	mxArray* vL_mx;
	mxArray* vC_mx;
	vL_mx = mxGetField(vsEdges, i,"vL2");
	pvL[i]     = mxGetPr(vL_mx);
	vC_mx = mxGetField(vsEdges, i,"vCtE");
	pvC[i]     = mxGetPr(vC_mx);
      }

    //refine by NL
#define LWA 100
    int info;
    int lwa = LWA;
    double tol=0.0001  ;
    double wa[LWA];
    int iwa[2];
    //double fvec[3];
    int m = nbIdx;
    int n = 2;
    double* fvec = (double *) mxMalloc((nbIdx)*sizeof(double));
    
    double  vAngle[2];
    vAngle[0] = atan2(vVP[1], vVP[0]);
    vAngle[1] = acos(vVP[2]);

    //printf("%f %f %f\n", vVP[0], vVP[1], vVP[2]);
    

    if (nin>=2)
      lmdif1_(&fcn_W, &m, &n, vAngle, fvec, &tol, &info, iwa, wa, &lwa);
    else
      lmdif1_(&fcn_noW, &m, &n, vAngle, fvec, &tol, &info, iwa, wa, &lwa);
    
    vVP[0] = cos(vAngle[0])*sin(vAngle[1]);
    vVP[1] = sin(vAngle[0])*sin(vAngle[1]);
    vVP[2] = cos(vAngle[1]);
    //printf("%f %f %f\n", vVP[0], vVP[1], vVP[2]);
    
    

    double nrm = sqrt(vVP[0]*vVP[0]+ vVP[1]*vVP[1] + vVP[2]*vVP[2]) ;
    vVP[0] /= nrm;
    vVP[1] /= nrm;
    vVP[2] /= nrm;
    
    
    mxFree(fvec);
  }


  return;
  
}


