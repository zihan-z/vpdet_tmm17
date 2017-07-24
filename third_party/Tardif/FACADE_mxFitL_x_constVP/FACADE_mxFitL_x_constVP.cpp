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

#define END_POINTS

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define fabs_(a) ((a) <  (0) ? -(a) : (a))




void
mexFunction(int nout, mxArray *out[], 
            int nin, const mxArray *in[])
{
  enum {vsEdges_i=0,  vVP_i} ;
  enum {vMaxErr_i=0};


  // int i =cvUseOptimized(1);
  // printf("optimized -> %d\n", i);
  /* ------------------------------------------------------------------
  **                                                Check the arguments
  ** --------------------------------------------------------------- */ 
  if (nin != 2) {
    mexErrMsgTxt(" 2 arguments required");
  } 
  if (nout > 1) {
    mexErrMsgTxt("Too many output arguments");
  }


  double* vVP  = (double*) mxGetData(in[vVP_i]) ;

  int nVP  = mxGetN(in[vVP_i]); 
  if (nVP != 1) {
    mexErrMsgTxt(" vP must be a 3-vector");
  } 
  //structure
  const mxArray* vsEdges = in[vsEdges_i];
  
  int nbEdges  = mxGetN(vsEdges);
  mxArray* vMaxErr_mx = (mxArray*) mxCreateDoubleMatrix(nbEdges,1,mxREAL) ; 
  double*  vMaxErr    = (double*) mxGetPr(vMaxErr_mx);
  double vL[3];

  for (int i = 0; i<nbEdges; i++)
    {
#ifdef END_POINTS
      mxArray* skMat_mx = mxGetField(vsEdges, i,"skmVctE");
#else
      mxArray* skMat_mx = mxGetField(vsEdges, i,"skmVct");
#endif
      mxArray* vP1_mx   = mxGetField(vsEdges, i,"vPoint1");
      mxArray* vP2_mx   = mxGetField(vsEdges, i,"vPoint2");
      double* skMat     = mxGetPr(skMat_mx);
      double* vP1       = mxGetPr(vP1_mx);
      double* vP2       = mxGetPr(vP2_mx);
      
      //Computing line
      vL[0] = skMat[0]*vVP[0] + skMat[3]*vVP[1] + skMat[6]*vVP[2];
      vL[1] = skMat[1]*vVP[0] + skMat[4]*vVP[1] + skMat[7]*vVP[2];
      vL[2] = skMat[2]*vVP[0] + skMat[5]*vVP[1] + skMat[8]*vVP[2];

      //normlization
      double nrm12 = sqrt(vL[0]*vL[0] + vL[1]*vL[1]);
      
      vL[0] = vL[0] / nrm12;
      vL[1] = vL[1] / nrm12;
      vL[2] = vL[2] / nrm12;
      
      //computing error
      double err = vL[0]*vP1[0] + vL[1]*vP1[1] +  vL[2];
      vMaxErr[i] = fabs_(err);
      //double err1 = fabs(vL[0]*vP1[0] + vL[1]*vP1[1] +  vL[2]);
      //double err2 = fabs(vL[0]*vP2[0] + vL[1]*vP2[1] +  vL[2]);
      //double maxErr = max(err1,err2);
      //vMaxErr[i] = maxErr;

    }

  out[vMaxErr_i]  =  vMaxErr_mx;

  // double* vPh  = (double*) mxGetData(in[vPh_i]) ;

  // int vPh_nrow  = mxGetM(in[vPh_i]) ;
  // int vPh_ncol  = mxGetN(in[vPh_i]) ;

  // if (vPh_nrow != 2) {
  //   mexErrMsgTxt("3rd arg must be 3xN");
  // } 

  // double meanErr = 0;
  // maxErr = 0;
  // if (nout ==3)
  //   {
  //     int j =0;
  //     for (int i = 0; i< vPh_ncol; i++)
  // 	{
  // 	  double err = fabs(vL[0]*vPh[j] + vL[1]*vPh[j+1] +  vL[2]);
  // 	  maxErr = max(maxErr, err);
  // 	  meanErr += err;
  // 	  j+=2;
  // 	}
  //   }
  // else
  //   {
  //     mxArray* vErr_mx = (mxArray*) mxCreateDoubleMatrix(1,vPh_ncol,mxREAL) ; 
  //     double* vErr     = (double*) mxGetPr(vErr_mx);
      
  //     int j =0;
  //     for (int i = 0; i< vPh_ncol; i++)
  // 	{
  // 	  vErr[i] = fabs(vL[0]*vPh[j] + vL[1]*vPh[j+1] +  vL[2]);
  // 	  maxErr = max(maxErr, vErr[i]);
  // 	  meanErr += vErr[i];
  // 	  j+=2;
  // 	}
  //     out[vErr_i]    = vErr_mx; 
      
  //   }

  // meanErr /= vPh_ncol;

  //stuff to return
  //out[maxErr_i]  = mxCreateDoubleScalar(maxErr);  
  //out[meanErr_i] = mxCreateDoubleScalar(meanErr);

  
  


}
