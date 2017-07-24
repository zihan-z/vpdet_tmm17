/*******************************************************
% Calculates pairwise jaccard distance
%
% Authors: R.Toldo A.Fusiello, department of computer science - University of Verona.
% Reference Paper: R. Toldo, A. Fusiello. Robust Multiple Structures Estimation with J-linkage. Proceeding of the European Conference on Computer Vision, 2008.
*******************************************************/

#include "mex.h"
#include <math.h>
#include <string.h>

/* Jaccard distance */
void pdist(double *x, mwSize m, mwSize n, double *d)
{

    mwIndex i,j,k;
    int   theSum;
    double  *XI, *XJ, *XI0;
	
    /* Determine which rows of X have missing data.
     */
    mexPrintf("calculating distances...");
    
    XI = x;
    for (i=0; i<m; i++) {
      
      XI0 = XI;
      XJ = XI+n;
      for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
	XI = XI0;
	
	theSum = 0;
	for (k=0;k<n;k++,XI++,XJ++){
	  theSum += (*XI) * (*XJ);
	}
	*(d++) = (double)theSum/(double)n;
      }
    }
   
    mexPrintf("done\n");
}


/* the dispatcher function */
void distfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize  numCoords,numPoints;
  double *x;
  double *d;
  
  /*  create a pointer to the input matrix y */
  x = (double *)mxGetData(prhs[0]);

  /*  get the dimensions of the matrix input y */
    numCoords = mxGetM(prhs[0]);
    numPoints = mxGetN(prhs[0]);

    /* make sure that the distance matrix can be created, then create it.  doing
     * this in double remains exact except in the cases where we error out anyway. */
    double numDists = ((double)numPoints * (double)(numPoints-1)) / 2;
    if (numDists >= (double)MWSIZE_MAX) {
        mexErrMsgIdAndTxt("stats:pdistmex:OutputTooLarge",
                          "Distance matrix has more elements than the maximum allowed size in MATLAB.");
    }
     plhs[0] = mxCreateDoubleMatrix(1, (mwSize)numDists, mxREAL);
    
  /*  create a pointer to a copy of the output matrix */
    d = (double *)mxGetData(plhs[0]);
    pdist(x,numPoints,numCoords,d);
 }

/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /*  check for proper number of arguments */
    if (nrhs<1) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooFewInputs",
                          "Two input arguments required.");
    } else if(nlhs>1) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooManyOutputs",
                          "Too many output arguments.");
    }

  /* Check the type of the input array */
  /* Currently only works with double or single(float) */
    if (mxIsDouble(prhs[0])) {
        distfun(nlhs, plhs, nrhs, prhs);
    } else {
        mexErrMsgIdAndTxt("stats:pdistmex:BadInputType",
                          " only supports real data.");
    }
}
