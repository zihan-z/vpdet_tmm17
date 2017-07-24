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
void jacdist(bool *x, mwSize m, mwSize n, double *d)
{

    mwIndex i,j,k;
    int   theSum,nz;
    bool  *XI, *XJ, *XI0;
	
    /* Determine which rows of X have missing data.
     */

	// Create the waitbar
// 	mxArray *hw,*title,*percent;
// 	double *percentDoublePtr;
// 	mxArray *input_array[2];
// 	mxArray *output_array[1];

	mexPrintf("calculating distances...");
	// title = mxCreateString("Calculating distances ...");
// 	percent = mxCreateDoubleMatrix(1,1, mxREAL);
// 	percentDoublePtr = (double *)mxGetPr(percent);
// 	*percentDoublePtr = 0.0;

// 	input_array[0] = percent;
// 	input_array[1] = title;

	// mexCallMATLAB(1, output_array, 2, input_array, "waitbar");
	// hw = output_array[0];
// 	input_array[1] = hw;
	
    XI = x;
    for (i=0; i<m; i++) {
		// *percentDoublePtr = (double)i/(double)m;
		// update the waitbar
		// mexCallMATLAB(0, output_array, 2, input_array, "waitbar");
		

        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;

			theSum = 0;
			nz = 0;
			for (k=0;k<n;k++,XI++,XJ++){
				if ((*XI) || (*XJ)) {
					nz++;
					if ((*XI)!=(*XJ)) {
						theSum++;
					}
				}
			}
			if (nz) {
				*(d++) = (double)theSum/(double)nz;
			} else {
				*(d++) = 1.0;
			}
        }
    }

	// close the waitbar
// 	input_array[0] = hw;
// 	mexCallMATLAB(0, output_array, 1, input_array, "close");
// 	mxDestroyArray(percent);
// 	mxDestroyArray(title);
	mexPrintf("done\n");
}


/* the dispatcher function */
void distfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize  numCoords,numPoints;
    bool    *x;
	double	*d;

  /*  create a pointer to the input matrix y */
    x = (bool *)mxGetLogicals(prhs[0]);

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
    jacdist(x,numPoints,numCoords,d);
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
    if (mxIsLogical(prhs[0])) {
        distfun(nlhs, plhs, nrhs, prhs);
    } else {
        mexErrMsgIdAndTxt("stats:pdistmex:BadInputType",
                          "PDISTJaccardMEX only supports real Logicals data.");
    }
}
