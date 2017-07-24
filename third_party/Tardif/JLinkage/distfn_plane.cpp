// distfn_plane.cpp : mex-function interface implentation file



#include "mex.h"
#include <math.h>

/* the dispatcher function */
void distfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize  numPoints;
    double    *points,*parameters;
	double	*d;
	double	*plane_normal = new double[3];
	double	plane_d;

  /*  create a pointer to the input matrix */
    parameters = (double *)mxGetData(prhs[0]);
	plane_normal[0] = parameters[0];
	plane_normal[1] = parameters[1];
	plane_normal[2] = parameters[2];
	plane_d = parameters[3];

    points = (double *)mxGetData(prhs[1]);

  /*  get the dimensions of the matrix input */
    numPoints = mxGetN(prhs[1]);

    /* make sure that the distance matrix can be created, then create it.  doing
     * this in double remains exact except in the cases where we error out anyway. */
     plhs[0] = mxCreateDoubleMatrix((mwSize)numPoints,1, mxREAL);
    
  /*  create a pointer to a copy of the output matrix */
    d = (double *)mxGetData(plhs[0]);
	
	double t = 0;
	double norm = 0;
	

	for (mwSize i = 0; i < numPoints; i++){
		// Project point on plane normal	
		d[i] = fabs(plane_normal[0] * points[i * 3 + 0] + plane_normal[1] * points[i * 3 + 1] + plane_normal[2] * points[i * 3 + 2] + plane_d);
	}



}


/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /*  check for proper number of arguments */
    if (nrhs<2) {
        mexErrMsgIdAndTxt("stats:distfn_cylindermex:TooFewInputs",
                          "Two input arguments required.");
    } else if(nlhs>2) {
        mexErrMsgIdAndTxt("stats:distfn_cylindermex:TooManyOutputs",
                          "Too many output arguments.");
    }

  /* Check the type of the input array */
  /* Currently only works with double or single(float) */
    if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])) {
        distfun(nlhs, plhs, nrhs, prhs);
    } else {
        mexErrMsgIdAndTxt("stats:pdistmex:BadInputType",
                          "distfn_cylindermex only supports double data.");
    }
}




