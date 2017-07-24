// distfn_cylinder.cpp : mex-function interface implentation file


#include "mex.h"
#include <math.h>

/* the dispatcher function */
void distfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize  numPoints;
    double    *points,*parameters;
	double	*d;
	double	*cylinder_axis_start = new double[3];
	double	*cylinder_axis_direction = new double[3];
	double	cylinder_radius;
	double *projected_point = new double[3];

  /*  create a pointer to the input matrix */
    parameters = (double *)mxGetData(prhs[0]);
	cylinder_axis_start[0] = parameters[0];
	cylinder_axis_start[1] = parameters[1];
	cylinder_axis_start[2] = parameters[2];
	cylinder_axis_direction[0] = parameters[3];
	cylinder_axis_direction[1] = parameters[4];
	cylinder_axis_direction[2] = parameters[5];
	cylinder_radius = parameters[6];

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
    norm = sqrt(pow(cylinder_axis_direction[0],2) + pow(cylinder_axis_direction[1],2) + pow(cylinder_axis_direction[2],2));
	cylinder_axis_direction[0] = cylinder_axis_direction[0]/norm;
	cylinder_axis_direction[1] = cylinder_axis_direction[1]/norm;
	cylinder_axis_direction[2] =cylinder_axis_direction[2]/norm;
	
	for (mwSize i = 0; i < numPoints; i++){
		// Project point on cylinder axis	
		t = (points[i * 3 + 0] - cylinder_axis_start[0]) * cylinder_axis_direction[0];
		t += (points[i * 3 + 1] - cylinder_axis_start[1]) * cylinder_axis_direction[1];
		t += (points[i * 3 + 2] - cylinder_axis_start[2]) * cylinder_axis_direction[2];
		

		projected_point[0] =  cylinder_axis_start[0] + cylinder_axis_direction[0] * t;
		projected_point[1] =  cylinder_axis_start[1] + cylinder_axis_direction[1] * t;
		projected_point[2] =  cylinder_axis_start[2] + cylinder_axis_direction[2] * t;
		norm = sqrt(pow(points[i * 3 + 0]-projected_point[0],2) + pow(points[i * 3 + 1]-projected_point[1],2) + pow(points[i * 3 + 2]-projected_point[2],2));
		
		// project to Cylinder surface
		projected_point[0] = projected_point[0] + (cylinder_radius * (points[i * 3 + 0]-projected_point[0]) / norm);
		projected_point[1] = projected_point[1] + (cylinder_radius * (points[i * 3 + 1]-projected_point[1]) / norm);
		projected_point[2] = projected_point[2] + (cylinder_radius * (points[i * 3 + 2]-projected_point[2]) / norm);

		// Compute distance beetween the two points
		norm = sqrt(pow(points[i * 3 + 0]-projected_point[0],2) + pow(points[i * 3 + 1]-projected_point[1],2) + pow(points[i * 3 + 2]-projected_point[2],2));

		d[i] = norm;
	}

	delete(cylinder_axis_start);
	delete(projected_point);


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

