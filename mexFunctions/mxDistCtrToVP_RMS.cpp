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

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define fabs_(a) ((a) <  (0) ? -(a) : (a))
#define pow2(a) ((a)*(a))

// Modified 12/27/2015: L1 to Least-square line fitting (RMS error)
void mexFunction(int nout, mxArray *out[],
        int nin, const mxArray *in[])
{
    enum {vsEdges_i=0,  vVP_i} ;
    enum {vMaxErr_i=0};
    
    int i, j, k, t1, t2;
    
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
    double minErr, curErr;
    
    for (i = 0; i<nbEdges; i++)
    {
        mxArray* vPts_mx   = mxGetField(vsEdges, i,"vPts0");
        int nPts  = mxGetN(vPts_mx);
        double* vPts     = (double*) mxGetPr(vPts_mx);
        
        minErr = 100000;
        t1 = 0;
        for (j = 0; j<nPts; j++)
        {
            // note: vVP[1] = -vVP[1]
            vL[0] = 1*vVP[1] + vPts[t1+1]*vVP[2];
            vL[1] = vVP[0] - vPts[t1]*vVP[2];
            vL[2] = -vPts[t1+1]*vVP[0] - vPts[t1]*vVP[1];
            
            //normlization
            double nrm12 = sqrt(vL[0]*vL[0] + vL[1]*vL[1]);
            
            vL[0] = vL[0] / nrm12;
            vL[1] = vL[1] / nrm12;
            vL[2] = vL[2] / nrm12;
            
            curErr = 0;
            t2 = 0;
            for (k = 0; k<nPts; k++)
            {
                curErr += pow2(vL[0] * vPts[t2] + vL[1] * vPts[t2+1] + vL[2]);
                t2 += 2;
            }
            curErr /= nPts;
            curErr = sqrt(curErr);
            minErr = min(minErr, curErr);
            
            t1 += 2;           
        }
        
        //computing error
        vMaxErr[i] = minErr;
    }    
    out[vMaxErr_i]  =  vMaxErr_mx;    
}
