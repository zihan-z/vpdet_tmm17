// Jlinkage algorithm, 
// implemented by Jean-Philippe Tardif (tardifj@gmail.com)
// 19/05/2009
// This implementation is faster and less complxe than Toldo's
//
typedef unsigned int uint;
#include "mex.h"
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <algorithm>

#define DBL_MAX 1.7976931348623158e+308

#define MODULE "mexJLinkage"
#define MPRINTF(format, ...) {                  \
    mexPrintf ("[%s] ", MODULE);                \
    mexPrintf (format, ##__VA_ARGS__);          \
}


//#define DEBUG


//structure of a cluster
typedef struct {
  std::vector<uint>   vData;
  //std::vector<double> vDist;
  double* vDist;
  double minDist;
  uint minDistIdx;
  bool*   pBin;
  uint pLength;
  bool disable;


} Cluster;


//For debugging
void
checkDist(std::vector<Cluster> &vCluster)
{
  
  for (uint i=0; i < vCluster.size()-1; i++)
    {
      Cluster c = vCluster[i];

      if (c.disable) continue;

      double minDist = DBL_MAX;
      uint c2 = INT_MAX;
      for (uint j=0; j < vCluster.size(); j++)
	{	  
	  if (j<=i){
	    if (c.vDist[j] < DBL_MAX)
	      mexErrMsgIdAndTxt("stats:linkagemex:error",
				"non infinit value");
	    continue;
	  }
	  
	  //MPRINTF("\t(%d) %lf\n", j, c.vDist[j]);
	  if (c.vDist[j] < minDist)
	    {
	      c2 = j;
	      minDist = c.vDist[j];
	    }	  

	}
      if (minDist == DBL_MAX &&   c.minDist == DBL_MAX)
	continue;

      if (minDist != c.minDist){
	MPRINTF("%d) %d -> %lf, %d -> %lf \n", i, c2, minDist, c.minDistIdx, c.minDist);
	mexErrMsgIdAndTxt("stats:linkagemex:error",
			  "wrong minDist");
      }
      if (c2 != c.minDistIdx)	{
	MPRINTF("%d) %d -> %lf, %d -> %lf \n", i, c2, minDist, c.minDistIdx, c.minDist);
	mexErrMsgIdAndTxt("stats:linkagemex:error",
			  "wrong minDistIdx");
      }

    }
}



//compute jaccard distance between two clusters
double 
getClusterJaccardDist(Cluster &c1, Cluster &c2)
{
  double p=0;
  double s=0;
  
  bool* b1 = c1.pBin;
  bool* b2 = c2.pBin;

  for (uint i = 0; i < c1.pLength; i++){
    //Jlinkage
    p+= (double) (b1[i] & b2[i]);
    s+= (double) (b1[i] | b2[i]);
  }
  if (s>0)
    p/= s;
  else
    p = 0;

  return (1.0 - p);  
}

//intersects the two clusters, 
//put results in c1
inline void 
intersectClusters(Cluster &c1, Cluster &c2)
{
  
  for (uint j = 0; j < c1.pLength; j++)
    {
      c1.pBin[j] &= c2.pBin[j];
    }

}


//Find the two "closests" clusters 
double
getClosestClusters(std::vector<Cluster> &vCluster, uint &c1, uint &c2)
{
#ifdef DEBUG
  MPRINTF("getting closest clusters (%d)...", vCluster.size() );
  uint nbDisable =0;
  for (uint i=0; i < vCluster.size(); i++)
    if (vCluster[i].disable) nbDisable++;
  MPRINTF("nb disable= %d\n", nbDisable);
#endif
  double minDist = DBL_MAX;


  for (uint i=0; i < vCluster.size()-1; i++)
    {
      Cluster c = vCluster[i];
      
      if (c.disable) continue;
      
       if (c.minDist >= minDist) continue; //nothing closer to that cluster, skipping check

      minDist = c.minDist;
      c1 = i;
      c2 = c.minDistIdx;

      //below is useless
      // 	{
      // 	  //MPRINTF("\t(%d) %lf\n", j, c.vDist[j]);
      // 	  if (c.vDist[j] < minDist)
      // 	    {
      // 	      c1 = i;
      // 	      c2 = j;
      // 	      minDist = c.vDist[j];
      // 	    }	  
      // 	}
    }

  assert(c1<c2);
  

#ifdef DEBUG
  printf("%d %d -> %lf\n", c1,c2, minDist);
#endif
  return minDist;

}

//merges cluster c2 into c1
void
mergeClusters(std::vector<Cluster> &vCluster, uint c1, uint c2)
{
#ifdef DEBUG
  MPRINTF("Merging %d %d\n", c1,c2);
#endif

#ifdef DEBUG
  assert(c1<c2);
  
  if(vCluster[c2].disable )
    mexErrMsgIdAndTxt("stats:linkagemex:error",
		      "already disable");
  if(vCluster[c1].disable )
    mexErrMsgIdAndTxt("stats:linkagemex:error",
		      "Should not be disable");
#endif
  
  //merging:
  //insert data of c2 into c1
  vCluster[c1].vData.insert(vCluster[c1].vData.end(), 
			    vCluster[c2].vData.begin(), 
			    vCluster[c2].vData.end()  );

  //erase c2's data and disable this cluster
  //we don't remove from vCluster for efficiency
  vCluster[c2].vData.clear();
  vCluster[c2].disable = true;

#ifdef DEBUG
  MPRINTF("Cluster %d contains %d point\n", c1, 
	 vCluster[c1].vData.size() );
#endif

  //Compute new value of c1
  intersectClusters( vCluster[c1],  vCluster[c2] );

  //min Dist of c1 is now invalid, recompute below
  vCluster[c1].minDist    = DBL_MAX;
  vCluster[c1].minDistIdx = INT_MAX;
  
  //recompute distance of all cluster to c1
  for (uint i=0; i< vCluster.size(); i++)
    {
      //clear distance for removed (merged into) cluster
      vCluster[i].vDist[c2] = DBL_MAX;
      vCluster[c2].vDist[i] = DBL_MAX;

      if (i==c1) continue;
      if (vCluster[i].disable) continue;


      //distance is saved in "upper diagonal"
      uint imin = std::min(i,c1);
      uint imax = std::max(i,c1);

      //get distance
      double d = getClusterJaccardDist(vCluster[imin], vCluster[imax]);
      vCluster[imin].vDist[imax] = d;

      //save as min distance if it smaller than the current one
      if (d< vCluster[imin].minDist)
      	{
      	  vCluster[imin].minDist = d;
      	  vCluster[imin].minDistIdx = imax;
      	}
    }

  //clusters with minimal distance to c2 and c1, min dist must be reinitialized
  for (uint i=0; i< vCluster.size(); i++)
    {
      if (vCluster[i].minDistIdx != c2 &&
	  vCluster[i].minDistIdx != c1) continue; //nothing to do
      
      //recompute minimal distance
      vCluster[i].minDist    =  DBL_MAX;
      for (uint j=i+1;j<vCluster.size(); j++)
	{
	  double d = vCluster[i].vDist[j];// = d;
	  if (d< vCluster[i].minDist)
	    {
	      vCluster[i].minDist = d;
	      vCluster[i].minDistIdx = j;
	    }
	}

    }


}




void mexLinkageTEMPLATE(int nlhs,             /* Number of left hand side (output) arguments */
			mxArray *plhs[],      /* Array  of left hand side (output) arguments */
			int nrhs,             /* Number of right hand side (input) arguments */
			const mxArray *prhs[]) /* Array  of right hand side (input) arguments */
{

  /* get the dimensions of inputs */
  bool* totbin = (bool   *) mxGetData(prhs[0]);
  
  int nbIterMax = (int) mxGetScalar(prhs[1]);

  uint nMod = mxGetM(prhs[0]); //# of columns
  uint nPts = mxGetN(prhs[0]); //# of rows

  std::vector< Cluster > vCluster;
  
  //create initial clusters == nbPts
  for (uint i=0; i< nPts; i++)
    {
      Cluster c;
      c.disable = false;
      c.vData.push_back(i);                           //cluster contains only one point
      c.vDist = (double*) malloc(nPts*sizeof(double));

      c.pBin = (bool*) malloc(nMod*sizeof(bool));
      memcpy(c.pBin, totbin+i*nMod, nMod*sizeof(bool) ); //copy point in model space
      
      c.pLength = nMod;

      vCluster.push_back(c); //add to list
    }

#ifdef DEBUG
  MPRINTF("%d clusters\n", vCluster.size() );
#endif

  //initializing distances
  for (uint i=0; i< vCluster.size(); i++)
    for (uint j=0; j< vCluster.size(); j++)
      {
	vCluster[i].vDist[j] = DBL_MAX;
      }
  
  //compute distances between all clusters
  for (uint i=0; i< vCluster.size(); i++)
    {

      //minimal distance info to speed up computation
      double  minDist    = DBL_MAX;
      uint    minDistIdx = INT_MAX;

      for (uint j=i+1; j< vCluster.size(); j++)
	{
	  double d = getClusterJaccardDist(vCluster[i], vCluster[j]);
	  vCluster[i].vDist[j] = d;
	  
	  if (d<minDist)
	    {
	      minDist = d;
	      minDistIdx = j;
	    }
	}
      //min distance to speed up computation later
      vCluster[i].minDist    = minDist;
      vCluster[i].minDistIdx = minDistIdx;
    }


#ifdef DEBUG
  MPRINTF("Beginining linkage\n");
#endif

#ifdef DEBUG
  checkDist(vCluster);
#endif

  int nbIter=0;
  while (true)
    {
      uint c1 = 0,c2 = 0;
      double d = getClosestClusters(vCluster, c1, c2);
      
      if (d>=1)
	break;

      mergeClusters(vCluster, c1,c2);
      nbIter++;
      //if (nbIter>=nbIterMax)
      //break;


#ifdef DEBUG
      checkDist(vCluster);
#endif
    }
  MPRINTF("Done in %d iteration\n", nbIter);
  
  //plhs[0] = mxCreateNumericMatrix(nPts,1,mxGetClassID(prhs[0]),mxREAL);
  plhs[0] = mxCreateDoubleMatrix(nPts,1, mxREAL);
 /*  create pointers to the output matrix */
 double* vClusterIdx = (double *) mxGetData(plhs[0]);   /*leftmost  column */
 uint idx = 1;
 for (uint i=0; i< vCluster.size(); i++)
   {
     if (vCluster[i].disable) continue;
     
     std::vector<uint> vData = vCluster[i].vData;
     for (uint j=0;j<vData.size();j++)
       vClusterIdx[ vData[j] ] = (double) idx;
     idx++;
     
   }

#ifdef DEBUG
  MPRINTF("clean up\n");
#endif
  //clear stuff
  for (uint i=0; i< vCluster.size(); i++)
    {
      free(vCluster[i].vDist);
      free(vCluster[i].pBin);
      
    }
}



/* GATEWAY FUNCTION */
void mexFunction(int nlhs,               /* Number of left hand side (output) arguments */
		 mxArray *plhs[],        /*  Array  of left hand side (output) arguments */
		 int nrhs,               /* Number of right hand side (input) arguments */
		 const mxArray *prhs[] ) /* Array  of right hand side (input) arguments */
		 
{

//#warning init agument unchecked, FIXME
  /*  check for proper number of arguments */
  // if((nrhs!=1))
  //   mexErrMsgIdAndTxt("stats:linkagemex:OneInputsRequired",
  // 		      "One");
  if(nlhs>1)
    mexErrMsgIdAndTxt("stats:linkagemex:TooManyOutputArguments",
		      "Too many output arguments for linkagemex.");  

  /* check input type */
  if (!mxIsLogical(prhs[0]))
    mexErrMsgIdAndTxt("stats:linkagemex:UndefinedFunctionOnlyLogical",
		      "Function linkagemex is only defined for values of class 'logical'.");

  mexLinkageTEMPLATE(nlhs,plhs,nrhs,prhs);

} /* void mexFunction  */


