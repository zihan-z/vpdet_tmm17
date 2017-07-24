/*******************************************************
%  Create hierarchical cluster tree Z using the J-Linkage algorith,
%
% Authors: R.Toldo A.Fusiello, department of computer science - University of Verona.
% Reference Paper: R. Toldo, A. Fusiello. Robust Multiple Structures Estimation with J-linkage. Proceeding of the European Conference on Computer Vision, 2008.
*******************************************************/


#include "mex.h"
#include <math.h>
#include <string.h>


#define ISNAN(a) (a != a)
#define MAX_NUM_OF_INPUT_ARG_FOR_PDIST 50



double jaccard(bool *x1, bool *x2, mwSize size) 
{ int countIn, countUn; 
  mwSize i;
  
  for (countUn = 0, countIn = 0, i = 0; i < size; i++){
    if ( x1[i] & x2[i])
      countIn++;
    if ( x1[i] | x2[i])
      countUn++;
  }
  if(countUn > 0)
    return (1.0 - ((double)countIn/(double)countUn));
  else
    return 1.0;
}

void intersect(bool *x1, bool *x2, mwSize size) 
{ mwSize i; 
  for (i = 0; i < size; i++){
    x1[i] = x2[i] = x1[i] & x2[i];
  }
}


template<class TEMPL>
void mexLinkageTEMPLATE(
			int nlhs,             /* Number of left hand side (output) arguments */
			mxArray *plhs[],      /* Array  of left hand side (output) arguments */
			int nrhs,             /* Number of right hand side (input) arguments */
			const mxArray *prhs[],/* Array  of right hand side (input) arguments */
			TEMPL classDummy
			)
{

  static TEMPL  inf;
  mwSize        m,m2,m2m3,m2m1,n,i,j,bn,bc,bp,p1,p2,q,q1,q2,h,k,l;
  mwSize        nk,nl,nkpnl,sT,N;
  mwSize        *obp,*scl,*K,*L,*obsLogi;
  TEMPL         *y,*yi,*s,*b1,*b2,*T;
  TEMPL         t1,t2,t3;
  mwSize	      nMod,nPts;
  bool		  *logi,*logimt;
  int		      b1int, b2int;

  /* get the dimensions of inputs */
  n  = mxGetN(prhs[0]);                /* number of pairwise distances --> n */
  m = (mwSize) ceil(sqrt(2*(double)n));   /* size of distance matrix --> m = (1 + sqrt(1+8*n))/2 */

  /*  create a pointer to the input pairwise distances */
  yi = (TEMPL *) mxGetData(prhs[0]);

  /* set space to copy the input */
  y =  (TEMPL *) mxMalloc(n * sizeof(TEMPL));
  memcpy(y,yi,n * sizeof(TEMPL));

  nMod = mxGetN(prhs[1]);

  nPts = mxGetM(prhs[1]);

  logimt =  (bool *) mxGetLogicals(prhs[1]);

  logi =  (bool *) mxMalloc(nMod * nPts * sizeof(bool));

  for(mwSize i=0;i<nMod;i++)
    for(mwSize j=0;j<nPts;j++)
      logi[j*nMod + i] = logimt[i*nPts + j];

  obp = (mwSize *) mxMalloc(m * sizeof(mwSize));

  obsLogi = (mwSize *) mxMalloc(m * sizeof(mwSize));

  /* calculate some other constants */
  bn   = m-1;                        /* number of branches     --> bn */
  m2   = m * 2;                      /* 2*m */
  m2m3 = m2 - 3;                     /* 2*m - 3 */
  m2m1 = m2 - 1;                     /* 2*m - 1 */

  inf  = (TEMPL) mxGetInf();     /* inf */

  scl = (mwSize *) mxMalloc(m * sizeof(mwSize));
  /* initialize obp and scl */
  for (i=0; i<m; obsLogi[i] = i, obp[i]=i, scl[i++]=1);

  /*  allocate space for the output matrix  */
  plhs[0] = mxCreateNumericMatrix(bn,3,mxGetClassID(prhs[0]),mxREAL);

  /*  create pointers to the output matrix */
  b1 = (TEMPL *) mxGetData(plhs[0]);   /*leftmost  column */
  b2 = b1 + bn;                        /*center    column */
  s  = b2 + bn;                        /*rightmost column */
  b1int = 0; b2int = 1;


  /* find the best value for N (size of the temporal vector of  */
  /* minimums) depending on the problem size */
  if      (m>1023) N = 512;
  else if (m>511)  N = 256;
  else if (m>255)  N = 128;
  else if (m>127)  N = 64;
  else if (m>63)   N = 32;
  else             N = 16;

  /* set space for the vector of minimums (and indexes) */
  T = (TEMPL *) mxMalloc(N * sizeof(TEMPL));
  K = (mwSize *) mxMalloc(N * sizeof(mwSize));
  L = (mwSize *) mxMalloc(N * sizeof(mwSize));

  /* set space for the obs-branch pointers  */
  obp = (mwSize *) mxMalloc(m * sizeof(mwSize));
  for (i=0; i<m; i++) obp[i]=i;

  sT = 0;  t3 = inf;

  // Create the waitbar
  //mxArray *hw,*title,*percent;
  //double *percentDoublePtr;
  //mxArray *input_array[2];
  //mxArray *output_array[1];

  mexPrintf("Clustering...");
  //title = mxCreateString("Clustering ...");
  //percent = mxCreateDoubleMatrix(1,1, mxREAL);
  //percentDoublePtr = (double *)mxGetPr(percent);
  //*percentDoublePtr = 0.0;

  //input_array[0] = percent;
  //input_array[1] = title;

  // mexCallMATLAB(1, output_array, 2, input_array, "waitbar");
  // hw = output_array[0];
  // input_array[1] = hw;


  for (bc=0,bp=m;bc<bn;bc++,bp++){

    // update the waitbar
    //*percentDoublePtr = (double)bc/(double)bn;
    // 	mexCallMATLAB(0, output_array, 2, input_array, "waitbar");

    /* *** MAIN LOOP ***
       bc is a "branch counter" --> bc = [ 0:bn-1]
       bp is a "branch pointer" --> bp = [ m:m+bc-1 ], it is used to point
       branches in the output since the values [0:m-1]+1 are reserved for
       leaves.
    */

    /*  NEW METHOD: Keeps a sorted vector "T" with the N minimum distances,
	at every branch iteration we only pick the first entry. Now the whole
	"y" is not searched at every step, we will search it again only when
	all the entries in "T" have been used or invalidated. However, we need
	to keep track of invalid distances already sorted in "y", and also
	update the index vectors "K" and "L" with permutations occurred in the
	half matrix "y"
    */

    /* cuts "T" so it does not contain any distance greater than any of the
       new distances computed when joined the last clusters ("t3" contains
       the minimum new distance computed in the last iteration). */
    for (h=0;((T[h]<t3) && (h<sT));h++);
    sT = h; t3 = inf;
    /* ONLY when "T" is empty it searches again "y" for the N minimum
       distances  */
    if (sT==0) {
      for (h=0; h<N; T[h++]=inf);
      p1 = ((m2m1 - bc) * bc) >> 1; /* finds where the matrix starts */
      for (j=bc; j<m; j++) {
	for (i=j+1; i<m; i++) {
	  t2 = y[p1++];
	  /*  this would be needed to solve NaN bug in MSVC*/
	  /*  if (!mxIsNaN(t2)) { */
	  if (t2 <= T[N-1]) {
	    for (h=N-1; ((t2 <= T[h-1]) && (h>0)); h--) {
	      T[h]=T[h-1];
	      K[h]=K[h-1];
	      L[h]=L[h-1];
	    } /* for (h=N-1 ... */
	    T[h] = t2;
	    K[h] = j;
	    L[h] = i;
	    sT++;
	  } /* if (t2<T[N-1]) */
	  /*}*/
	} /*  for (i= ... */
      } /* for (j= ... */
      if (sT>N) sT=N;
    } /* if (sT<1) */

    /* if sT==0 but bc<bn then the remaining distances in "T" must be
       NaN's ! we break the loop, but still need to fill the remaining
       output rows with linkage info and NaN distances
    */
    if (sT==0) break;


    /* the first entry in the ordered vector of distances "T" is the one
       that will be used for this branch, "k" and "l" are its indexes */
    k=K[0]; l=L[0]; t1=T[0];

    /* some housekeeping over "T" to inactivate all the other minimum
       distances which also have a "k" or "l" index, and then also take
       care of those indexes of the distances which are in the leftmost
       column */
    for (h=0,i=1;i<sT;i++) {
      /* test if the other entries of "T" belong to the branch "k" or "l"
	 if it is true, do not move them in to the updated "T" because
	 these distances will be recomputed after merging the clusters */
      if ( (k!=K[i]) && (l!=L[i]) && (l!=K[i]) && (k!=L[i]) ) {
	T[h]=T[i];
	K[h]=K[i];
	L[h]=L[i];
	/* test if the preserved distances in "T" belong to the
	   leftmost column (to be permutated), if it is true find out
	   the value of the new indices for such entry */
	if (bc==K[h]) {
	  if (k>L[h]) {
	    K[h] = L[h];
	    L[h] = k;
	  } /* if (k> ...*/
	  else K[h] = k;
	} /* if (bc== ... */
	h++;
      } /* if k!= ... */
    } /* for (h=0 ... */
    sT=h; /* the new size of "T" after the shifting */
 
    /* Update output for this branch, puts smaller pointers always in the
       leftmost column  */

    if (obp[k]<obp[l]) {
      *b1++ = (TEMPL) (obp[k]+1); /* +1 since Matlab ptrs start at 1 */
      *b2++ = (TEMPL) (obp[l]+1);
    } else {
      *b1++ = (TEMPL) (obp[l]+1);
      *b2++ = (TEMPL) (obp[k]+1);
    }
    *s++ =  t1;

    /* Updates obs-branch pointers "obp" */
    obp[k] = obp[bc];        /* new cluster branch ptr */
    obp[l] = bp;             /* leftmost column cluster branch ptr */

    /*
       Merges two observations/clusters ("k" and "l") by re-calculating new
       distances for every remaining observation/cluster and place the
       information in the row/col "l" */

    /*

      example:  bc=2  k=5  l=8   bn=11   m=12

      0
      1    N                             Pairwise
      2    N   N                         Distance
      3    N   N   Y                     Half Matrix
      4    N   N   Y   Y
      5    N   N  p1*  *   *
      6    N   N   Y   Y   Y   +
      7    N   N   Y   Y   Y   +   Y
      8    N   N  p2*  *   *   []  +   +
      9    N   N   Y   Y   Y   o   Y   Y   o
      10   N   N   Y   Y   Y   o   Y   Y   o   Y
      11   N   N   Y   Y   Y   o   Y   Y   o   Y   Y

      0   1   2   3   4   5   6   7   8   9   10   11


      p1 is the initial pointer for the kth row-col
      p2 is the initial pointer for the lth row-col
      *  are the samples touched in the first loop
      +  are the samples touched in the second loop
      o  are the samples touched in the third loop
      N  is the part of the whole half matrix which is no longer used
      Y  are all the other samples (not touched)

    */

    /* computing some limit constants to set up the 3-loops to
       transverse Y */
    q1 = bn - k - 1;
    q2 = bn - l - 1;

    /* initial pointers to the "k" and  "l" entries in the remaining half
       matrix */
    p1 = (((m2m1 - bc) * bc) >> 1) + k - bc - 1;
    p2 = p1 - k + l;

    /* Get the cluster cardinalities  */
    nk     = scl[k];
    nl     = scl[l];
    nkpnl  = nk + nl;

    /* Updates cluster cardinality "scl" */
    scl[k] = scl[bc];        /* letfmost column cluster cardinality */
    scl[l] = nkpnl;          /* new cluster cardinality */


    intersect(&logi[nMod * (obsLogi[k])], &logi[nMod * (obsLogi[l])], nMod);

    mwSize col = bc;
    mwSize row = l;
			
    for (q=bn-bc-1; q>q1; q--) {
      t2 = jaccard(&logi[nMod * (obsLogi[col])], &logi[nMod * (obsLogi[row])] ,nMod);
      if (t2 < t3) t3 = t2 ;
      y[p2] = t2;
      p1 = p1 + q;
      p2 = p2 + q;
      col++;
    }
    col++;
    p1++;
    p2 = p2 + q;
    for (q=q1-1;  q>q2; q--) {
      t2 = jaccard(&logi[nMod * (obsLogi[col])], &logi[nMod * (obsLogi[row])] ,nMod);
      if (t2 < t3) t3 = t2 ;
      y[p2] = t2;
      p1++;
      p2 = p2 + q;
      col++;
    }
    p1++;
    p2++;
    col = l;
    row++;
    for (q=q2+1; q>0; q--) {
      t2 = jaccard(&logi[nMod * (obsLogi[col])], &logi[nMod * (obsLogi[row])] ,nMod);
      if (t2 < t3) t3 = t2 ;
      y[p2] = t2;
      p1++;
      p2++;
      row++;
    }
	     
    obsLogi[k] = obsLogi[bc];
	     
    /* 
       moves the leftmost column "bc" to row/col "k" */
    if (k!=bc) {
      q1 = bn - k;
	       
      p1 = (((m2m3 - bc) * bc) >> 1) + k - 1;
      p2 = p1 - k + bc + 1;

      for (q=bn-bc-1; q>q1; q--) {
	p1 = p1 + q;
	y[p1] = y[p2++];
      }
      p1 = p1 + q + 1;
      p2++;
      for ( ; q>0; q--) {
	y[p1++] = y[p2++];
      }
    } /*if (k!=bc) */

  } /*for (bc=0,bp=m;bc<bn;bc++,bp++) */

  /* loop to fill with NaN's in case the main loop ended prematurely */
  for (;bc<bn;bc++,bp++) {
    k=bc; l=bc+1;
    if (obp[k]<obp[l]) {
      *b1++ = (TEMPL) (obp[k]+1);
      *b2++ = (TEMPL) (obp[l]+1);
    } else {
      *b1++ = (TEMPL) (obp[l]+1);
      *b2++ = (TEMPL) (obp[k]+1);
    }
    obp[l] = bp;
    *s++ = (TEMPL) mxGetNaN();
  }

  mxFree(y);                               /* ... delete y mem */
  mxFree(logi);
  mxFree(scl);
  mxFree(obp);
  mxFree(L);
  mxFree(K);
  mxFree(T);
  // close the waitbar
  // input_array[0] = hw;
  // mexCallMATLAB(0, output_array, 1, input_array, "close");
  mexPrintf("end\n");
  //mxDestroyArray(percent);
  //mxDestroyArray(title);
}

void mexFunction(         /* GATEWAY FUNCTION */
		 int nlhs,             /* Number of left hand side (output) arguments */
		 mxArray *plhs[],      /* Array  of left hand side (output) arguments */
		 int nrhs,             /* Number of right hand side (input) arguments */
		 const mxArray *prhs[] /* Array  of right hand side (input) arguments */
			  )
{

  /*  check for proper number of arguments */
  if((nrhs!=2) && (nrhs!=3))
    mexErrMsgIdAndTxt("stats:linkagemex:TwoInputsRequired",
		      "Two");
  if(nlhs>1)
    mexErrMsgIdAndTxt("stats:linkagemex:TooManyOutputArguments",
		      "Too many output arguments for linkagemex.");

  /* check input type */
  if (!mxIsDouble(prhs[0]))
    mexErrMsgIdAndTxt("stats:linkagemex:UndefinedFunctionOnlyDouble",
		      "Function linkagemex is only defined for values of class 'double'.");

  /* check input type */
  if (!mxIsLogical(prhs[1]))
    mexErrMsgIdAndTxt("stats:linkagemex:UndefinedFunctionOnlyLogical",
		      "Function linkagemex is only defined for values of class 'logical'.");

  mexLinkageTEMPLATE(nlhs,plhs,nrhs,prhs,(double)(1.0));

} /* void mexFunction  */



