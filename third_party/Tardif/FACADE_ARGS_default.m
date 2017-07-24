function ARGS = FACADE_ARGS_default(ARGS)
    
    
    const = FACADE_const();
    
    %ARGS.ALL_pixelThreshold =2;
    ARGS.ALL_calibrated = false;
    ARGS.REF_remNonManhantan = false; %remove , to used with Calibrated JLinkage
         
    ARGS.ALL_samplingAlgo = const.SAMPLING_EDGELENGTH;
    %--------------------------------------------------------------------
    ARGS.ALL_sampleSize  = 2; %# of edges to compute candidates
    ARGS.ALL_minInliers  = 5;%5; %min # of inliers to keep a solutio
    
    %--------------------------------------------------------------------
    %error functions
    
    ARGS.ERROR = const.ERROR_DIST; %regular distance
    ARGS.ERROR_DIST = 2;
  
    %ARGS.ERROR = const.ERROR_NDIST; %normalized distance
    ARGS.ERROR_NDIST = 0.020;%*2;
                             
    %ARGS.ERROR = const.ERROR_ANGLE; %angle, buggy?
    ARGS.ERROR_ANGLE = 1.5;
    %--------------------------------------------------------------------
    %J-Linkage
    ARGS.JL_minNModel          = 500;
    ARGS.JL_solveVP            = const.SOLVE_VP_MAX;  %don't touch  
    ARGS.JL_ERROR_DIST         = 2;    %overide the above threshold
    ARGS.JL_GRICselect         = false;%true;
           
    %--------------------------------------------------------------------
    %GRIC
    ARGS.GRIC_ERROR = const.ERROR_DIST; %normalized distance
    ARGS.GRIC_minNModel   = 500;
    ARGS.GRIC_sigmaThresh  = 1;
    ARGS.GRIC_K            = 2; %# of parameters of manifold
    ARGS.GRIC_D            = 0; %Dimension of manifold: 0, since this is a point
    if ARGS.GRIC_ERROR == const.ERROR_DIST
      %dimension of data: two end points
      ARGS.GRIC_R            = 4; 
    elseif  ARGS.GRIC_ERROR == const.ERROR_DIST;
      ARGS.GRIC_R            = 4; 
    end
    ARGS.GRIC_useTSSE      = false;

    %--------------------------------------------------------------------
    %randomHOUGH
    ARGS.RH_pixelThreshold     = 1;
    ARGS.RH_sampleSize         = 2;
    ARGS.RH_nbIter             = 5;
    
    %--------------------------------------------------------------------
    %EM
    ARGS.EM_imgSigma       = 1^2;%(2/ARGS.imgS)^2;
    ARGS.EM_gaussSphSigma  = 3.6637e-05;%0.004;%cos(89.9*pi/180)
    ARGS.EM_useThreshold   = false;%true;%false;
    ARGS.EM_pixelThreshold = 2;

   
    %--------------------------------------------------------------------
    %tests: runTest function
    ARGS.manhattanVP       = false;
    ARGS.useGaussSphere    = false;  %default, GS is turned off
    ARGS.refine            = true;   %refine using e.g. EM
    ARGS.selfCalib_nonlin  = false;  %not fair
    ARGS.fixedrand         = false;  %initialize random gen with the same val, for DEBUGGING
    
    %--------------------------------------------------------------------
    %edges
    ARGS.minEdgeLength     = 20;    %in pixels
    ARGS.linesegTOL        = 2;     %pixel tolerance when selecting straight edges    
    ARGS.edgeCache         = false; %precomputing and saving, or online
    ARGS.edgeGT            = false; %use ground truth for Yord DB
    ARGS.plot              = 0;     %show plot
    ARGS.savePlot          = false; %save plot
