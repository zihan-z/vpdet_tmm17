% Get vanishing point using JLinkage
function [vsVP,vCluster] = getVP_JLinkage(vsEdges, im, ARGS)

const = FACADE_const();

%set threshol
switch ARGS.ERROR
    case const.ERROR_DIST
        THRESH = ARGS.JL_ERROR_DIST / ARGS.imgS; %overide ARGS.ERROR_DIST
    case const.ERROR_NDIST
        THRESH = ARGS.ERROR_NDIST*2;
    case const.ERROR_ANGLE
        THRESH = ARGS.ERROR_ANGLE*2;
end

%get hypothesis matrix
totd =genHyp(ARGS,THRESH,vsEdges, ARGS.imgS, im, ARGS.JL_minNModel*10, ARGS.JL_minNModel);
if isempty(totd)
    vsVP = [];
    vCluster = zeros(1,length(vsEdges));
    return;
end

% Perform J-Linkage clusterization
[T, totdbin] = clusterPoints(totd, THRESH, ARGS.JL_ALGO);
vCluster = T(:)';
FCprintf('nb Class %d\n', max(vCluster) );

%fit VP to each cluster
for c=1:max(vCluster)
    vb = (vCluster==c);% & ~vbOutliers;
    vsEdges__ = [vsEdges(vb)];
    vsVP(c).VP=FACADE_mxFitVP_x_MAX( vsEdges__, [vsEdges__.nb] );
end

%---------------------------------------------------
%ordering based on nb of edges corresponding to vp
[vsVP, vCluster] = FACADE_orderVP_nbPts(vsVP, vsEdges, vCluster);
end

%------------------------------------------------------------------------------
function [totd] =genHyp(ARGS, THRESH, vsEdges,  imgS, im, nMaxTrials, minNSol)

const = FACADE_const();
sampleSize = ARGS.ALL_sampleSize;

%Pre-generate a big bunch of numbers, because its faster
vNUSamples = FACADE_getSamples(vsEdges, ARGS.ALL_samplingAlgo);
vIdxAll = randInt(1,length(vNUSamples),1,2*nMaxTrials*sampleSize);
vIdxAll = vNUSamples(vIdxAll);

vL = [vsEdges.halfLength]';
vL2 = vL*2;
%assert(all(vS<1));

vp_all = zeros(2,minNSol);
%------------------------------
cErr = [];
ii=1;
for i = 1:nMaxTrials
    
    vIdx = vIdxAll(3*i:3*i+sampleSize-1);
    vsEdges__ = [vsEdges(vIdx)];
    %find vanishing point
    vsVP_VP = FACADE_mxFitVP_x_MAX(vsEdges__);
    
    %error
    vErr = mxDistCtrToVP_RMS(vsEdges, vsVP_VP);
    
    %refine solution with inliers
    %-->not recommended/necessary
    
    switch ARGS.ERROR
        case const.ERROR_DIST
            %nothing
        case const.ERROR_NDIST
            vErr = vErr ./ vL2;
        case const.ERROR_ANGLE
            vErr = atan2(vErr,vL)*180/pi;
    end
    
    %min number of inlier
    if sum(vErr <= THRESH ) <  ARGS.ALL_minInliers
        continue;
    end
    
    vp_all(1,ii) = vsVP_VP(1)/vsVP_VP(3) * imgS;
    vp_all(2,ii) = -vsVP_VP(2)/vsVP_VP(3) * imgS;   
    
    %save that for later
    cErr{ii} = vErr;
    ii=ii+1;
    if ii==minNSol
        break;
    end
end

%join all cells with vectors
if isempty(cErr)
    totd = [];
else
    totd = [cErr{:}];
end
FCprintf('done (%d trial, %d sol)\n', i,ii);
end
    
