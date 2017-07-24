function vNUSamples = FACADE_getSamples(vsEdges, algo)

const = FACADE_const();


switch algo
 
 case const.SAMPLING_UNIFORM
  vNUSamples = 1:length(vsEdges);
  
 case const.SAMPLING_EDGELENGTH
  %generate non-uniform sampling, based on the length of the edges
  nbPts = sum([vsEdges.nb]);
  vNUSamples = zeros(1,nbPts); %non-uniform sampling
  ii=1;
  for i=1:length(vsEdges)
    nb =vsEdges(i).nb;
    vNUSamples(ii:ii+nb-1) = i;
    ii=ii+nb;
  end
  
 case const.SAMPLING_EDGEANGLE
  
  %orientation histogram
  vTheta = [vsEdges.theta];
  %sfigure(100);clf
  %hist(vTheta, 10);
  [N,X] = hist(vTheta, 10); %hirtogram of 
  [vH,vG] =  min(abs(repmat(X(:),1,length(vTheta)) - repmat(vTheta,10,1)));
  vG = vG(:)';
  mxH = max(N);
  tmp = {};
  for i=1:length(X)
    vI = vG==i;
    if sum(vI)==0, continue;end
    vIedge = find(vI);
    rep = ceil(mxH/sum(vI));
    tmp{i} = repmat(vIedge,1,rep);
  end
  vNUSamples = [tmp{:}];
  %vIdx2 = randInt(1,length(vIdx),1,2*nMaxTrials*sampleSize);
  %vNUSamples = vIdx(vIdx2);
  %vTheta = [vsEdges(vIdxAll).theta];
  %sfigure(101);clf
  %hist(vTheta, 10);
  %pause
  
end