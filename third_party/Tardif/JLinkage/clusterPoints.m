%-----------------------------------------------------------
% Cluster Points using the J-Linkage algorithm
%
% Usage:
%
% [T, totdbin] = clusterPoints(Points,inliersThreshold,totd)
%
% Arguments:
%     totd             - Point-Model distance matrix
%     inliersThreshold - Ransac inlier Threshold value.
%
% Returns:
%     T         - Clustered points Labels
%     Z         - Hierarchical clustering tree
%     Y         - Pairwise Jaccard Distance
%     totdbin   - Consensus/Preference set matrix
%
% Authors: R.Toldo A.Fusiello, department of computer science - University of Verona.
% Reference Paper: R. Toldo, A. Fusiello. Robust Multiple Structures Estimation with J-linkage. Proceeding of the European Conference on Computer Vision, 2008.
%
%Modified by Jean-Philippe Tardif
%
%-----------------------------------------------------------
function [T, totdbin] = clusterPoints(totd, inliersThreshold, ALGO)
    
if nargin <=2
  ALGO = 2;
end

NBITERMAX = 2000000;
%ALGO = 3;

totdbin = totd < inliersThreshold;



%my own implementation
t=tic;
  switch ALGO
   case 1
    %P linkage
    totdIn=totd(totdbin);
    v = var([totdIn(:); -totdIn(:)]); %computing sigma
% $$$     %P linkage
    %totp = 1/sqrt(2*pi*v)*exp(-(totd.^2)/2/v);
    totp = exp(-(totd.^2)/2/v);
   
    
% $$$     sfigure(104);
% $$$     subplot(1,3,1)
% $$$     totp(totdbin);
% $$$     hist(ans(:),20)
% $$$     subplot(1,3,2)
% $$$     totp(~totdbin);
% $$$     hist(ans(:),20)
% $$$     subplot(1,3,3)
% $$$     totd(totdbin);
% $$$     hist(ans(:),20)
% $$$     pause
    
    
    totdbin = totd < inliersThreshold;
    totp(~totdbin) =0;
    
    T = mexPLinkage(totp', totdbin', NBITERMAX);
   case 2
    %standard Jlinkage
    totd = double(totdbin);
    T = mxJLinkage(totdbin', NBITERMAX);
  case 3
    %Toldo's implementation
    Y = pDistJaccard(totdbin');
    Z = linkageIntersect(Y, totdbin);
    T = cluster(Z,'cutoff',1-(1/(size(totdbin,1)))+eps,'criterion','distance');
  end
  
  
  Z = [];%linkageIntersect(Y, totdbin);
fprintf('done jp: %f \n', toc(t) );
return

%original implementation of The authors
t2=tic;
  Y = pDistJaccard(totdbin');
  Z = linkageIntersect(Y, totdbin);
  T2 = cluster(Z,'cutoff',1-(1/(size(totdbin,1)))+eps,'criterion','distance');
  fprintf('done orig: %f \n', toc(t2) );
  
