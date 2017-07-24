function edgeList = divideContour(ucm)

% parameters
ucmThres = 0.1;
devRatio = 0.05;
pixelUpper = devRatio * 200;
pixelLower = 2;

ucm(ucm<=ucmThres) = 0;
vals = sort(unique(ucm(:)));

h = size(ucm,1);
w = size(ucm,2);

% find contour edges
edgeList = [];
for k = 2:length(vals)
    cnv2 = vals(k);
    tempMap = (ucm==cnv2);
    
    CC = bwconncomp(tempMap,4);
    for ii = 1:CC.NumObjects
        cnList = CC.PixelIdxList{ii};
        
        tempList = findLineSeg(cnList,w,h,devRatio,pixelUpper,pixelLower);
        edgeList = [edgeList tempList];
    end
end




