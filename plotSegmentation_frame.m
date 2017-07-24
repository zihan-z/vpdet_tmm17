function [mClr] = plotSegmentation_frame(im,  vsEdges,  vClass, vbOutliers, framesize)

cx = round(size(im,2)/2);
cy = round(size(im,1)/2);

nbClass = max(vClass);
mClr = FACADE_colorMap(nbClass);

hold on
for c=1:nbClass
    vsEdges_s = [vsEdges(vClass==c)];
    for i=1:length(vsEdges_s)
        
        vPts =  [vsEdges_s(i).vPts_un];
        vPts = vPts - repmat([cx;cy]-[framesize(2);framesize(1)]/2,1,size(vPts,2));
        
        if isempty(vPts), continue; end;
        plot(vPts(1,:), vPts(2,:), 'Color',  mClr(c,:), 'LineWidth', 2 );
    end
end
hold off

if true
    STEP=1;
    COL = 0;
    MARKER = '.';
else
    STEP=4;
    COL = 0.3;
    MARKER = 'x';
end

hold on
%outliers
if  ~isempty(vbOutliers) && sum(vbOutliers)>=1
    vb = (vbOutliers == 1);
    %vsEdges_s = [vsEdges(vbOutliers)];
    vsEdges_s = [vsEdges(vb)];
    
    vPts =  [vsEdges_s.vPts_un];
    vPts = vPts - repmat([cx;cy]-[framesize(2);framesize(1)]/2,1,size(vPts,2));
    
    plot(vPts(1,1:STEP:end), vPts(2,1:STEP:end),'.', 'Color', [0,0,0], 'LineWidth', 2 );
end
hold off

drawnow
return
