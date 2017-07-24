function [vsEdges] = edges2vsEdges (lines, ARGS, im)
for i = 1:size (lines, 1)
    p1 = zeros (3, 1);
    p2 = zeros (3, 1);
    
    p1(1:2, 1) = lines (i, 1:2);
    p2(1:2, 1) = lines (i, 3:4);
    p1(3, 1) = 1;
    p2(3, 1) = 1;
    
    l = cross (p1, p2);
    if (abs (p1(1) - p2(1)) > abs(p1(2) - p2(2)))
        x = min (p1 (1), p2(1)) : max (p1 (1), p2(1));
        y = (-l(3) -l(1) * x)/l(2);
    else
        y = min (p1 (2), p2(2)) : max (p1 (2), p2(2));
        x = (-l(3) -l(2) * y)/l(1);
    end
    vPts_un{i} = zeros (2, size (x, 2));
    vPts_un{i}(1, :) = x;
    vPts_un{i}(2, :) = y;
    
end
ARGS.imgW =size(im,2);
ARGS.imgH =size(im,1);
vsEdges = FACADE_build_vsEdges(vPts_un, ARGS, false);
end