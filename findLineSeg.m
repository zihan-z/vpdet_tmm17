function outList = findLineSeg(cnList,w,h,devRatio,pixelUpper,pixelLower)
% break a contour into multiple line segments

nu = [-1 0 1 0];
nv = [0 -1 0 1];

eMap = zeros(h,w);
eMap(cnList) = 1;

np = length(cnList);
cnListSorted = [];

% find the vertex
have_vertex = 0;
for k = 1:np
    [I,J] = ind2sub([h,w],cnList(k));
    cc = 0;
    for t = 1:4   
        if (I+nu(t)) > 0 && (I+nu(t)) <= h && (J+nv(t)) > 0 && (J+nv(t)) <= w...
                && eMap(I+nu(t),J+nv(t)) == 1
            cc = cc+1;
        end
    end
    if cc == 1
        cnListSorted{1}(1) = cnList(k);
        have_vertex = 1;
        break;
    end
end
if have_vertex == 0 % randomly pick a vertex
    cnListSorted{1}(1) = cnList(1);
end
    
% sort the list
k = 0;
while k < length(cnListSorted{1}) % more points to explore
    k = k+1;
    [I,J] = ind2sub([h,w],cnListSorted{1}(k));
    for t = 1:4
        if (I+nu(t)) > 0 && (I+nu(t)) <= h && (J+nv(t)) > 0 && (J+nv(t)) <= w...
                && eMap(I+nu(t),J+nv(t)) == 1
            cnListSorted{1}(k+1) = sub2ind([h,w],I+nu(t),J+nv(t));
            eMap(I,J) = 0; % cannot go back
            break;
        end
    end   
end

% break recursively
ptr0 = 0;
ptr1 = 1;
outList =[];
outListCount = 0;
while ptr0 < ptr1
    ptr0 = ptr0 + 1;
    
    tempList = cnListSorted{ptr0};
    [cnValid, bkIdx] = validList(tempList,w,h,devRatio,pixelUpper,pixelLower);
    if cnValid == 1
        outListCount = outListCount + 1;
        outList{outListCount} = tempList;
    else
        ptr1 = ptr1 + 1;
        cnListSorted{ptr1} = tempList(1:bkIdx);
        ptr1 = ptr1 + 1;
        cnListSorted{ptr1} = tempList(bkIdx+1:end);
    end
end
end

function [cnValid, bkIdx] = validList(inList,w,h,devRatio,pixelUpper,pixelLower)

[I,J] = ind2sub([h,w],inList);

x1 = J(1);
x2 = J(end);
y1 = I(1);
y2 = I(end);

dist = abs((y2-y1)*J - (x2-x1)*I + x2*y1 -y2*x1) / ...
            sqrt((y2-y1)^2 + (x2-x1)^2);
        
[val,bkIdx] = max(dist); % find the max distance
if (val / sqrt((y2-y1)^2 + (x2-x1)^2) > devRatio && val > pixelLower) ||...
        val > pixelUpper
    cnValid = 0;
else
    cnValid = 1;
end
end
