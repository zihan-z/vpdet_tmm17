function vsEdges = getEdgeFromContour(img, ucm, ARGS)

ARGS.imgW =size(img,2);
ARGS.imgH =size(img,1);

edgelets = divideContour(ucm);
ne = length(edgelets);

% convert to subscripts
vn = zeros(1,ne);
edgelist = cell(1,ne);
for k = 1:ne
    [I,J] = ind2sub([ARGS.imgH,ARGS.imgW],edgelets{k});
    if size(I,1) == 1
        edgelist{k} = [J; I];
    else
        edgelist{k} = [J'; I'];
    end
    vn(k) = length(edgelets{k});
end
% remove short segments
edgelist = edgelist(vn>ARGS.minEdgeLength);
FCprintf('Found %d lines\n', length(edgelist));

vPts_un = edgelist;
vsEdges = build_vsEdges(vPts_un, ARGS, false);

return

