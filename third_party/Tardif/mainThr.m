function [err, VPnum, thr] = mainThr(imgStr, algoNum, vp_association, lines, drawerror, iter)

ind = 0;
if (nargin < 1)
    %imgStr = 'P1020833.jpg';
    algoNum = 1;
end

if (nargin < 6)
         Thr = 0;
end

%includes
addpath(genpath('JLinkage'));
addpath 'lineSegDetect/'

ARGS = FACADE_ARGS_default();
const = FACADE_const();
ARGS.mKinv = [];

%arguments for Vanishing point detection
ARGS.plot = 0; 
ARGS.savePlot = false;%true;

ARGS.manhattanVP = true;
ARGS.manhattanVP = false;
ARGS.REF_remNonManhantan = true;
ARGS.ALL_calibrated = false;
ARGS.ERROR = const.ERROR_DIST;
%------------------------------------------------------------------------------
%                                Edges and VP
%------------------------------------------------------------------------------
%get data (edges)
%read image
im = imread(imgStr);

im = rgb2gray(im);
if ARGS.plot
  f1 = sfigure(1); clf; imshow(im);
end
ARGS.imgS = max(size(im));


%getting edges
[vsEdges,ARGS.imE] = FACADE_getEdgelets2(im, ARGS);

%get vp
ARGS.JL_ALGO=2;
[vsVP,vClass] = FACADE_getVP_JLinkage(vsEdges, im,ARGS);
%vsVP = vsVP(1:3);
vbOutliers = FACADE_getOutliers(ARGS,vsEdges, vClass, vsVP);           


%%%OUTLIER THRESHHOLD

thr = 5:5:70;
for Thr = 5:5:70
     f1 = sfigure(1); clf; imshow(im);
    
classThr = 0;
while (sum (vClass == classThr + 1) >= Thr)
    classThr = classThr + 1;
end
vbOutliers = vbOutliers | vClass > classThr;
vsVP = vsVP(1:classThr );


%draw_cumerr(vsVP, vClass, vsEdges, vbOutliers, algoNum);

if (drawerror)
    vsVP = vsVP(1:3);
    vbOutliers = vbOutliers | vClass > 3;
    draw_cumerr(vsVP, vp_association, lines, algoNum);
end

%ploting
point = zeros (size(vsVP, 2), 3);
for i = 1:size(vsVP, 2)
    point(i, 1:3) = vsVP(i).VP;
    point (i, 1) = point (i, 1) / point (i, 3) * ARGS.imgS;
    point (i, 2) = - point (i, 2) / point (i, 3) * ARGS.imgS; 
    point (i, 3) = 1;
end

classClr = FACADE_plotSegmentation(im,  vsEdges,  vClass, [], vbOutliers, algoNum); %-1->don't save

%draw_horizont (point, size(im, 2));

hor = drawHorizont2 (point, size(im, 2), algoNum);

drawLinesToVP (point, size(im, 2), size(im, 1), classClr);

        filename= sprintf( 'fig_hor/hor%02diter%dThr%d.jpg', algoNum, iter, Thr);
        %saveas(pFigEdges,filename, 'jpg'); 
        print('-noui' ,'-r0', '-djpeg',  filename);
        %print('-r0', '-djpeg',  filename);
        FCprintf('Wrote %s\n', filename);
    
%hold on
%figure;
%plot( 50:100, 900:950, 'Color',  'Yellow', 'LineWidth', 2 );
%hold on;
%plot( 1:640, 1:640, 'Color',  'Yellow', 'LineWidth', 2 );


width = size (im, 2);
st = load (sprintf ('%d.mat', algoNum));

l = st.l;

dist = zeros (1, 4);
x = 1;
p1 = [x; (-l(1)*x - l(3))/ l(2); 1];
x = width;
p2 = [x; (-l(1)*x - l(3))/ l(2); 1];

dist(1) = abs(dot(hor, p1));
dist(2) = abs(dot(hor, p2));

x = 1;
p1 = [x; (-hor(1)*x - hor(3))/ hor(2); 1];
x = width;
p2 = [x; (-hor(1)*x - hor(3))/ hor(2); 1];

dist(3) = abs(dot(l, p1));
dist(4) = abs(dot(l, p2));
ind = ind +1;
err(ind) = max (dist(:));
sprintf ('Thr= %d Err = %d ', Thr, err(ind));
VPnum(ind) = classThr;
end
