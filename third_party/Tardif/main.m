function [hor, err] = main(imgStr, algoNum, vp_association, lines, drawerror, Thr)

ind = 0;
if (nargin < 1)
    %imgStr = 'P1020833.jpg';
    
    imgStr = '008.jpg';
    %imgStr = 'P1040870.jpg';
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
tic
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

    
classThr = 0;
while ( (sum (vClass == classThr + 1) >= Thr) && (classThr <size(vsVP, 2)  ) )
    classThr = classThr + 1;
end
vbOutliers = vbOutliers | vClass > classThr;
vsVP = vsVP(1:classThr );


%   vsVP = vsVP(1:3);
 %   vbOutliers = vbOutliers | vClass > 3;
 
%draw_cumerr(vsVP, vClass, vsEdges, vbOutliers, algoNum);

toc
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


drawLinesToVP (point, size(im, 2), size(im, 1), classClr);


hor = drawHorizont2 (point, size(im, 2));

        filename= sprintf( 'fig_hor/hor%02d_2.jpg', algoNum);
        %saveas(pFigEdges,filename, 'jpg'); 
        print('-noui' ,'-r0', '-djpeg',  filename);
        %print('-r0', '-djpeg',  filename);
        FCprintf('Wrote %s\n', filename);
 
        
%hold on
%figure;
%plot( 50:100, 900:950, 'Color',  'Yellow', 'LineWidth', 2 );
%hold on;
%plot( 1:640, 1:640, 'Color',  'Yellow', 'LineWidth', 2 );


