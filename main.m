clear;
clc;

addpath(genpath(('third_party/Tardif')));
addpath('third_party/BSR/grouping/lib');
addpath('mexFunctions');

imgFolder = 'data';
outFigureFolder = 'results';

% ------------------------------------------------------------------------------
%                                Parameters
%------------------------------------------------------------------------------
ARGS = FACADE_ARGS_default();
const = FACADE_const();
ARGS.mKinv = [];

framesize = [1000, 1000]; % range in which a VP is considered; default image size is 500 pixels
Thr = 2; % minimum points in a cluster

ARGS.minEdgeLength = 40;
ARGS.ALL_minInliers  = 2;

ARGS.JL_minNModel = 2000;
ARGS.JL_ERROR_DIST = 3;

display = true;
SAVE_FIGURE = 1;
SAVE_RESULT = 1;

% ------------------------------------------------------------------------------
%                                Edges and VP
%------------------------------------------------------------------------------

D = dir(fullfile(imgFolder, '*.jpg'));
for imgIdx = 1:length(D)
    imgName = D(imgIdx).name;
    disp(imgName);
    
    im = imread(fullfile(imgFolder, imgName)); 
    ARGS.imgS = max(size(im));   
    width = size(im, 2);
    height = size(im, 1);
    
    colorImg = im;
    if size(im,3)>1
        im = rgb2gray(im);
    end
    
    % contour detection
    disp('Detecting contour...');
    gPb_orient = globalPb(fullfile(imgFolder, imgName));
    ucm = contours2ucm(gPb_orient, 'imageSize');
    
    % contour to line segments
    vsEdges = getEdgeFromContour(im, ucm, ARGS);   
    if isempty(vsEdges) % no line segment detected
        continue;
    end
    
    % J-Linkage clustering
    ARGS.JL_ALGO=2;
    [vsVP, vClass] = getVP_JLinkage(vsEdges, im,ARGS);
    
    % remove outliers
    vbOutliers = FACADE_getOutliers(ARGS,vsEdges, vClass, vsVP);
    
    classThr = 0;
    while ( (sum (vClass == classThr + 1) >= Thr) && (classThr <size(vsVP, 2)  ) )
        classThr = classThr + 1;
    end
    vbOutliers = vbOutliers | vClass > classThr;
    vsVP = vsVP(1:classThr );
    
    % dominant VP selection
    % VPs already sorted by EDGE NUM
    im_coords = zeros(2, size(vsVP,2));
    IN = zeros(1,size(vsVP,2));
    weights = zeros(1,size(vsVP,2));
    
    xo = floor(framesize(2)/2 - width/2);
    yo = floor(framesize(1)/2 - height/2);
    
    outIdx = 0;
    for i = 1:size(vsVP, 2)
        x = vsVP(i).VP(1)/vsVP(i).VP(3);
        y = vsVP(i).VP(2)/vsVP(i).VP(3);
        x = x * ARGS.imgS;
        y = - y * ARGS.imgS;
        im_coords(:,i) = [x;y];
        
        x_ = x+xo;
        y_ = y+yo;
        
        if x_ > 0 && x_ <= framesize(2) && y_ > 0 && y_ <= framesize(1)
            IN(i) = 1;
            if outIdx == 0
                outIdx = i;
            end
        end
        weights(i) = sum(vClass == i);
    end
    
    % display the results
    if (display)
        if size(colorImg,3)>1
            colorImg_frameR = ones(framesize);
            colorImg_frameG = ones(framesize);
            colorImg_frameB = ones(framesize);
            
            colorImg_frameR(yo+1:yo+height,xo+1:xo+width)= double(colorImg(:,:,1))/255;
            colorImg_frameG(yo+1:yo+height,xo+1:xo+width)= double(colorImg(:,:,2))/255;
            colorImg_frameB(yo+1:yo+height,xo+1:xo+width)= double(colorImg(:,:,3))/255;
            
            colorImg_frame = cat(3,colorImg_frameR,colorImg_frameG,colorImg_frameB);
        else
            colorImg_frame = zeros(framesize);
            colorImg_frame(yo+1:yo+height,xo+1:xo+width)= double(colorImg)/255;
        end
          
        figure(2);clf;
        imshow(colorImg_frame);
        
        % display line segments
        hold on;
        classClr = plotSegmentation_frame(im,  vsEdges,  vClass, vbOutliers, framesize);
        hold off;

        % display VPs
        hold on;
        for i = 1:size(vsVP, 2)
            plot (im_coords(1,i)+xo, im_coords(2,i)+yo, 'x', 'Color',  classClr(i, :), 'LineWidth', 3, 'MarkerSize',16);
        end
        if outIdx ~= 0
            plot (im_coords(1,outIdx)+xo, im_coords(2,outIdx)+yo, 'o', 'Color',  classClr(outIdx, :), 'LineWidth', 3, 'MarkerSize',16);
        end
        hold off;
        
        if SAVE_FIGURE
            saveas(2,fullfile(outFigureFolder, imgName));
        end
    end
        
    if SAVE_RESULT
        vp_detected = [];
        if outIdx ~= 0
            vp_detected = im_coords(:,outIdx);
        end
        save(fullfile(outFigureFolder, [imgName(1:end-4) '.mat']),...
            'vp_detected', 'im_coords', 'vsEdges', 'vClass');
    end
    
    pause(.01);
end