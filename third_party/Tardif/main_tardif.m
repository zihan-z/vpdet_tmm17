function [hor_tard, zenith_tard] = main_tardif(imgStr, pictNum, Thr, outputFolder, orth, pp_imcoords)

addpath ('../lineSegDetect/')
addpath ('../JLinkage/')

rgb_palette = [255, 0, 0;...
    0, 255, 0; ...
    0, 0, 255; ...
    128, 128, 0; ...
    128, 0, 128; ...
    0, 128, 128; ...
    255, 128, 0; ...
    255, 0, 128; ...
    128, 255, 0; ...
    128, 0, 255; ...
    0, 128, 255; ...
    0, 255, 128; ...
    0, 0, 0];

init = 'Tardif';
%init = 'our';
readfromfile = false;

splitP = 'thr';
%splitP = 'kmeans';

refineTardifMethod = 'EM';
%refineTardifMethod = 'threshold';
%refineTardifMethod = 'none';

detectTardifTriplet = orth;

horizonMethod = 'orthogonal';
%horizonMethod = 'lsfit';

orthogonalizePoints = 'tardif';
%orthogonalizePoints = 'our';

%maxImageSize = 840;


ind = 0;
if (nargin < 1)
    %imgStr = 'P1020833.jpg';
    %imgStr = '1156n.jpg';
    %imgStr = 'P1020177.jpg';
    %imgStr = 'MarkedBase/1033n.jpg';
    %
    imgStr = 'Urban/003.jpg';
    %imgStr = 'UrbanImagesBase/1001.jpg';
    
    %imgStr = 'test.jpg';
    pictNum = 003;
    Thr = 50;

    %algoNum = 1;
end



FCprintf('\n picture ¹%d\n', pictNum);
%includes
path = genpath('JLinkage');
addpath(path);
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

rgbim = imread(imgStr);
%if (max(size(rgbim)) > maxImageSize) 
%    rgbim = imresize(rgbim, maxImageSize/max(size(rgbim)));
%end;

im = rgb2gray(rgbim);
if true %ARGS.plot
  f1 = sfigure(1); clf; imshow(im);
end
ARGS.imgS = max(size(im));

width = size(im, 2);
height = size(im, 1);

principal_point = [width/2 -height/2] / max(size(im));


pointsfile = sprintf('%s\\VP%d.mat', outputFolder, pictNum);


if (~readfromfile)
    [vsEdges,ARGS.imE] = FACADE_getEdgelets2(im, ARGS);
    ARGS.JL_ALGO=2;
    tic;        

    [vsVP,vClass] = FACADE_getVP_JLinkage(vsEdges, im,ARGS);

    toc;
    vbOutliers = FACADE_getOutliers(ARGS,vsEdges, vClass, vsVP);               
    %save (pointsfile, 'vsVP', 'vsEdges', 'vClass', 'vbOutliers');

else
    st = load (pointsfile);
    vsVP = st.vsVP;
    vsEdges = st.vsEdges;
    vClass = st.vClass;

end

%output tardif
%st = load (pointsfile);
%vsVP = st.vsVP;
%vsEdges = st.vsEdges;
%vClass = st.vClass;
%vbOutliers = st.vbOutliers;

if (strcmp (refineTardifMethod, 'threshold'))

    classThr = 0;
    while ( (sum (vClass == classThr + 1) >= Thr) && (classThr <size(vsVP, 2)  ) )
        classThr = classThr + 1;
    end
    vbOutliers = vbOutliers | vClass > classThr;
    vsVP = vsVP(1:classThr );

    %ploting
  

    
else
    if (strcmp (refineTardifMethod, 'EM'))
     
        Lines = zeros(length(vsEdges), 6);

        for iLine = 1:length(vsEdges)
            x1 = vsEdges(iLine).vPts_un(1, 1);
            y1 = vsEdges(iLine).vPts_un(2, 1);
            x2 = vsEdges(iLine).vPts_un(1, end);
            y2 = vsEdges(iLine).vPts_un(2, end);
            theta = atan2(y2-y1, x2-x1);
            r =  sqrt((x1-x2)^2 + (y1-y2)^2);
            Lines(iLine, :) = [x1 x2 y1 y2 theta r];
        end;

        imSize = size(rgbim);
        normLines = normalizeLines(Lines, imSize(1:2));

        [VPS, sigma, p, vClass] = APPestimateVpInitialized(normLines(:, :), vClass, rgbim, Thr);
        idx = 1:size(Lines, 1);
        vsVP = [];
        for i = 1:size (VPS, 1)
               vsVP(i).VP(1, 1) = (VPS (i, 1)) / ARGS.imgS;
               vsVP(i).VP(2, 1) = - ((VPS (i, 2)) / ARGS.imgS); 
               vsVP(i).VP(3, 1) = 1;
        end
        %output Kosecka      
        
        %vClass(vClass == 0) = 13;
        
        %VanDirImg = draw_line_image_whole(rgbim, Lines(idx, :)', rgb_palette(vClass, :), 2);

        %for i = 1:size(VPS, 1)
        %   VanDirImg = SimpleDrawLine(VanDirImg, [imSize(2) / 2, ...
        %       VPS(i, 1), imSize(1) / 2, VPS(i, 2)]', ...
        %       rgb_palette(i, :), 2);
        %end;

        %sfigure(1); clf;
        %imshow(VanDirImg);
      
    end;
        
end;
       
%for i = 1:size (vsVP, 2)
%    VPS(i, 1:3) = vsVP(i).VP(1:3, 1);
%end;

if (detectTardifTriplet && size(vsVP, 2) >= 3)

    %vClass(vClass == 13) = 0;

    if (strcmp (orthogonalizePoints, 'our'))
    best_ang = 100;
    bestzenith = 1;
    besthor1 = 2;
    besthor2 = 3;            
    for vp1 = 1:size(VPS, 1)
        for vp2 = 1:size(VPS, 1)
            for vp3 = 1:size(VPS, 1)

                if (vp1 == vp2 || vp1 == vp3 || vp2 == vp3)
                    continue;
                end;

                zenith_direction = VPS(vp1, 1:2) - [width/2 height/2];
                cur_horizon = VPS(vp2, 1:2) - VPS(vp3, 1:2);

                norm1 = norm(zenith_direction); 
                norm2 = norm(cur_horizon);

                cross_product = sum(zenith_direction .* cur_horizon);

                ang = acos(cross_product / (norm1*norm2));                    
                ang = abs(ang-pi/2);

                if (ang < best_ang)
                    best_ang = ang;
                    bestzenith = vp1;
                    besthor1 = vp2;
                    besthor2 = vp3;
                end;

            end;
        end;
    end;

    if (best_ang < pi/6)
        VPS = VPS( [bestzenith besthor1 besthor2], :);
    else  
        
        support = [];
        for ivp = 1:size(VPS, 1)
            support(ivp) = sum(vClass == ivp);
        end;
        [i_max_support1 max_support] = max(support);
        
        if (i_max_support1 <= size(VPS, 1))         
            support(i_max_support1) = [];        
            [i_max_support2 max_support] = max(support);
            VPS = VPS( [i_max_support1 i_max_support2], :);                    
        end;
    end;
    else
      [vsVP, vClass] = FACADE_orderVP_Mahattan(vsVP, vsEdges, vClass);
    
    end
end;
        
%output result 

  VPS = zeros (size(vsVP, 2), 3);
    for i = 1:size(vsVP, 2)
        VPS(i, 1:3) = vsVP(i).VP;
        VPS (i, 1) = VPS (i, 1) / VPS (i, 3) * ARGS.imgS;
        VPS (i, 2) = - VPS (i, 2) / VPS (i, 3) * ARGS.imgS; 
        VPS (i, 3) = 1;
    end

    classClr = FACADE_plotSegmentation(im,  vsEdges,  vClass, [], vbOutliers, pictNum); %-1->don't save

    %draw_horizont (point, size(im, 2));

    addpath ('../CODE_visualization')
    drawLinesToVP (VPS, size(im, 2), size(im, 1), classClr);


if (strcmp (horizonMethod, 'lsfit'))

    hor_tard = drawHorizont2 (VPS, size(im, 2), 'Magenta');   
    
else
    if (strcmp (horizonMethod, 'orthogonal'))     
 
        if (strcmp (splitP, 'thr'))
            
            %vp = zeros(size(VPS, 1), 2);
            vp_imcoord = VPS(:, 1:2);
            [hor_points, zen_points] = splitPointsThrAngleAndY (vp_imcoord, width, height);
           
%             for i = 1:size(VPS, 1)
%                 vp (i, 1) = (VPS (i, 1)-width/2) / ARGS.imgS;
%                 vp (i, 2) = - ((VPS (i, 2)-height/2) / ARGS.imgS); 
%             end
%             [hor_points, zen_points] = splitPointsThrAngleAndY (vp);
%             
%             for i = 1:size(hor_points, 1)
%                 hor_points (i, 1) = hor_points (i, 1)*ARGS.imgS + width/2;
%                 hor_points (i, 2) = - (hor_points (i, 2)*ARGS.imgS) + height/2; 
%             end
%             
%             for i = 1:size(zen_points, 1)
%                 zen_points (i, 1) = zen_points (i, 1)*ARGS.imgS + width/2;
%                 zen_points (i, 2) = - (zen_points (i, 2)*ARGS.imgS) + height/2; 
%             end
        else
            [hor_points, zen_points] = splitPoints (VPS);
        end

                       
        weights = zeros(size(VPS, 1), 1);
        weights_hor = zeros(size(hor_points, 1), 1);
        
        for i=1:size(VPS, 1)
            weights(i) = sum (vClass == i);
        end;
        
        for i = 1:size (hor_points, 1)
            weights_hor(i) = weights (find( (VPS(:, 1)==hor_points(i, 1)) .* (VPS(:, 2)==hor_points(i, 2)) ));
        
        end
        %if ~isempty(zen_points)
        %    zen_index = find( (VPS(:, 1)==zen_points(1, 1)) .* (VPS(:, 2)==zen_points(1, 2)) );        
        %    weights(zen_index) = [];
        %end;
        
        width = size(rgbim, 2);
        height = size(rgbim, 1);
        
        hor_tard = drawHorizont_orthogonal(hor_points, zen_points, weights_hor, width, height, 'Magenta', pp_imcoords);
        zenith_tard = zen_points;
    end;
end;

    filename= sprintf( '%s\\%03djusttardthr%f.jpg', outputFolder, pictNum, Thr);
    %saveas(pFigEdges,filename, 'jpg'); 
    print('-noui' ,'-r0', '-djpeg',  filename);
    %print('-r0', '-djpeg',  filename);
    FCprintf('Wrote %s\n', filename);

    horfile = sprintf ('%s\\%03dhor_tardif.mat', outputFolder, pictNum);
    save (horfile, 'hor_tard', 'zenith_tard');







