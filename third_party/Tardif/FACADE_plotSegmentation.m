function [mClr] = FACADE_plotSegmentation(im,  vsEdges,  vClass, vP_v, vbOutliers, algoNum, noIm)
    
    %only3 = true;false;
    only3 = false;
    if nargin <= 5
        algoNum = 0; %no save
    end 
    if nargin <= 6
        noIm = false; %no save
    end
    
    
% $$$     if false %ARGS.only3                
% $$$         vSupport = hist(vClass, 1:max(vClass));
% $$$         [dum, vIdx] = sort(vSupport, 'descend');
% $$$         vClassCpy = vClass;
% $$$         for c=1:length(vIdx)
% $$$             vClassCpy(vClass==vIdx(c)) =c;
% $$$         end
% $$$         vClass = vClassCpy;
% $$$         %vClass(vClass>3) = 0; %eliminate outliers
% $$$     end
    %if ARGS.only3 
    %    vClass(vClass>3) = 0;
    %end
    
    %FCprintf('%d outliers\n', sum(vbOutliers))
    
    nbClass = max(vClass);
    mClr = FACADE_colorMap(nbClass);
    
    %ploting this
    sfigure(1); %   clf;
    %if ~noIm
    %    imshow(im);
    %end
    %edges with color/association
    hold on
    for c=1:nbClass        
        vsEdges_s = [vsEdges(vClass==c)];
        for i=1:length(vsEdges_s)
            if noIm
                vPts =  [vsEdges_s(i).vPts];            
            else
                vPts =  [vsEdges_s(i).vPts_un];            
            end
            if isempty(vPts), continue; end;        
            plot(vPts(1,:), vPts(2,:), 'Color',  mClr(c,:), 'LineWidth', 3 );
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
    
    
    if only3
        %flag last class as outliers
        hold on;
        vsEdges_s = [vsEdges(vClass>3)];
        if noIm
            vPts =  [vsEdges_s.vPts];      
        else
            vPts =  [vsEdges_s.vPts_un];      
        end
        
        if ~isempty(vPts)            
            plot(vPts(1,1:STEP:end), vPts(2,1:STEP:end),MARKER, 'Color', [COL,COL,COL], 'LineWidth', 3  );
        end
        hold off
    end
    
    hold on
    %outliers
    if  ~isempty(vbOutliers) && sum(vbOutliers)>=1
        vb = (vbOutliers == 1);
        %vsEdges_s = [vsEdges(vbOutliers)];
        vsEdges_s = [vsEdges(vb)];
        if noIm
            vPts =  [vsEdges_s.vPts];      
        else
            vPts =  [vsEdges_s.vPts_un];      
        end
        plot(vPts(1,1:STEP:end), vPts(2,1:STEP:end),'.', 'Color', [0,0,0], 'LineWidth', 3 );
    end    
    hold off
    
    
    if false %algoNum > 0 %save plot
        filename= sprintf( 'fig/res_class_%02d.jpg', algoNum);
        %saveas(pFigEdges,filename, 'jpg'); 
        print('-noui' ,'-r0', '-djpeg',  filename);
        %print('-r0', '-djpeg',  filename);
        FCprintf('Wrote %s\n', filename);
    end
    
 
    %done with this plot
    
    drawnow
    return
    