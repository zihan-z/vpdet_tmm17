function [vsEdges,imgEdge] = FACADE_getEdgelets2(img, ARGS)


addpath ('../../CPP/mexFunctions/');
addpath ('../../CPP/mxGetEdges/');

%ARGS.savePlot = false;
    pFigEdges = 100;
    iptsetpref('ImshowBorder', 'tight');
    
    %--------------------------------------------------
    FCprintf('find connected components\n');
    %[edgelist, labelededgeim] = edgelink(imE, 10, imgCplxEdge);
    [edgelist, imgEdge] = mxGetEdges(uint8(img), ARGS.minEdgeLength);
    %edgelist = cleanedgelist(edgelist, ARGS.minEdgeLength);
    
    %--------------------------------------------------
    FCprintf('segment connected components in straight lines\n');
    edgelist = lineseg2(edgelist, ARGS.linesegTOL);
    vn = cellfun(@length, edgelist);
    edgelist = {edgelist{vn>ARGS.minEdgeLength}};
    FCprintf('Found %d lines\n', length(edgelist));
    
    if ARGS.plot
        sfigure(pFigEdges); clf(pFigEdges);
        set(pFigEdges, 'PaperPositionMode', 'auto')
        set(pFigEdges, 'inverthardcopy',    'off');
        imshow(img, 'InitialMagnification', 100)
        %imshow(imgEdge, 'InitialMag', 100)
        hold on
        drawedgelist(edgelist, size(img), 3, 'rand'); %axis off        
        hold off
        set(pFigEdges, 'PaperPositionMode', 'auto')
        set(pFigEdges, 'inverthardcopy',    'off');
        
        drawnow 
        if ARGS.savePlot
            filename= sprintf( 'fig/__tmp_edge_overlay.jpg');
            %saveas(pFigEdges,filename, 'jpg'); 
            print('-noui' ,'-r0', '-djpeg',  filename);
            %print('-r0', '-djpeg',  filename);
            FCprintf('Wrote %s\n', filename);
        end
    end
    %-----------------------------------------------------------------
    ARGS.imgW =size(img,2);
    ARGS.imgH =size(img,1);
               
    nbLine =  length(edgelist);
       
    %-----------------------------------------------------------------
    FCprintf('normalization and subpixel accuracy by line fitting\n');
    
    m=[[0,1];...
       [1,0]];
    for j=1:nbLine
        vPts_un{j} =m*[edgelist{j}]'+1;
    end
    vsEdges = FACADE_build_vsEdges(vPts_un, ARGS, false);
    
    return
    
