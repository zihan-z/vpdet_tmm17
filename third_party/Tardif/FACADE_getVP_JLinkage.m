%
% Get vanishing point using JLinkage
%
% ARGS.JL_GRICselect: experimental selection of vp based on information criteria
%
function [vsVP,vCluster] = FACADE_getVP_JLinkage(vsEdges, im, ARGS)

    addpath('../../CPP/mexFunctions/');

    const = FACADE_const();
        
    %set threshold
    switch ARGS.ERROR
      case const.ERROR_DIST
       THRESH = ARGS.JL_ERROR_DIST / ARGS.imgS; %overide ARGS.ERROR_DIST
     case const.ERROR_NDIST
        THRESH = ARGS.ERROR_NDIST*2;
      case const.ERROR_ANGLE
        THRESH = ARGS.ERROR_ANGLE*2;
    end

    %get hypothesis matrix
    if ~ARGS.ALL_calibrated
        FCprintf('Not calibrated| SOLVE_VP: %d\n', ARGS.JL_solveVP )
    elseif ARGS.ALL_calibrated                        
        FCprintf('Calibrated | SOLVE_VP: %d\n', ARGS.JL_solveVP );
    end
    totd =genHyp(ARGS,THRESH,vsEdges, ARGS.imgS, im, ARGS.JL_minNModel*10, ARGS.JL_minNModel);

    % Perform J-Linkage clusterization        
    [T, totdbin] = clusterPoints(totd, THRESH, ARGS.JL_ALGO);
    vCluster = T(:)';         
    FCprintf('nb Class %d\n', max(vCluster) );
                
    %fit VP to each cluster
    for c=1:max(vCluster)
        vb = (vCluster==c);% & ~vbOutliers;
        vsEdges__ = [vsEdges(vb)];
        
        switch ARGS.JL_solveVP
          case const.SOLVE_VP_MLE 
            vsVP(c).VP=FACADE_mxFitVP_x_MLE( vsEdges__, [vsEdges__.nb]);
          case const.SOLVE_VP_MAX
            vsVP(c).VP=FACADE_mxFitVP_x_MAX( vsEdges__, [vsEdges__.nb] );
          case const.SOLVE_VP_GS
            vsVP(c).VP=FACADE_mxFitVP_x_GS(  vsEdges__, [vsEdges__.nb] );
            %vsVP(c).VP=FACADE_mxFitVP_x_GS(  vsEdges__ );
        end
    end


    %---------------------------------------------------
    if ARGS.JL_GRICselect
        %GRIC     
        vModel = getGricData(ARGS, vsEdges, ARGS.imgS, vsVP, vCluster);
        %building GRIC matrix
        R = ARGS.GRIC_R; %dimensions of data
        
        P = vsEdges(1).Pout;      %identical for all
        Nt = length(vsEdges);            %page 4 bottom-right   
        [mQ,vBmodel] =mxBuildGRIC_Q(vModel, Nt, P, R);
        vModel = [vModel(vBmodel)];
        mQ(isnan(mQ)) = -1; 
        mQ(isinf(mQ)) = -1;
        
        if false 
            FCprintf('removed %d model\n', sum(~vBmodel));
            %diag(mQ)<0
            try 
                assert(false);
                [bestval, bestind] = mxSubmod_search_max(mQ);
            catch
                %using matlab file, 
                %warning('mxSubmod_search_max.cpp should be compiled\n');
                [bestval, bestind] = submod_search_max(mQ);
            end
            
            bestval
            mQ
            mQ(bestind,bestind)                
        end
        
        %keep ones with positive, otherwise, diminish entropy of system
        vQdiag = diag(mQ); vQdiag(:)'
        [vsVP, vCluster] = FACADE_sortClass(vsVP, vCluster, vQdiag);
        nbVp = sum(vQdiag>0);
        FCprintf('GRIC criterion: %d -> %d\n', length(vsVP), nbVp);
        vsVP = vsVP(1:nbVp);
        vCluster(vCluster>nbVp) = -1; %eliminate those line
        
        %DEBUG: check ordering
        %vModel = getGricData(ARGS, vsEdges, ARGS.imgS, vsVP, vCluster);        
        %[mQ,vBmodel] =mxBuildGRIC_Q(vModel, Nt, P, R);
        %assert(all(diag(mQ)>0));
    end
    
    
    %---------------------------------------------------
    %order VP 
    if ARGS.ALL_calibrated
        %need to use both  in this order
        [vsVP, vCluster] = FACADE_orderVP_nbPts(vsVP, vsEdges, vCluster);
        %look for manhattan frame
        [vsVP, vCluster] = FACADE_orderVP_Mahattan(vsVP, vsEdges, vCluster);
    else
        %ordering based on nb of edges corresponding to vp
        [vsVP, vCluster] = FACADE_orderVP_nbPts(vsVP, vsEdges, vCluster);
    end
    

if false %ARGS.plot
    sfigure(103);clf;
    res = showModelsPreferenceSet(totdbin, T);
    spy(res, 'b', 1);
    FACADE_printGraphics('', 'Model hypotheses', 'Edges',{}, 'fig/totdbin_clust.eps')
end

return

%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
function vModel = getGricData(ARGS, vsEdges, imgS, vsVP, vCluster)
    
    nbModel =0;
    for i=1:length(vsVP)
        
        %residuals
        vErr =FACADE_mxFitL_x_constVP(vsEdges,...
                                      vsVP(i).VP) * imgS;            
        %row vectors
        vErr = (vErr(:))';
        %vMaskIn = (vMaskIn(:))';
        %vMaskIn = vCluster==i & vErr <= ARGS.JL_pixelThreshold;
        vMaskIn = vErr <= ARGS.JL_pixelThreshold;
        vMaskIn = vMaskIn';
        sigma = sqrt(sum(vErr(vMaskIn).^2)/(sum(vMaskIn)-2));

        %vMaskInErr  = vErr <= ARGS.GRIC_pixelThreshold;
        %vMaskIn = vMaskInErr;
        %sigmaErr = sqrt(sum(vErr(vMaskInErr).^2)/(sum(vMaskInErr)-2));
        %disp([sum(vMaskIn), sigma, sum(vMaskInErr), sigmaErr])
        %[sigma,vMaskIn,b1,b2] = tsse(vErr,2); 
        %sigma
        %sfigure(1000)
        %plot(sort(vErr))
        
        %compute error        
        nbModel = nbModel+1;
        vModel(nbModel).model    = vsVP.VP;
        vModel(nbModel).sigma    = sigma;
        vModel(nbModel).sigmaSq  = sigma^2;
        vModel(nbModel).K        = ARGS.GRIC_K;
        vModel(nbModel).D        = ARGS.GRIC_D;
        vModel(nbModel).N        = sum(vMaskIn);
        vModel(nbModel).vResidSq = vErr .^2;
        vModel(nbModel).E        = sum(vModel(nbModel).vResidSq(vMaskIn))/ sigma^2;  %eq 15
        vModel(nbModel).V        = vMaskIn;
        %likelihood of each point
        vModel(nbModel).vLpV     = exp(-vModel(nbModel).vResidSq/vModel(nbModel).sigmaSq) / (sigma*sqrt(2*pi)); %eq 8
    end

    return
    
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
function [totd] =genHyp(ARGS, THRESH, vsEdges,  imgS, im, nMaxTrials, minNSol)
    
    const = FACADE_const();
    sampleSize = ARGS.ALL_sampleSize;
    
    %Pre-generate a big bunch of numbers, because its faster    
    vNUSamples = FACADE_getSamples(vsEdges, ARGS.ALL_samplingAlgo);
    vIdxAll = randInt(1,length(vNUSamples),1,2*nMaxTrials*sampleSize);
    vIdxAll = vNUSamples(vIdxAll);

    vL = [vsEdges.halfLength]';
    vL2 = vL*2;
    %assert(all(vS<1));
  
    %------------------------------
    ii=1;
    for i = 1:nMaxTrials
        
        vIdx = vIdxAll(3*i:3*i+sampleSize-1);
        vsEdges__ = [vsEdges(vIdx)];        
        %find vanishing point
        vsVP_VP = FACADE_mxFitVP_x_MAX(vsEdges__);

        %error
        vErr =FACADE_mxFitL_x_constVP(vsEdges,...
                                      vsVP_VP);

        %refine solution with inliers
        %-->not recommended/necessary
        
        switch ARGS.ERROR
          case const.ERROR_DIST
            %nothing
          case const.ERROR_NDIST
            vErr = vErr ./ vL2;
          case const.ERROR_ANGLE
            vErr = atan2(vErr,vL)*180/pi;
        end

        %min number of inlier
        if sum(vErr <= THRESH ) <  ARGS.ALL_minInliers 
            continue;
        end

        %save that for later
        cErr{ii} = vErr;
        ii=ii+1;
        if ii==minNSol
            break;
        end
    end
    
    %join all cells with vectors
    totd = [cErr{:}];
    
    FCprintf('done (%d trial, %d sol)\n', i,ii);
    
    
    return
    
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------

    
 
% $$$ function doPlot(im, vsEdges, vsEdges_a, vsEdges_b, vsEdges_c)
% $$$     
% $$$     sfigure(2);clf(2);
% $$$     %subplot(1,2,1);
% $$$     imshow(im);
% $$$     hold on
% $$$     vPts =  [vsEdges.vPts_un];    
% $$$     plot(vPts(1,:), vPts(2,:),'.', 'Color',  [0,0,0] );
% $$$ 
% $$$     vPts =  [vsEdges_c.vPts_un];    
% $$$     plot(vPts(1,:), vPts(2,:),'.', 'Color',  [0,0,1] );
% $$$ 
% $$$     
% $$$     vPts =  [vsEdges_b.vPts_un];    
% $$$     plot(vPts(1,:), vPts(2,:), 'Color',  [0,1,0] );
% $$$     
% $$$     vPts =  [vsEdges_a.vPts_un];    
% $$$     plot(vPts(1,:), vPts(2,:), 'Color',  [1,0,0] );    
% $$$         
% $$$     
% $$$  
% $$$     
% $$$     hold off
% $$$ 
% $$$     
% $$$         %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ function [totd] =genHyp_MLE(ARGS, vsEdges, imgS, im, nMaxTrials, minNSol)
% $$$     
% $$$     sampleSize = ARGS.JL_sampleSize;
% $$$     
% $$$     %Pre-generate a big bunch of numbers, because its faster    
% $$$     vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vsEdges))+1;
% $$$         
% $$$     ii=1;
% $$$     for i = 1:nMaxTrials
% $$$         
% $$$         vIdx = vIdxAll(3*i:3*i+sampleSize-1);
% $$$         vsEdges__ = [vsEdges(vIdx)];        
% $$$         %find vanishing point
% $$$ 
% $$$         
% $$$         vsVP_VP = FACADE_mxFitVP_x_MAX(vsEdges__);
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges__,...
% $$$                                       vsVP_VP);
% $$$         if max(vErr*imgS)>ARGS.JL_pixelThreshold 
% $$$             continue;
% $$$         end
% $$$         
% $$$         %maximal error per line
% $$$         vErr{ii} = FACADE_mxFitL_x_constVP_MLE(vsEdges,...
% $$$                                                 vsVP_VP);
% $$$         if sum([vErr{ii}]*imgS < ARGS.JL_pixelThreshold ) < 5
% $$$             continue;
% $$$         end
% $$$         
% $$$         ii=ii+1;
% $$$         if ii==minNSol
% $$$             break;
% $$$         end
% $$$     end
% $$$     totd = [vErr{:}];
% $$$     
% $$$     FCprintf('done (%d trial, %d sol)\n', i,ii);
% $$$     
% $$$     totd = totd*imgS; %de-normalize
% $$$     
% $$$     
% $$$     return
% $$$     
% $$$ %-----------------------------------------------------------------
% $$$ %----------------- -Error in Gaussian sphere----------------------
% $$$ %-----------------------------------------------------------------
% $$$ function [totd] =genHyp_GS(ARGS, vsEdges, imgS, im, nMaxTrials, minNSol)
% $$$     
% $$$     sampleSize = ARGS.JL_sampleSize;
% $$$ 
% $$$     
% $$$     %Pre-generate a big bunch of numbers, because its faster    
% $$$     vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vsEdges))+1;
% $$$         
% $$$     vL = [vsEdges.vL];
% $$$     assert(size(vL,1)==3);
% $$$     
% $$$     ii=1;
% $$$     for i = 1:nMaxTrials
% $$$         
% $$$         vIdx = vIdxAll(3*i:3*i+sampleSize-1);
% $$$         vsEdges__ = [vsEdges(vIdx)];        
% $$$         %find vanishing point
% $$$         %VP = FACADE_mxFitVP_x_MAX(vsEdges__)
% $$$         VP=FACADE_mxFitVP_x_noConst(vsEdges__);
% $$$         %pause
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges__,...
% $$$                                       VP);
% $$$         if max(vErr*imgS)>ARGS.JL_pixelThreshold 
% $$$             continue;
% $$$         end
% $$$         %abs((vsVP_VP'* vL(:,vIdx)))
% $$$         %pause
% $$$         
% $$$         %VP = VP/norm(VP);
% $$$         %maximal error per line
% $$$         vErr{ii} = abs((VP'* vL)');
% $$$         
% $$$         ii=ii+1;
% $$$         if ii>=minNSol
% $$$             break;
% $$$         end
% $$$     end
% $$$     totd = [vErr{:}];
% $$$     
% $$$     FCprintf('done (%d trial, %d sol)\n', i,ii);
% $$$     
% $$$        
% $$$     
% $$$     return
% $$$     
% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ function [totd] =genHypCalib(ARGS,vsEdges, imgS, im, nMaxTrials, minNSol)
% $$$ 
% $$$     mW = eye(3,3);
% $$$     
% $$$     sampleSize = 4;
% $$$     
% $$$     %Pre-generate a big bunch of numbers, because its faster    
% $$$     vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vsEdges))+1;
% $$$         
% $$$     vErr = {};
% $$$     ii=1;
% $$$     for i = 1:nMaxTrials
% $$$         
% $$$         vIdx = vIdxAll(3*i:3*i+sampleSize-1);
% $$$         
% $$$         vsEdges__a = [vsEdges(vIdx(1:2))];
% $$$         vsEdges__b = [vsEdges(vIdx(3:4))];
% $$$         
% $$$         %find vanishing point
% $$$         %vsVPa=FACADE_fitVP_x_MAX(vsEdges__a, vW(1:2));
% $$$         %vsVPb=FACADE_fitVP_x_MAX(vsEdges__b, vW(1:2));
% $$$         vsVPa_VP = FACADE_mxFitVP_x_MAX(vsEdges__a);
% $$$         vsVPb_VP = FACADE_mxFitVP_x_MAX(vsEdges__b);
% $$$         
% $$$         %useless if sampleSize is 4= 2 per line
% $$$ % $$$         if false
% $$$ % $$$             vErrA =FACADE_mxFitL_x_constVP(vsEdges__a,...
% $$$ % $$$                                            vsVP_VPa);
% $$$ % $$$             vErrB =FACADE_mxFitL_x_constVP(vsEdges__b,...
% $$$ % $$$                                            vsVP_VPb);
% $$$ % $$$             
% $$$ % $$$             if max(vErrA*imgS)>2 | max(vErrB*imgS)>2
% $$$ % $$$                 continue;
% $$$ % $$$             end
% $$$ % $$$         end
% $$$ % $$$         
% $$$         %project 2nd vp on space perpendicular to 1st
% $$$         vLb = mW*vsVPa_VP;
% $$$         %vLb = vLb / sqrt(vLb(1)^2+vLb(2)^2);
% $$$         vp_onvLb = FACADE_closest_x_on_l(vLb, vsVPb_VP);
% $$$         if true
% $$$             vErr =FACADE_mxFitL_x_constVP(vsEdges__b,...
% $$$                                           vp_onvLb);
% $$$             if max(vErr*imgS)>2
% $$$                 continue;
% $$$             end
% $$$         else
% $$$             vErr =FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                           vp_onvLb);
% $$$             if sum(vErr*imgS<=2)<5
% $$$                 continue;
% $$$             end
% $$$             vErr =FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                           vsVPa_VP);
% $$$             if sum(vErr*imgS<=2)<5
% $$$                 continue;
% $$$             end
% $$$         end
% $$$         
% $$$         %maximal error per line
% $$$         vErr{ii} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                              vsVPa_VP);
% $$$         ii=ii+1;
% $$$         %vErr{ii+1} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$         %                                        vsVPb_VP);
% $$$         %ii=ii+1;
% $$$         if ii>=minNSol
% $$$             break;
% $$$         end
% $$$     end
% $$$     totd = [vErr{:}];
% $$$ 
% $$$     if isempty(vErr)
% $$$         warning('failed to get hypothesis, falling on uncalibrated\n');        
% $$$         totd =genHyp(ARGS, vsEdges, imgS, im, nMaxTrials, minNSol);
% $$$         return;
% $$$     end
% $$$     
% $$$     FCprintf('done (%d trial, %d sol)\n', i,ii);
% $$$     
% $$$     totd = totd*imgS; %de-normalize
% $$$     
% $$$     
% $$$     return
% $$$     
% $$$     
% $$$     %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ function [totd] =genHypCalib2(ARGS,vsEdges, imgS, im, nMaxTrials, minNSol)
% $$$ 
% $$$     mW = eye(3,3);
% $$$     
% $$$     sampleSize = 4;
% $$$     
% $$$     %Pre-generate a big bunch of numbers, because its faster    
% $$$     vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vsEdges))+1;
% $$$         
% $$$     vErr = {};
% $$$     ii=1;
% $$$     for i = 1:nMaxTrials
% $$$         
% $$$         vIdx = vIdxAll(3*i:3*i+sampleSize-1);
% $$$         
% $$$         vsEdges__a = [vsEdges(vIdx(1:2))];
% $$$         vsEdges__b = [vsEdges(vIdx(3:4))];
% $$$         
% $$$         %find vanishing point
% $$$         vsVPa_VP = FACADE_mxFitVP_x_MAX(vsEdges__a);
% $$$         
% $$$         %project 2nd vp on space perpendicular to 1st
% $$$         vLb = mW*vsVPa_VP;
% $$$         vp_onvLb = vsEdges__b(1).skmvL* vLb;
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges__b,...
% $$$                                       vp_onvLb);
% $$$         if max(vErr*imgS)>2
% $$$             continue;
% $$$         end
% $$$         vp_onvLb = skew_mat(vsEdges__b(2).vL) * vLb;
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges__b,...
% $$$                                       vp_onvLb);
% $$$         if max(vErr*imgS)>2
% $$$             continue;
% $$$         end
% $$$         
% $$$         %maximal error per line
% $$$         vErr{ii} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                              vsVPa_VP);
% $$$         ii=ii+1;
% $$$         %vErr{ii+1} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$         %                                        vsVPb_VP);
% $$$         %ii=ii+1;
% $$$         if ii>=minNSol
% $$$             break;
% $$$         end
% $$$     end
% $$$     totd = [vErr{:}];
% $$$ 
% $$$     if isempty(vErr)
% $$$         warning('failed to get hypothesis, falling on uncalibrated\n');        
% $$$         totd =genHyp(ARGS, vsEdges, imgS, im, nMaxTrials, minNSol);
% $$$         return;
% $$$     end
% $$$     
% $$$     FCprintf('done (%d trial, %d sol)\n', i,ii);
% $$$     
% $$$     totd = totd*imgS; %de-normalize
% $$$     
% $$$     
% $$$     return
% $$$     


% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ %-----------------------------------------------------------------
% $$$ function [totd] =genHypCalib3(ARGS,vsEdges, vNUSamples, imgS, im, nMaxTrials, minNSol)
% $$$ 
% $$$     mW = eye(3,3);
% $$$     
% $$$     sampleSize = 4;
% $$$     
% $$$     %Pre-generate a big bunch of numbers, because its faster    
% $$$     if true
% $$$         %per points
% $$$         vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vNUSamples))+1;
% $$$         vIdxAll = vNUSamples(vIdxAll);
% $$$     else
% $$$         %per line
% $$$         vIdxAll = randInt(1,2*nMaxTrials*sampleSize,length(vsEdges))+1;
% $$$     end
% $$$     
% $$$     mL = [vsEdges.vL];
% $$$     
% $$$     vErr = {};
% $$$     ii=1;
% $$$     for i = 1:nMaxTrials
% $$$         
% $$$         vIdx = vIdxAll(3*i:3*i+sampleSize-1);
% $$$         vsEdges__a = [vsEdges(vIdx(1:2))];
% $$$         if false
% $$$             vTheta =[vsEdges__a.theta];
% $$$             diffTheta = min(abs(vTheta(1)-vTheta(2)),...
% $$$                             abs(180+vTheta(1)-vTheta(2)));
% $$$             if diffTheta > 80, continue; end
% $$$             %diffTheta
% $$$         end
% $$$         ARGS.plot= false;
% $$$ 
% $$$         vsEdges__b = [vsEdges(vIdx(3))];                
% $$$         
% $$$         %find vanishing point
% $$$         %vsVPa_VP = FACADE_mxFitVP_x_multiL(vsEdges__a);
% $$$         vsVPa_VP = FACADE_mxFitVP_x_MAX(vsEdges__a);
% $$$ 
% $$$                         
% $$$         %look for a line not consistent with vsVPa_VP
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                       vsVPa_VP);
% $$$         vb = vErr*ARGS.imgS >  3*ARGS.JL_pixelThresholdVPselect;
% $$$         if sum(vb) == 0, continue; end
% $$$         vIdxB =  find(vb);        
% $$$         idxB = randInt(1,1,sum(vb))+1;
% $$$         idxB = vIdxB(idxB);
% $$$          
% $$$         vLb = mW*vsVPa_VP; %vanishing line
% $$$         vp_onvLb = skew_mat(vLb) * mL(:,idxB); %orthogonal vanishing point
% $$$                                                 
% $$$         %look for consisten lines with this new VP
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                       vp_onvLb);
% $$$         vbConst = vErr*imgS<ARGS.JL_pixelThresholdVPselect;
% $$$         vbConst(vIdx(1:2))=false;
% $$$         vbConst(idxB) = false;
% $$$         if sum(vbConst )<3
% $$$             continue;
% $$$         end
% $$$ 
% $$$         if ARGS.plot
% $$$             vsEdges__b = [vsEdges(idxB)];
% $$$             vsEdges__c = [vsEdges(vbConst)];
% $$$             doPlot(im, vsEdges,  vsEdges__a ,  vsEdges__b, vsEdges__c);
% $$$             pause
% $$$         end
% $$$         
% $$$         
% $$$         
% $$$         %maximal error per line
% $$$         vErr{ii} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$                                              vsVPa_VP);
% $$$         ii=ii+1;
% $$$         %vErr{ii+1} = FACADE_mxFitL_x_constVP(vsEdges,...
% $$$         %                                        vsVPb_VP);
% $$$         %ii=ii+1;
% $$$         if ii>=minNSol
% $$$             break;
% $$$         end
% $$$     end
% $$$     totd = [vErr{:}];
% $$$ 
% $$$     if isempty(vErr)
% $$$         warning('failed to get hypothesis, falling on uncalibrated\n');        
% $$$         totd =genHyp(ARGS, vsEdges, imgS, im, nMaxTrials, minNSol);
% $$$         return;
% $$$     end
% $$$     
% $$$     FCprintf('done (%d trial, %d sol)\n', i,ii);
% $$$     
% $$$     totd = totd*imgS; %de-normalize
% $$$     
% $$$     
% $$$     return