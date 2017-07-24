%----------------------------------------------------
%----------------------------------------------------
%----------------------------------------------------
function vsEdges = FACADE_build_vsEdges(vPts_un, ARGS, noFit)
    const = FACADE_const();
    
    nbLine = length(vPts_un);
    
    %-------------------------------------------------
    %area of the image
    if isempty(ARGS.mKinv)
        ARGS.mKinv = diag([1/ARGS.imgS, 1/ARGS.imgS, 1]);
        ARGS.mK    = diag([  ARGS.imgS,   ARGS.imgS, 1]);
    end
    vCrn = [0 ARGS.imgW 0          ARGS.imgW;...
            0 0         ARGS.imgH  ARGS.imgH;...
            1 1         1          1];
    vCrnN = mToUh(ARGS.mKinv*vCrn);
    imgAreaN = abs(vCrnN(1,2) - vCrnN(1,1))* abs(vCrnN(2,2) - vCrnN(2,4));
    imgArea  = abs(vCrn(1,2)  - vCrn(1,1))*  abs(vCrn(2,2)  - vCrn(2,4));
    imgDiag  = norm(vCrn(1:2,1)  - vCrn(1:2,4),2);
    imgDiagN = norm(vCrnN(1:2,1) - vCrnN(1:2,4),2);
    
    %-------------------------------------------------
    
    for j=1:nbLine 
        
        vsEdges(j).vPts_un = vPts_un{j};
        vPts = vPts_un{j};
        
        if isempty(ARGS.mKinv)
            %normalization
            vPts(1,:) = (vPts(1,:) - ARGS.imgW/2)/ARGS.imgS;
            vPts(2,:) = (vPts(2,:) - ARGS.imgH/2)/ARGS.imgS;
            ARGS.mKinv = diag([1/ARGS.imgS, 1/ARGS.imgS, 1]);
        else
            vPts = mToUh(ARGS.mKinv*mToH(vPts));
        end
        
               
        %line fitting and subpixel accuracy
        vL = FACADE_fitL_x(vPts);  vL = vL/norm(vL(1:2));       
        if ~noFit
            vPts =  FACADE_closest_x_on_l(vL, vPts);
            %update unnormalized values
            if isempty(ARGS.mKinv)
                vsEdges(j).vPts_un(1,:) = vPts(1,:)*ARGS.imgS + ARGS.imgW/2;
                vsEdges(j).vPts_un(2,:) = vPts(2,:)*ARGS.imgS + ARGS.imgH/2;
            else
                vsEdges(j).vPts_un = mToUh(ARGS.mK*mToH(vPts));
            end
        end
        %invert Y
        vPts(2,:) = -vPts(2,:);
        vsEdges(j).vPts = vPts;
        vL = FACADE_fitL_x(vPts);
        vsEdges(j).vL = vL/norm(vL);
        vsEdges(j).vLtL = vL*vL';
        
        
        vsEdges(j).nb   = size(vPts,2);
        vsEdges(j).vPointUn1 = vsEdges(j).vPts_un(:,1);
        vsEdges(j).vPointUn2 = vsEdges(j).vPts_un(:,end);
        
        vsEdges(j).vPoint1 = vsEdges(j).vPts(:,1);
        vsEdges(j).vPoint2 = vsEdges(j).vPts(:,end);
        dir = vsEdges(j).vPoint1 - vsEdges(j).vPoint2;
        dir = dir/norm(dir);
        vsEdges(j).theta = atan2(dir(2),dir(1))*180/pi;
        %for FACADE_mxFiVP_x_MAX
        vsEdges(j).vCtE = (vsEdges(j).vPoint1  + vsEdges(j).vPoint2)/2;
        vsEdges(j).skmVctE = skew_mat([vsEdges(j).vCtE;1]);
        vsEdges(j).vL2 = vsEdges(j).skmVctE * [vsEdges(j).vPoint1;1];

        vsEdges(j).halfLength =  norm(vPts(1:2,1) - vPts(1:2,end),2)/2;

        
        %usefull stuff to precompute
        vsEdges(j).vCt     = mean(vPts')';       %centroid
        vsEdges(j).vCt_un     = mean(vsEdges(j).vPts_un')';       %centroid
        vsEdges(j).vCth    = [vsEdges(j).vCt;1]; %centroid homogenous coordinate
        v1 = ones(size(vPts,2),1);        
        vsEdges(j).mX = [vPts', v1];
        vsEdges(j).mXtmX = (vsEdges(j).mX)'*vsEdges(j).mX ;
        if (numel (vsEdges(j).vCth) == 2)
            fprintf ('');
        end
        vsEdges(j).skmVct = skew_mat(vsEdges(j).vCth);
        msk = vsEdges(j).skmVct;
        nbPts = size(vPts,2);
        vsEdges(j).mAtA = msk'*vsEdges(j).mXtmX*msk;% * nbPts^2;
        

        
        %set threshold
        switch ARGS.ERROR
          case const.ERROR_DIST
            vsEdges(j).Pout  = 1/imgArea^2 ;
          case const.ERROR_NDIST
            vsEdges(j).Pout  = 1/imgArea^2 ;
          case const.ERROR_ANGLE
            %warning('Pout is incorrect');
            vsEdges(j).Pout  = 1/imgArea^2 ;
        end
        
        %Probability of being an outlier is inv prop to lenght
        %  is this identical for all?, I think so
        %1) parametrization account for line length and position
        %vsEdges(j).Pout  = 1/imgArea^2 ;
        
        %2) for lines
        %vsEdges(j).Pout  = 1/ (4*pi*1^2); %surface of sphere
                                          
        %3) for lines and length
        %vsEdges(j).Pout  = 1/(4*pi*1^2*imgDiag); %volume of sphere of radius longhest possible line
        
        
        

    end
    
    