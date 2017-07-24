function vSupport = FACADE_get_ManhattanSupport(vsVP, vsEdges, vClass)
    
%classify according to # of edges
%[vsVP, vClass] = FACADE_orderVP_nbPts(vsVP, vsEdges, vClass);
    ARGS.minNbEdges = 3;
    ARGS.forceClassOfMaxNbEdge = true;
    
    nbClass = max(vClass);
    
    if nbClass <=3 
        vSupport = ones(1,length(vsVP));
        return;
    end
    
    
    nbPtTotal = sum([vsEdges.nb]);
    for c =1:nbClass
        vb = (vClass==c);
        vSupportNpts(c) = sum(vb);%sum([vsEdges(vb).nb])/nbPtTotal;
                                  %vSupportNpts(c) = sum([vsEdges(vb).nb])/nbPtTotal;
    end 
    
    mTriplets = pick(1:nbClass,3,'')';
    if ARGS.forceClassOfMaxNbEdge
        %contains first class
        vb = mTriplets(1,:)==1 | ...
             mTriplets(2,:)==1 | ...
             mTriplets(3,:)==1;
        mTriplets = mTriplets(:,vb);
    end
    
    
    for i=1:size(mTriplets,2)
    
        vsVP__    = [vsVP(   mTriplets(:,i))];
        vsEdges_1 = [vsEdges(vClass==mTriplets(1,i))];
        vsEdges_2 = [vsEdges(vClass==mTriplets(2,i))];
        vsEdges_3 = [vsEdges(vClass==mTriplets(3,i))];
        
        if any([length(vsEdges_1), length(vsEdges_2), length(vsEdges_3)] < ARGS.minNbEdges),
            vErr(i) = 1000;
            continue;
        end
        
        vErr(i) = getOrthogonalityError(vsVP__,vsEdges_1, vsEdges_2, vsEdges_3);
    end
    
    %select 3 best vp
    [mx,idxTr] = min(vErr);

% $$$     
% $$$     [mTriplets;...
% $$$      zeros(1,size(mTriplets,2));...
% $$$      vErr]
% $$$     
% $$$      idxTr
% $$$     pause
    
    %our choice
    tr = mTriplets(:,idxTr);    
    vSupport = zeros(1,length(vsVP));    
    vSupport(tr) = vSupportNpts(tr); %support of 3 orthogonal vp is nb of edge
                                     %rest are given zeros
    
    
    
    return
    
function err = getOrthogonalityError(vsVP, vsEdges_1, vsEdges_2, vsEdges_3)
    vL1 = vsVP(1).VP;
    vL2 = vsVP(2).VP;
    vL3 = vsVP(3).VP;
    
    vP1 = null([vL2';vL3']);        
    vP1 = vP1/norm(vP1);
    vP2 = null([vL1';vL3']);        
    vP2 = vP2/norm(vP2);
    vP3 = null([vL1';vL2']);        
    vP3 = vP3/norm(vP3);
    
    vErr1 =FACADE_mxFitL_x_constVP(vsEdges_1, vP1);
    vErr2 =FACADE_mxFitL_x_constVP(vsEdges_2, vP2);
    vErr3 =FACADE_mxFitL_x_constVP(vsEdges_3, vP3);
    
    err = mean([vErr1;vErr2;vErr3]);
    return
% $$$     
% $$$     
% $$$     %---------------------------------------------------
% $$$     %if calibrated, select 3rd vanishing point according to calibration
% $$$     if ARGS.JL_CALIB >=1
% $$$         vL1 = vsVP(1).VP;
% $$$         vL2 = vsVP(2).VP;
% $$$         
% $$$         vP3 = null([vL1';vL2']);        
% $$$         vP3 = vP3/norm(vP3);
% $$$         vErr =FACADE_mxFitL_x_constVP(vsEdges, vP3);
% $$$         vFit(1:2) = 0;
% $$$         for c =3:max(vClass)
% $$$             vb = vClass==c;
% $$$             if sum(vb)<=2
% $$$                 vFit(c) = 1000;
% $$$             else
% $$$                 vFit(c) = mean(vErr(vb));
% $$$             end
% $$$         end
% $$$         
% $$$         %vFit*ARGS.imgS
% $$$         vSupport = -vFit + max(vFit); %reverse of fit
% $$$         [vsVP, vClass] = FACADE_sortClass(vsVP, vClass, vSupport);
% $$$         
% $$$         %check
% $$$         %vErr =FACADE_mxFitL_x_constVP(vsEdges, vP3);
% $$$         %for c =1:max(vClass)
% $$$         %    vFit(c) = mean(vErr(vClass==c));
% $$$         %end
% $$$         %vFit*ARGS.imgS
% $$$     end