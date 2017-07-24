function vbOutliers = FACADE_getOutliers(ARGS,vsEdges, vClass, vsVP)

const = FACADE_const();

%set threshold
 switch ARGS.ERROR
  case const.ERROR_DIST
   THRESH = ARGS.JL_ERROR_DIST / ARGS.imgS; %overide ARGS.ERROR_DIST
  case const.ERROR_NDIST
   THRESH = ARGS.ERROR_NDIST;
  case const.ERROR_ANGLE
   THRESH = ARGS.ERROR_ANGLE;
 end
    
    
vbOutliers = logical(zeros(1,length(vsEdges)));

vL = [vsEdges.halfLength]';
vL2 = vL*2;

for i =1:length(vsVP)
  vb = (vClass==i);
  vErr = FACADE_mxFitL_x_constVP([vsEdges(vb)],...
                                 vsVP(i).VP);
  
  switch ARGS.ERROR
   case const.ERROR_DIST
    %nothing
   case const.ERROR_NDIST
    vErr = vErr ./ vL2(vb);
   case const.ERROR_ANGLE
    vErr = atan2(vErr,vL)*180/pi;
  end
  
  vbOutliers(vb) = vErr>THRESH;
end
    
