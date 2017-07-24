function vL = FACADE_fitL_x(vp)
    
    if size(vp,1)~=3
        v1 = ones(1,size(vp,2));
        vp = [vp; v1];
    end
    
    mAtA = vp*vp';       
    [v,d]= eig(mAtA);
    [dum, idx] = min(diag(d));
    vL = v(:,idx);