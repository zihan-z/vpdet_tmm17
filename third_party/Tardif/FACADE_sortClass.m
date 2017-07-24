function [vsVP, vClass] = FACADE_sortClass(vsVP, vClass, vSupport)
    
%highest support are first
    
    
    [dum,vIdx] = sort(vSupport, 'descend');
    vsVP = vsVP(vIdx);
    vClass_cpy = vClass;
    for c=1:length(vsVP)
        vClass_cpy(vClass==vIdx(c)) = c;
    end
    vClass =vClass_cpy;
    

    
    
    