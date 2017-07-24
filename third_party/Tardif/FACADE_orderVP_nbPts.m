function [vsVP, vClass] = FACADE_orderVP_nbPts(vsVP, vsEdges, vClass)
    
    
    nbPtTotal = sum([vsEdges.nb]);
    for c =1:max(vClass)
        vb = (vClass==c);
        vSupport(c) = sum(vb);%sum([vsEdges(vb).nb])/nbPtTotal;
                              %vSupport(c) = sum([vsEdges(vb).nb])/nbPtTotal;
    end 
    
    [vsVP, vClass] = FACADE_sortClass(vsVP, vClass, vSupport);