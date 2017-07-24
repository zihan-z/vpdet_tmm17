function  [vsVP, vClass] = FACADE_orderVP_Mahattan(vsVP, vsEdges, vClass)
    
    
    vSupport = FACADE_get_ManhattanSupport(vsVP, vsEdges, vClass);
    [vsVP, vClass] = FACADE_sortClass(vsVP, vClass, vSupport);
    %take just 3 points 
    vsVP = vsVP(1:3);
    vb = vClass > 3;
    vClass (vb) = 0;
        
    
    
    
    
    