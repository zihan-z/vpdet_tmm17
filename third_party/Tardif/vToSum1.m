%normalize each column of m
function m=vToSum1(m)
    
    v = sum(m);
    
    for i=1:size(m,1)
        m(i,:) = m(i,:)  ./ v;
    end

    m(isnan(m)) = 0;
    
    
    return