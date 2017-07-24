    
%normalize each column of m
function m=vToNorm(m)
    
    v = vNorm(m);
    
    for i=1:size(m,1)
        m(i,:) = m(i,:)  ./ v;
    end
    
    return