%normalize each column of m by the corresponding component of v
function m=mvScale(m, v)
    
    for i=1:size(m,1)
        m(i,:) = m(i,:)  .* v;
    end
    
    return