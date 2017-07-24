
    
%norm of each column of matrix m
function v=vNorm(m)
    
    v = sqrt(sum(m .* m));
    
    return
    
