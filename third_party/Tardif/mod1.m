%modulo between 1 and md
%as opposed to between 0 and md-1, as usual

function v=mod1(v, md)
    
    
    v= mod(v -1, md) +1;
