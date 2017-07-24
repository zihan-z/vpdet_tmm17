function b = isInList(vStr, str)
    
    b=false;
    for i=1:length(vStr)
        
        if strcmp(vStr{i},str)
            b=true; 
            return;
        end
    end