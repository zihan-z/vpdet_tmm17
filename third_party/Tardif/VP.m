%stucture for vanishing point

function sVP = VP(vp)
    
    sVP.VP          = vp(:); %make sure its a column
    %sVP.mOrthoSpace = [];

    %build orthospace
    [u,d,v] = svd([sVP.VP';...
                   [0,0,0];...
                   [0,0,0]]);
    sVP.mOrthoSpace = v(:,2:end);