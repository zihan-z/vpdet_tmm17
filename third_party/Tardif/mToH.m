%in homogenious coord
function mM=mToH(mM)
    v1 = ones(1,size(mM,2));
    mM = [mM;v1];
    return;