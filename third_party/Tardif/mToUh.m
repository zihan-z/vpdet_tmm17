%put back un euclidean space
function mM=mToUh(mM)
    for i=1:size(mM,1)-1
        mM(i,:) = mM(i,:)./mM(end,:);
    end
    mM = mM(1:2,:);
    