function rScalar=randInt(min,max, row,col);

if (nargin <2)
    error('[randInt] you must specify min and max');
elseif nargin <=2
    row=1;
    col=1;
elseif nargin<=3
    col = row;
end;


rScalar =round((max-min)*rand(row,col)+min);

return
