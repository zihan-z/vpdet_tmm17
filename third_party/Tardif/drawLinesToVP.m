function drawLinesToVP (vsVP, W, H, classClr)

%middlePoint = [W/2; H/2+126; 1];
middlePoint = [W/2; H/2; 1];

for i = 1:size(vsVP, 1)
    try
    point(1:3) = vsVP(i, 1:3)';
    %point = mToH (point);
    
    l = cross (middlePoint, point);
    if (abs (point(1) - middlePoint(1)) > abs(point(2) - middlePoint(2)))
        x = min (middlePoint (1), point(1)) : max (middlePoint (1), point(1));
        y = (-l(3) -l(1) * x)/l(2);
    else
        y = min (middlePoint (2), point(2)) : max (middlePoint (2), point(2));
        x = (-l(3) -l(2) * y)/l(1);
    end
    hold on
    
    
    
    plot( x, y, 'Color',  classClr(i, :), 'LineWidth', 4 );
    catch
        fprintf ('smth wrong in drawing lines');
    end
end


end

