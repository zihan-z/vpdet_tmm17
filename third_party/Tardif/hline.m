% HLINE - Plot 2D lines defined in homogeneous coordinates.
%
% Function for ploting 2D homogeneous lines defined by 2 points
% or a line defined by a single homogeneous vector
%
% Usage:   hline(p1,p2)   where p1 and p2 are 2D homogeneous points.
%          hline(p1,p2,'colour_name')  'black' 'red' 'white' etc
%          hline(l)       where l is a line in homogeneous coordinates
%          hline(l,'colour_name')
%

%  Peter Kovesi
%  School of Computer Science & Software Engineering
%  The University of Western Australia
%  pk @ csse uwa edu au
%  http://www.csse.uwa.edu.au/~pk
%
%  April 2000

function hline(a,varargin)

a = a./a(3);   % ensure line in z = 1 plane (not needed??)

if abs(a(1)) > abs(a(2))   % line is more vertical
    
    ylim = get(get(gcf,'CurrentAxes'),'Ylim');
    p1 = hcross(a, [0 -1/ylim(1) 1]');      
    p2 = hcross(a, [0 -1/ylim(2) 1]');
    %p2 = hcross(a, [0 1/SC -1]');
else                       % line more horizontal
    xlim = get(get(gcf,'CurrentAxes'),'Xlim')   ;
    
    p1 = hcross(a, [-1/xlim(1) 0 1]');
    p2 = hcross(a, [-1/xlim(2) 0 1]');
end

line([p1(1) p2(1)], [p1(2) p2(2)], varargin{:});


% HCROSS - Homogeneous cross product, result normalised to s = 1.
%
% Function to form cross product between two points, or lines,
% in homogeneous coodinates.  The result is normalised to lie
% in the scale = 1 plane.
% 
% Usage: c = hcross(a,b)
%

% Copyright (c) 2000-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

%  April 2000

function c = hcross(a,b)
c = cross(a,b);
c = c/c(3);
