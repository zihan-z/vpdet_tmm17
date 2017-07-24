function h = invfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
% Modified by Jean-Philippe Tardif
%
% See also figure

if nargin>=1 
    if ishandle(h)
        set(0, 'CurrentFigure', h);
        set(h,'Visible','off'); 
    else
        if isempty(h)
            h=sfigure();
        else
            h = figure(h);
            set(h,'Visible','off') 
        end
    end
else
    h = figure('Visible','off');
end
