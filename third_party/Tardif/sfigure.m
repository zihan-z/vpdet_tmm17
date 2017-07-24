function h = sfigure(h, winTitle)
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
    else
        if isempty(h)
            h=sfigure();
        else
            h = figure(h);
        end
    end
else
    h = figure;
end

if nargin >= 2
    set(gcf,'Name', winTitle);
end