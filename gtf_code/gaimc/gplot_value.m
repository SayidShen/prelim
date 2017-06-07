function [ Xout, Yout ] = gplot_value( A,xy,v,yrange)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[i,j] = find(A);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';

% if isfloat(xy) || nargout ~= 0
%     X = [X; NaN(size(i))'];
%     Y = [Y; NaN(size(i))'];
% end


maxval=yrange(2);%max(v);
minval=yrange(1);%min(v);
v(v>maxval)=maxval;
v(v<minval)=minval;
colorrange=maxval-minval;
colormap1=colormap('jet');


sz=length(i);
hold off;
for ii=1:sz
    val_here=mean([v(i(ii)),v(j(ii))]);
    cc=colormap1(floor((val_here-minval)/colorrange*63)+1,:);
    plot(X(:,ii),Y(:,ii),'Color',cc,'linewidth',2)
    hold on;
end
hold off;
    Xout = X(:);
    Yout = Y(:);
end

function [lsty, csty, msty] = gplotGetRightLineStyle(ax, lc)
%  gplotGetRightLineStyle
%    Helper function which correctly sets the color, line style, and marker
%    style when plotting the data above.  This style makes sure that the
%    plot is as conformant as possible to gplot from previous versions of
%    MATLAB, even when the coordinates array is not a floating point type.
co = get(ax,'ColorOrder');
lo = get(ax,'LineStyleOrder');
holdstyle = getappdata(gca,'PlotHoldStyle');
if isempty(holdstyle)
    holdstyle = 0;
end
lind = getappdata(gca,'PlotLineStyleIndex');
if isempty(lind) || holdstyle ~= 1
    lind = 1;
end
cind = getappdata(gca,'PlotColorIndex');
if isempty(cind) || holdstyle ~= 1
    cind = 1;
end
nlsty = lo(lind);
ncsty = co(cind,:);
nmsty = 'none';
%  Get the linespec requested by the user.
[lsty,csty,msty] = colstyle(lc);
if isempty(lsty)
    lsty = nlsty;
end
if isempty(csty)
    csty = ncsty;
end
if isempty(msty)
    msty = nmsty;
end

end

