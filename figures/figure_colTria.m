% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

clear

T1 = @(e) 8*(e(1)-e(2))*(e(2)-e(3))/(sum(e).^2);
T2 = @(e) min(e(1)-e(2), e(2)-e(3)) / (0.5*(e(1)-e(3)) + eps);

clf

cm = flip(fix.cmap.blackredwhite(600)).^1;
cm = cm(1:500,:);

figure(1)
clf

btensor_triangle(T1, 'cmap', cm, 'showCornerLabels', 0)
clim([0 1])

fix.colorbar.position([], [1.15 1 .5 .3], [0 .5 0 0])

fix.plot.axis
fix.figure

% fix.figure.save('tria', '', [6 6], 600)