% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

clear

w = 0.25;
xl = linspace(0, 1-w, 7);
yl = linspace(0, 1-w, 4);

l1 = 1;
l2 = [0 .3 .3 .6 .6 .6 .9 .9 .9 .93]+0.1;
l3 = [0  0 .3 0 .2 .6 0 .3 .7 .9]+0.1;

l1 = 1;
l2 = [0.02 .27 .27 .6 .6  .6 .94 .94 .94 .94]+0.06;
l3 = [0.02  0  .27  0 .27 .6 0.04 .27 .6 .94]+0.06;

% e1 = 1;
% e2 = 0.7;
% el = [e1 e2 e2 e2 e2 e2 e1 e2 e2 e1];

e = [ones(10,1) l2' l3'];

% e = e./sum(e,2)

el = 1-8*(e(:,1)-e(:,2)).*(e(:,2)-e(:,3))./(sum(e,2).^2);

el = (el).^1.5;

x = xl([4 3 5 2 4 6 1 3 5 7]);
y = yl([4 3 3 2 2 2 1 1 1 1]);
rvec = [1 1 1]*0.7;
clf


for i = 1:numel(l2)

    B = diag([l1 l2(i) l3(i)]);

    axes('position', [x(i) y(i) w w])
    fix.plot.superQuadratic(B, 'col', abs(rvec).^2, 'level', 1, ...
        'eps1', el(i), 'eps2', el(i));

    title ''
    axis off

    axis([-1 1 -1 1 -1 1])
    lighting phong;
    lightangle(-15,30);
    camlight;

end

fix.figure

set(gcf, 'color', 'k')

fix.figure.save('tensor_triangle', '', [6.2 5.3], 600)
