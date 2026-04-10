% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

clear

w = 0.25;
xl = linspace(0, 1-w, 7);
yl = linspace(0, 1-w, 4);

l1 = 1;
l2 = [0 .3 .3 .6 .6 .6 .9 .9 .9 .9]+0.1;
l3 = [0  0 .3 0 .2 .6 0 .3 .7 .9]+0.1;

e1 = 1;
e2 = 0.7;
el = [e1 e2 e2 e2 e2 e2 e1 e2 e2 e1];

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

% fix.figure.save('tensor_triangle', '', [5 4], 600)
