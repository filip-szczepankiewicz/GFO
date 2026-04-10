% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

clear

nRot = 30;
R = library.loadSet(nRot, 'GFOD2');
B = diag([0 1 2])/2;

%% Plot glyphs
switch 1
    case 1
        B = diag([0 1 2]+0.1);
        rvec = [.8 0.1 0.1];
        enam = 'OTE';

    case 2
        B = diag(sqrt([0.01 0.01 2]));
        rvec = [0 0 1];
        enam = 'LTE';
end

p = fix.subplot('tight');
p.nr = 1;
p.nc = 2;
p.wo = -0.1;

for n = 1:nRot
    Bp(:,:,n) = R(:,:,n)*B*R(:,:,n)';
    xp(n,:) = [1 0 0]*R(:,:,n);
    yp(n,:) = [0 1 0]*R(:,:,n);
    zp(n,:) = [0 0 1]*R(:,:,n);
end
%%
clf

fix.subplot(p, 1)
for i = 1:size(R,3)

    hold on
    Rt = R(:,:,i);

    hp  = fix.plot.superQuadratic(Bp(:,:,i), 'col', abs(rvec*Rt'), 'level', 1, ...
        'eps1', 0.2, 'eps2', 0.5, 'alpha', 0.5);

end

lighting phong;
lightangle(-45,30);
camlight headlight;

campos([1 1 1]*108)
axis off
campos([1 1 1]*15*2)

fix.subplot(p, 2)

hold on

plot_gdirs([xp; -xp], [xp; -xp]*0+fix.color.red)
plot_gdirs([yp; -yp], [xp; -xp]*0+fix.color.green)
plot_gdirs([zp; -zp], [xp; -xp]*0+fix.color.blue)

lighting phong;
material dull;
lightangle(-45,30);
camlight headlight;

campos([1 1 1]*15)
axis off

fix.figure

% fix_save_figure('fig_tensorAndDirDist', '', [6 4], 600)