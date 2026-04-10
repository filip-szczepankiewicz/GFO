% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

clear

nD = 30;
nR = 1600;

dspace = linspace(0,2,nD);

[d1,d2] = meshgrid(dspace,dspace);

d1 = d1(:); d2 = d2(:);

d3 = 2*ones(size(d1));

D = [d1, d2, d3];
% D = D./sum(D,2)*3;

f = [0:0.02:0.1 0.15:0.1:0.85 0.9:0.02:1]'; nB = numel(f);
B = [zeros(nB,1) f ones(nB,1)];
B = B./sum(B,2);

RR = safe_rotmat_random(nR);

Rgfo = library.loadSet(44, 'GFOD2');
Resr = library.loadSet(44, 'ESRD2');
Rz   = library.loadSet(44, 'ESRS2');
Rt   = library.loadSet(44, 'Std');

u = util.applyRotMatToVec([0 0 1], Rz);


for i = 1:size(u,1)
    Rx(:,:,i) = util.rotVec2Vec([1 0 0], u(i,:));
end


for i = 1:nB
    bg = rotSet(B(i,:), Rgfo);
    be = rotSet(B(i,:), Resr);
    bt = rotSet(B(i,:), Rt);
    bx = rotSet(B(i,:), Rx);
    bz = rotSet(B(i,:), Rz);

    for j = 1:size(D,1)

        d  = rotSet(D(j,:), RR);

        Sg = exp(-bg*d');
        Se = exp(-be*d');
        St = exp(-bt*d');
        Sx = exp(-bx*d');
        Sz = exp(-bz*d');

        Pg = mean(Sg,1);
        Pe = mean(Se,1);
        Pt = mean(St,1);
        Px = mean(Sx,1);
        Pz = mean(Sz,1);

        SDg(i,j) = std(Pg);
        SDe(i,j) = std(Pe);
        SDt(i,j) = std(Pt);
        SDx(i,j) = std(Px);
        SDz(i,j) = std(Pz);

        Ag(i,j) = mean(Pg);
        Ae(i,j) = mean(Pe);
        At(i,j) = mean(Pt);
        Ax(i,j) = mean(Px);
        Az(i,j) = mean(Pz);

    end
    i
end

CVg = SDg./Ag*100;
CVe = SDe./Ae*100;
CVt = SDt./At*100;
CVx = SDx./Ax*100;
CVz = SDz./Az*100;


%%
figure(1)
clf

plot(f, mean(Ag,2), f, mean(Ae,2), f, mean(Ax,2), f, mean(Az,2))
legend('GFO-D_2', 'ESR-D_2', 'ESR-PTE', 'ESR-LTE')
ylabel('Average signal')
xlabel('b2-factor [1]')

figure(2)
clf

semilogy(f, mean(SDg,2), f, mean(SDe,2), f, mean(SDx,2), f, mean(SDz,2))
legend('GFO', 'ESR', 'ESR-PTE', 'ESR-LTE')
ylabel('SD of signal')
xlabel('b2-factor [1]')


%%
figure(3)
clf

pp = fix.subplot();
pp.nr = 1;
pp.nc = 1;
pp.bo = 0.23;
pp.lo = 0.09;
pp.ro = 0.06;
pp.to = 0.01;

p2 = fix.subplot('tight');
p2.nr = 1;
p2.nc = 6;
p2.to = 0.8;
p2.bo = -0.08;
p2.lo = 0.01;
p2.ro = -0.025;

ut = 5;

fix.subplot(pp,1)

b2f = f./(f+1);
fix.plot.lineArea(b2f, median(CVg,2), prctile(CVg, 100-ut, 2), prctile(CVg, ut, 2));
fix.plot.lineArea(b2f, median(CVe,2), prctile(CVe, 100-ut, 2), prctile(CVe, ut, 2));
fix.plot.lineArea(b2f, median(CVt,2), prctile(CVt, 100-ut, 2), prctile(CVt, ut, 2));
fix.plot.lineArea(b2f, median(CVz,2), prctile(CVz, 100-ut, 2), prctile(CVz, ut, 2), 'color', [1 1 1]*0.5);
fix.plot.lineArea(b2f, median(CVx,2), prctile(CVx, 100-ut, 2), prctile(CVx, ut, 2), 'color', [1 1 1]*0.5, 'linestyle', ':');

hl = legend('GFO', 'ESR $SO(3)$', 't-design', 'ESR $\mathbb{S}^2$-LTE', 'ESR $\mathbb{S}^2$-PTE', 'Location', 'sw', 'interpreter', 'latex', 'box', 'off');

fix.legend.general(hl, 10, 15)

set(gca, 'yscale', 'log')
set(gca, 'XTick', [0:0.1:0.5])
set(gca, 'YTick', [1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])

grid on

ylabel('CV [%]')
xlabel('{\itb}_2/{\itb} [1]')
ylim([0.5e-4 20])
fix.plot.axis

e1 = 1;
e2 = 0.7;
el = [e1 e2 e2 e2 e2 e1];

fi = linspace(0, 1, 6);
for i = 1:6

    ind = round(fi(i)*(size(B,1)-1))+1;

    tmp = diag(B(ind, [3 2 1])+0.1);

    fix.subplot(p2, i)

    fix.plot.superQuadratic(tmp, 'col', [1 1 1]*0.5, 'level', 1, ...
        'eps1', el(i), 'eps2', el(i))

    title ''
    axis off
    axis([-1 1 -1 1 -1 1])
    lighting phong;
    lightangle(-15,30);
    camlight;

end

fix.figure

fix.figure.save('fig_cv_allD_v9', '', [6 4], 600)


%% Is GFO always best

figure(4)
clf

subplot(2,2,1)
semilogy(f, CVe./CVg, 'color', [0 0 0 .1])
title 'ESR/GFO'
yline(1, 'r:')

ylim([0.5 10000])
ylabel('CV(x)/CV(GFO)')

subplot(2,2,2)
semilogy(f, CVt./CVg, 'color', [0 0 0 .1])
title 't-design/GFO'
yline(1, 'r:')

ylim([0.5 10000])

subplot(2,2,3)
semilogy(f, CVx./CVg, 'color', [0 0 0 .1])
title 'ESR-PTE/GFO'
yline(1, 'r:')
ylabel('CV(x)/CV(GFO)')

ylim([0.5 10000])
xlabel('b2factor')

subplot(2,2,4)
semilogy(f, CVz./CVg, 'color', [0 0 0 .1])
title 'ESR-LTE/GFO'
yline(1, 'r:')

ylim([0.5 10000])
xlabel('b2factor')

fix.figure

fix.figure.save('GFO always best (N44)', '', [6 4], 300)


%% Show direction sets

figure(5)

clf
clear vg vx
for i = 1:size(Rgfo,3)
    vg(i,:) = [1 0 0] * Rgfo(:,:,i)';
    vx(i,:) = [1 0 0] * Rx(:,:,i)';
end

pp = fix.subplot('tight');
pp.nr = 2;
pp.nc = 1;

fix.subplot(pp, 1)
plot_gdirs([vx; -vx], abs([vx; -vx]).^2/1.2)

text(0, 1.2, 'ESR 𝕊^2-LTE', 'FontName', 'times', 'HorizontalAlignment','center')
axis off
axis([-1 1 -1 1 -1 1]*1.3)

fix.subplot(pp, 2)
plot_gdirs([vg; -vg], abs([vg; -vg]).^2/1.2, 0, 0, sqrt([61 125 158]/255)*0.95)
text(0, 1.2, 'GFO', 'FontName', 'times', 'HorizontalAlignment','center')
axis off
axis([-1 1 -1 1 -1 1]*1.3)

fix.figure

fix.figure.save('GFO always best (N44)', '', [2 4], 600)


%%
function b = rotSet(beig, R)

b = zeros(size(R,3), 6);

B = diag(beig);

for i = 1:size(R,3)
    b(i,:) = tm_3x3_to_1x6(R(:,:,i) * B * R(:,:,i)');
end

end