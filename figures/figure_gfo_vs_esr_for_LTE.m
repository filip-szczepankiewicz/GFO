% This code requires the fix matlab repository
% https://github.com/filip-szczepankiewicz/fix_matlab

% ESRS2 LTE is best for N = 6, 7, 11, 22, and for PTE it is N = 6 and 7

clear

n = [6:1:64];

nD = 30;
nR = 1600;

dspace = linspace(0,2,nD);

[d1,d2] = meshgrid(dspace,dspace);
d1 = d1(:); d2 = d2(:);
d3 = 2*ones(size(d1));
D  = [d1, d2, d3];

B = [0 0 1];

RR = util.randrotations(nR);


for i = 1:numel(n)

    Rg = library.loadSet(n(i), 'GFOD2');
    Re = library.loadSet(n(i), 'ESRD2');
    Rz = library.loadSet(n(i), 'ESRS2');

    u = util.applyRotMatToVec([0 0 1], Rz);

    for j = 1:size(u,1)
        Rx(:,:,j) = util.rotVec2Vec([1 0 0], u(j,:));
    end

    bg = rotSet(B, Rg);
    be = rotSet(B, Re);
    bx = rotSet(B, Rx);
    bz = rotSet(B, Rz);

    for j = 1:size(D,1)

        d  = rotSet(D(j,:), RR);

        Sg = exp(-bg*d');
        Se = exp(-be*d');
        Sx = exp(-bx*d');
        Sz = exp(-bz*d');

        Pg = mean(Sg,1);
        Pe = mean(Se,1);
        Px = mean(Sx,1);
        Pz = mean(Sz,1);

        SDg(i,j) = std(Pg);
        SDe(i,j) = std(Pe);
        SDx(i,j) = std(Px);
        SDz(i,j) = std(Pz);

        Ag(i,j) = mean(Pg);
        Ae(i,j) = mean(Pe);
        Ax(i,j) = mean(Px);
        Az(i,j) = mean(Pz);

    end
    i
end


CVg = SDg./Ag*100;
CVe = SDe./Ae*100;
CVx = SDx./Ax*100;
CVz = SDz./Az*100;


%%
figure(1)
clf

f = n;

plot(f, mean(Ag,2), f, mean(Ae,2), f, mean(Ax,2), f, mean(Az,2))
legend('GFO-D_2', 'ESR-D_2', 'ESR-PTE', 'ESR-LTE')
ylabel('Average signal')
xlabel('N [1]')

figure(2)
clf

semilogy(f, mean(SDg,2), f, mean(SDe,2), f, mean(SDx,2), f, mean(SDz,2))
legend('GFO', 'ESR', 'ESR-PTE', 'ESR-LTE')
ylabel('SD of signal')
xlabel('N [1]')


%%
figure(3)
clf

pp = fix.subplot();
pp.nr = 1;
pp.nc = 1;
pp.bo = 0.15;
pp.lo = 0.12;

ut = 5;

fix.subplot(pp,1)

x = n;
fix.plot.lineArea(x, median(CVg,2), prctile(CVg, 100-ut, 2), prctile(CVg, ut, 2));
fix.plot.lineArea(x, median(CVe,2), prctile(CVe, 100-ut, 2), prctile(CVe, ut, 2));
fix.plot.lineArea(x, median(CVz,2), prctile(CVz, 100-ut, 2), prctile(CVz, ut, 2), 'color', [1 1 1]*0.5);

hl = legend('GFO', 'ESR $SO(3)$', 'ESR $\mathbb{S}^2$-LTE', 'ESR $\mathbb{S}^2$-PTE', 'Location', 'sw', 'interpreter', 'latex', 'box', 'off');

fix.legend.general(hl, 10, 15)

set(gca, 'yscale', 'log')
set(gca, 'YTick', [1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])

grid on

ylabel('CV [%]')
xlabel('N [1]')
ylim([0.5e-4 40])
fix.plot.axis

fix.figure

fix.figure.save('fig_cv_vs_Nrot', '', [6 4], 600)


%% Is GFO always best

figure(4)
clf

subplot(2,2,1)
semilogy(f, CVe./CVg, 'color', [0 0 0 .1])
title 'ESR/GFO'
yline(1, 'r:')

ylim([0.01 10000])
ylabel('CV(x)/CV(GFO)')

% subplot(2,2,2)
% semilogy(f, CVt./CVg, 'color', [0 0 0 .1])
% title 't-design/GFO'
% yline(1, 'r:')
% 
% ylim([0.5 10000])

subplot(2,2,3)
semilogy(f, CVx./CVg, 'color', [0 0 0 .1])
title 'ESR-PTE/GFO'
yline(1, 'r:')
ylabel('CV(x)/CV(GFO)')

ylim([0.01 10000])
xlabel('N [1]')

subplot(2,2,4)
semilogy(f, CVz./CVg, 'color', [0 0 0 .1])
title 'ESR-LTE/GFO'
yline(1, 'r:')

ylim([0.01 10000])
xlabel('N [1]')

fix.figure

fix.figure.save('GFO is best after N23', '', [6 4], 300)


%%
function b = rotSet(beig, R)

b = zeros(size(R,3), 6);

B = diag(beig);

for i = 1:size(R,3)
    b(i,:) = tm_3x3_to_1x6(R(:,:,i) * B * R(:,:,i)');
end

end