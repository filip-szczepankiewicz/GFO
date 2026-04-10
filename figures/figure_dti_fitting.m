clear

nD = 15;
nR = 300;

dspace = linspace(0,2,nD);

[d1,d2] = meshgrid(dspace,dspace);

d1 = d1(:); d2 = d2(:);

d3 = 2*ones(size(d1));

D = [d1, d2, d3];
D = D./sum(D,2)*3;

% f = [0:0.01:0.2 0.25:0.1:0.75 0.8:0.01:1]'; nB = numel(f);
f = [0 0.1 0.2:0.2:0.8 0.9 1]'; nB = numel(f);
B = [zeros(nB,1) f ones(nB,1)];
B = B./sum(B,2);

RR = safe_rotmat_random(nR);

% Rgfo = GFO_generateSet(44, 'GFOD2');
% Resr = GFO_generateSet(44, 'electroD2');
load('gfod2_044.mat');
load('esrd2_044.mat');
load('forFilip.mat'); euA = rotQ{4};

for i = 1:size(euA,1)
    Rt(:,:,i) = util.Rzyz(euA(i,1),euA(i,2),euA(i,3));
end

[Rx, ~, Rz] = fsz_createSet(44);

opt = dti_lls_opt;

for i = 1:nB
    bg = [zeros(1,6); rotSet(B(i,:), Rgfo)];
    be = [zeros(1,6); rotSet(B(i,:), Resr)];
    bt = [zeros(1,6); rotSet(B(i,:), Rt)];
    bx = [zeros(1,6); rotSet(B(i,:), Rx)];
    bz = [zeros(1,6); rotSet(B(i,:), Rz)];

    xpsg = mdm_xps_from_bt(bg*1e9);
    xpse = mdm_xps_from_bt(be*1e9);
    xpst = mdm_xps_from_bt(bt*1e9);
    xpsx = mdm_xps_from_bt(bx*1e9);
    xpsz = mdm_xps_from_bt(bz*1e9);

    for j = 1:size(D,1)

        d  = rotSet(D(j,:), RR);

        rn = randn(45, nR)/20;

        Sg = exp(-bg*d')+rn;
        Se = exp(-be*d')+rn;
        St = exp(-bt*d')+rn;
        Sx = exp(-bx*d')+rn;
        Sz = exp(-bz*d')+rn;

        for k = 1:size(Sg,2)
            tmpg(k,:) = dti_lls_1d_data2fit(Sg(:,k), xpsg, opt);
            tmpe(k,:) = dti_lls_1d_data2fit(Se(:,k), xpse, opt);
            tmpt(k,:) = dti_lls_1d_data2fit(St(:,k), xpst, opt);
            tmpx(k,:) = dti_lls_1d_data2fit(Sx(:,k), xpsx, opt);
            tmpz(k,:) = dti_lls_1d_data2fit(Sz(:,k), xpsz, opt);
        end

        parg = tm_1x6_to_tpars(tmpg(:,2:end));
        mdg(i,j) = mean(parg.trace/3);
        fag(i,j) = mean(parg.fa);
        mdgv(i,j) = std(parg.trace/3);
        fagv(i,j) = std(parg.fa);


        pare = tm_1x6_to_tpars(tmpe(:,2:end));
        mde(i,j) = mean(pare.trace/3);
        fae(i,j) = mean(pare.fa);
        mdev(i,j) = std(pare.trace/3);
        faev(i,j) = std(pare.fa);


        part = tm_1x6_to_tpars(tmpt(:,2:end));
        mdt(i,j) = mean(part.trace/3);
        fat(i,j) = mean(part.fa);
        mdtv(i,j) = std(part.trace/3);
        fatv(i,j) = std(part.fa);


        parx = tm_1x6_to_tpars(tmpx(:,2:end));
        mdx(i,j) = mean(parx.trace/3);
        fax(i,j) = mean(parx.fa);
        mdxv(i,j) = std(parx.trace/3);
        faxv(i,j) = std(parx.fa);


        parz = tm_1x6_to_tpars(tmpz(:,2:end));
        mdz(i,j) = mean(parz.trace/3);
        faz(i,j) = mean(parz.fa);
        mdzv(i,j) = std(parz.trace/3);
        fazv(i,j) = std(parz.fa);

    end
    i
end



%% MEAN
figure(1)
clf
ut = 25;

subplot(2,2,1)
fix.plot.lineArea(f, median(mdg,2), prctile(mdg, 100-ut, 2), prctile(mdg, ut, 2));
fix.plot.lineArea(f, median(mde,2), prctile(mde, 100-ut, 2), prctile(mde, ut, 2));
fix.plot.lineArea(f, median(mdt,2), prctile(mdt, 100-ut, 2), prctile(mdt, ut, 2));
fix.plot.lineArea(f, median(mdx,2), prctile(mdx, 100-ut, 2), prctile(mdx, ut, 2), 'color', [1 1 1]*0.5);
fix.plot.lineArea(f, median(mdz,2), prctile(mdz, 100-ut, 2), prctile(mdz, ut, 2), 'color', [1 1 1]*0.5, 'linestyle', ':');
title 'mean MD'


subplot(2,2,3)
fix.plot.lineArea(f, median(mdgv,2), prctile(mdgv, 100-ut, 2), prctile(mdgv, ut, 2));
fix.plot.lineArea(f, median(mdev,2), prctile(mdev, 100-ut, 2), prctile(mdev, ut, 2));
fix.plot.lineArea(f, median(mdtv,2), prctile(mdtv, 100-ut, 2), prctile(mdtv, ut, 2));
fix.plot.lineArea(f, median(mdxv,2), prctile(mdxv, 100-ut, 2), prctile(mdxv, ut, 2), 'color', [1 1 1]*0.5);
fix.plot.lineArea(f, median(mdzv,2), prctile(mdzv, 100-ut, 2), prctile(mdzv, ut, 2), 'color', [1 1 1]*0.5, 'linestyle', ':');
title 'std MD'


subplot(2,2,2)
fix.plot.lineArea(f, median(fag,2), prctile(fag, 100-ut, 2), prctile(fag, ut, 2));
fix.plot.lineArea(f, median(fae,2), prctile(fae, 100-ut, 2), prctile(fae, ut, 2));
fix.plot.lineArea(f, median(fat,2), prctile(fat, 100-ut, 2), prctile(fat, ut, 2));
fix.plot.lineArea(f, median(fax,2), prctile(fax, 100-ut, 2), prctile(fax, ut, 2), 'color', [1 1 1]*0.5);
fix.plot.lineArea(f, median(faz,2), prctile(faz, 100-ut, 2), prctile(faz, ut, 2), 'color', [1 1 1]*0.5, 'linestyle', ':');
title 'mean FA'


subplot(2,2,4)
fix.plot.lineArea(f, median(fagv,2), prctile(fagv, 100-ut, 2), prctile(fagv, ut, 2));
fix.plot.lineArea(f, median(faev,2), prctile(faev, 100-ut, 2), prctile(faev, ut, 2));
fix.plot.lineArea(f, median(fatv,2), prctile(fatv, 100-ut, 2), prctile(fatv, ut, 2));
fix.plot.lineArea(f, median(faxv,2), prctile(faxv, 100-ut, 2), prctile(faxv, ut, 2), 'color', [1 1 1]*0.5);
fix.plot.lineArea(f, median(fazv,2), prctile(fazv, 100-ut, 2), prctile(fazv, ut, 2), 'color', [1 1 1]*0.5, 'linestyle', ':');
title 'std FA'


hl = legend('GFO', 'ESR $SO(3)$', 't-design', 'ESR  $\mathbb{S}^2$-PTE', 'ESR  $\mathbb{S}^2$-LTE', 'Location', 'nw', 'interpreter', 'latex', 'box', 'off');

fix.legend.general(hl, 10, 15)


%%
function b = rotSet(beig, R)

b = zeros(size(R,3), 6);

B = diag(beig);

for i = 1:size(R,3)
    b(i,:) = tm_3x3_to_1x6(R(:,:,i) * B * R(:,:,i)');
end

end