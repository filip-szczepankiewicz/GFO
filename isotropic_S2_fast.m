function [q, output] = isotropic_S2_fast(N, metric)
% function [q, output] = isotropic_S2_fast(N, metric)
%
% Electrostatic repulsion on RP2: N directions such that [q; -q] is
% uniformly distributed on S2. Uses fminunc with analytic gradient;
% sphere constraint is eliminated by normalising inside the objective.
%
% Note that this funciton is vastly faster than isotorpic_S2, but has a
% worse performance wrt the powder average CV.

if nargin < 2, metric = 'geodesic'; end

switch metric
    case 'euclid'
        energyFcn = @(x) energy_euclid(x, N);
    case 'geodesic'
        energyFcn = @(x) energy_geodesic(x, N);
    otherwise
        error('Invalid metric: %s', metric);
end

opts = optimoptions('fminunc', ...
    'Algorithm',                'quasi-newton', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations',   4e3 * N^2, ...
    'MaxIterations',            2e3 * N, ...
    'Display',                  'none');

nRestarts = 5;
bestFval  = inf;
bestX     = [];
output    = [];

for r = 1:nRestarts
    x0 = randn(N, 3);                          % random Cartesian init
    x0 = x0 ./ vecnorm(x0, 2, 2);             % project to sphere

    [x, fval, exitflag, out] = fminunc(energyFcn, x0(:), opts);

    if fval < bestFval
        bestFval        = fval;
        bestX           = x;
        output          = out;
        output.exitflag = exitflag;
        output.fval     = fval;
    end
    if exitflag <= 0
        fprintf('Restart %d: exitflag %d\n', r, exitflag);
    end
end

% Recover unit vectors
q = reshape(bestX, N, 3);
q = q ./ vecnorm(q, 2, 2);

end

% -------------------------------------------------------------------------
function [ii, jj] = rp2_pairs(N)
% All upper-triangle pairs for 2N points, excluding self-antipodal (i, i+N).
M = 2*N;
[ii_all, jj_all] = find(triu(true(M), 1));
keep = (jj_all ~= ii_all + N);
ii = ii_all(keep);
jj = jj_all(keep);
end

% -------------------------------------------------------------------------
function [E, dEdx] = energy_geodesic(x_flat, N)
% Geodesic energy on RP2. Sphere constraint handled by normalising x first.

X   = reshape(x_flat, N, 3);
n   = vecnorm(X, 2, 2);
XYZ = X ./ n;                                % project onto S2

R   = [XYZ; -XYZ];                           % 2N x 3
M   = 2*N;

C    = R * R';                               % 2N x 2N dot products
C    = max(-1+1e-9, min(1-1e-9, C));
absC = abs(C);

[ii, jj] = rp2_pairs(N);
nPairs   = numel(ii);
idx      = sub2ind([M M], ii, jj);

c_ij = absC(idx);
D_ij = acos(c_ij) / pi;
E    = mean(1 ./ max(D_ij, eps));

% Gradient w.r.t. XYZ
w_D  = -1 ./ (D_ij.^2) / nPairs;
sgn  = sign(C(idx));
dDdc = -sgn ./ (pi * sqrt(max(1 - c_ij.^2, eps)));
w    = w_D .* dDdc;                          % dE/d(c_ij)

dEdXYZ = zeros(M, 3);
for k = 1:3
    dEdXYZ(:,k) = dEdXYZ(:,k) ...
        + accumarray(ii, w .* R(jj,k), [M 1]) ...
        + accumarray(jj, w .* R(ii,k), [M 1]);
end
dEdXYZ = dEdXYZ(1:N,:) - dEdXYZ(N+1:end,:);  % chain rule through -XYZ

% Chain rule through normalisation: project onto tangent plane
% dE/dX_i = (I - u_i u_i') / ||x_i|| * dE/dXYZ_i
dEdX = zeros(N, 3);
for i = 1:N
    ui       = XYZ(i,:)';
    P        = (eye(3) - ui*ui') / n(i);      % tangent plane projector
    dEdX(i,:) = (P * dEdXYZ(i,:)')';
end

dEdx = dEdX(:);
end

% -------------------------------------------------------------------------
function [E, dEdx] = energy_euclid(x_flat, N)
% Euclidean energy on RP2.

X   = reshape(x_flat, N, 3);
n   = vecnorm(X, 2, 2);
XYZ = X ./ n;

R   = [XYZ; -XYZ];
M   = 2*N;

A  = permute(R, [1 3 2]);
B  = permute(R, [3 1 2]);
dR = A - B;                                  % M x M x 3
D2 = sum(dR.^2, 3);

[ii, jj] = rp2_pairs(N);
nPairs   = numel(ii);
idx      = sub2ind([M M], ii, jj);

d_ij = sqrt(D2(idx)) / 2;
E    = mean(1 ./ max(d_ij, eps));

w = -1 ./ (d_ij.^3) / nPairs / 4;

dEdR = zeros(M, 3);
for k = 1:3
    dr_k    = squeeze(dR(:,:,k));
    contrib = w .* dr_k(idx);
    dEdR(:,k) = dEdR(:,k) ...
        + accumarray(ii, contrib, [M 1]) ...
        - accumarray(jj, contrib, [M 1]);
end
dEdXYZ = dEdR(1:N,:) - dEdR(N+1:end,:);

dEdX = zeros(N, 3);
for i = 1:N
    ui        = XYZ(i,:)';
    P         = (eye(3) - ui*ui') / n(i);
    dEdX(i,:) = (P * dEdXYZ(i,:)')';
end

dEdx = dEdX(:);
end