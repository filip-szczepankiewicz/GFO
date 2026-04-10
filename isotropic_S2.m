function [x, f,flag,out] = isotropic_S2(N, mode, Lmax, Vdeg, x0)
% function [x, f,flag,out] = isotropic_S2(N, mode, Lmax, Vdeg, x0)

%anti-podal flag to identify anti-podal points
use_proj = true;

% fmincon options (interior-point; increase if you want tighter solves)
opt = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','none', ...   % 'iter' if you want logs
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 1e5, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

% Initial step: use spherical Fibonacci as a good seed
if nargin < 6 || isempty(x0)
    x0 = util.sph_fibonacci(N);
end

switch mode
    case {'electro', 'electroS2', 'ESR', 'ESRS2'}
        obj = @(x)electroS2_cost_fct(x, use_proj);
    case {'GFO', 'GFOS2'}
        obj = @(x)GFOS2_cost_fct(x, Lmax, Vdeg);
    otherwise
        error('Mode is not recognized!')
end

lb = -ones(size(x0));
ub =  ones(size(x0));
problem = struct('objective', obj, ...
    'x0', x0, ...
    'Aineq', [], 'bineq', [], ...
    'Aeq', [], 'beq', [], ...
    'lb', lb, 'ub', ub, ...
    'nonlcon', [], ...
    'solver', 'fmincon', ...
    'options', opt);

[x, f,flag,out] = fmincon(problem);
x = util.normalize_rows(x);
x = util.canon_antipodes(x);
end


