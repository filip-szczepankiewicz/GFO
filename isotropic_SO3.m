function [t,x,y,z,fval,exitflag,output] = isotropic_SO3(N, cost, Lmax, bandWeight,x0_in)
%   ISOTROPIC_SO3  Optimize an approximately isotropic set of rotations in SO(3).
%
%   [t,x,y,z,fval,exitflag,output] = ISOTROPIC_SO3(N, cost, Lmax, bandWeight)
%   [t,x,y,z,fval,exitflag,output] = ISOTROPIC_SO3(N, cost, Lmax, bandWeight, x0_in)
%
%   Produces N rotations (as unit quaternions) by minimizing a user-selected
%   isotropy cost on SO(3) using unconstrained optimization (fminunc). A
%   deterministic Hopf–Kronecker S^3 lattice initializer is used by default,
%   and the final quaternions are normalized and antipodally canonicalized
%   (q ≡ -q) for SO(3).
%
%   INPUTS
%     N          : Number of rotations/samples.
%     cost       : Cost-function selector (string):
%                    "GFO"         -> GFOSO3_cost_fct
%                    "GFOD2"       -> GFOSO3_cost_fct_D2
%                    "electro"     -> electroSO3_cost_fct   (base only)
%                    "electroD2"   -> electroSO3_cost_fct_D2
%     Lmax       : Maximum bandlimit for GFO-based costs (ignored by electro*).
%     bandWeight : Band weights for GFO-based costs. See GFOSO3_cost_fct for
%                  indexing/interpretation (weights for ℓ=2..Lmax are used).
%                  (ignored by electro*).    
%     x0_in      : (optional) Initial iterate as a (4N-by-1) vector or any
%                  shape convertible to a column vector. If omitted/empty,
%                  a deterministic Hopf–Kronecker initializer is used.
%
%   OUTPUTS
%     t,x,y,z    : N-by-1 quaternion components of the optimized rotations,
%                  in convention [t x y z] with t the scalar part.
%     fval       : Final objective value returned by fminunc.
%     exitflag   : fminunc termination flag.
%     output     : fminunc output structure (iterations, evaluations, etc.).
%
%   DETAILS
%     - Initialization: util.hopf_kronecker(N) provides a deterministic
%       S^3 lattice, then the components are stacked into x0.
%     - Optimization: uses fminunc with quasi-newton updates and central
%       finite-difference gradients (no analytic gradient provided).
%     - Postprocessing: reshapes xbest into N-by-4, normalizes rows, and
%       applies util.canon_antipodes to enforce q ≡ -q.
%
%   EXAMPLE
%     N = 24;
%     Lmax = 8;
%     lvals = 0:Lmax;
%     Vdeg = (2*lvals+1).^2.*util.sobolev(8,8,Lmax).^2;
%     [t,x,y,z,fval] = isotropic_SO3(N,"GFO",Lmax,Vdeg.^2);
%
%   SEE ALSO
%     GFOSO3_cost_fct, GFOSO3_cost_fct_D2, electroSO3_cost_fct,
%     util.hopf_kronecker, fminunc


% --- Deterministic initializer (Hopf–Kronecker on S³, antipodal-canonical) ---
if nargin >=5 && ~isempty(x0_in)
    x0 = x0_in(:);
else
    [t,x,y,z] = util.hopf_kronecker(N);
    x0 = cat(1,t(:),x(:),y(:),z(:));
end
if strcmp(cost,"GFO")
    % --- Objective handle (numeric, bounds-only; normalizes inside) ---
    obj = @(x) GFOSO3_cost_fct(x,N,Lmax,bandWeight);
elseif strcmp(cost,"GFOD2")
    % --- Objective handle (numeric, bounds-only; normalizes inside) ---
    obj = @(x) GFOSO3_cost_fct_D2(x,N,Lmax,bandWeight);
elseif strcmp(cost,"electro")
    obj = @(x) electroSO3_cost_fct(x,N,false); %no addpen
elseif strcmp(cost,"electroD2")
    obj = @(x) electroSO3_cost_fct_D2(x,N);
else
    error('Invalid cost function specified. Choose either "GFO", "GFOD2" or "electro" or "electroD2"');
end

opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'MaxIterations', 2000, ...             
    'MaxFunctionEvaluations', 2e6, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'Display','off', ...
    'FiniteDifferenceType','central', ...
    'UseParallel', false);             

                  
    [xbest,fval,exitflag,output] = fminunc(obj, x0', opts);

% --- Extract unit quaternions from best vector (normalize + canonicalize) ---
Q = reshape(xbest, [N,4]);
Q = util.normalize_rows(Q);
Q = util.canon_antipodes(Q);

t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);

end
