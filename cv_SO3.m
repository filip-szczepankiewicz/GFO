function [cv_unw, cv_w, s_unw, s_w] = cv_SO3(E, w, B, D, J, opts)
%CV_SO3  Coefficient of variation of SO(3) powder-averaged signals.
%
%   [cv_unw, cv_w, s_unw, s_w] = CV_SO3(E, w, B, D, J)
%   [cv_unw, cv_w, s_unw, s_w] = CV_SO3(E, w, B, D, J, opts)
%
%   Computes Monte-Carlo estimates of the coefficient of variation (CV) of
%   powder-averaged signals defined as
%       S_i(U) = exp( -trace( B_i * U * D * U.' ) ),
%   where B_i = R_i * B * R_i.' are rotated versions of B and U are random
%   rotations in SO(3). The CV is evaluated over J rotations of D, using both
%   unweighted and weighted averages over the B orientations.
%
%   INPUTS
%     E    : Either an N×3 array of ZYZ Euler angles, or a 3×3×N array of
%            rotation matrices defining orientations R_i.
%     w    : N×1 vector of weights associated with each orientation R_i.
%     B    : 3×3 symmetric matrix.
%     D    : 3×3 diagonal matrix.
%     J    : Number of random (or quasi-uniform) SO(3) rotations of D.
%     opts : (optional) struct with fields:
%            - batch     : Batch size for Monte-Carlo sampling (default 5000)
%            - useParfor : Logical flag to use PARFOR (default false)
%            - seeds     : Vector of RNG seeds, length ≥ number of batches
%            - sample    : 'random' (Haar, default) or 'uniform'
%            - rots      : Optional 3×3×J array of user-supplied rotations
%
%   OUTPUTS
%     cv_unw : Coefficient of variation using uniform (unweighted) average.
%     cv_w   : Coefficient of variation using weighted average w.
%     s_unw  : J×1 vector of unweighted powder-averaged signals.
%     s_w    : J×1 vector of weighted powder-averaged signals.
%
%   EXAMPLE
%     E  = rand(100,3);                 % Euler angles
%     w  = ones(100,1)/100;             % uniform weights
%     B  = eye(3);
%     D  = diag([1 0.5 0.2]);
%     J  = 1e4;
%     [cv_u, cv_w] = cv_SO3(E, w, B, D, J);
%
%   See also PAGEMTIMES, PARFOR.

clear so3_cv_batch;
if nargin < 6, opts = struct; end
if ~isfield(opts,'batch'),     opts.batch = 5000;      end
if ~isfield(opts,'useParfor'), opts.useParfor = false; end
if ~isfield(opts,'seeds'),     opts.seeds = [];        end
if ~isfield(opts,'sample'),     opts.sample = 'random';        end

if ismatrix(E) %Euler angles
    N     = size(E,1);
    Ri = zeros(3,3,N);
    for i = 1:N
        Ri(:,:,i) = util.Rzyz(E(i,1),E(i,2),E(i,3));
    end
elseif ndims(E) == 3  && size(E,1) == 3  && size(E,2) == 3 % rot matrices
    N     = size(E,3);
    Ri = E;
else
    error("wrong input");
end
w     = w(:);
w_unw = ones(N,1)/N;

% --- Powder: Precompute B_i pages and vectorize to 9 x N -------------------------
%B_i = R_i B R_i^T
Bi = pagemtimes(pagemtimes(Ri, repmat(B,1,1,N)), permute(Ri,[2 1 3])); % 3x3xN
Bi9 = reshape(Bi, 9, N);                                                % 9 x N
d   = diag(D).';                                                        % 1 x 3

% --- Split J into batches ------------------------------------------------
BATCH    = opts.batch;
nBatches = ceil(J / BATCH);
ranges   = arrayfun(@(b) ((b-1)*BATCH+1):min(b*BATCH, J), 1:nBatches, 'uni', 0);
if isfield(opts,'rots')
    assert(size(opts.rots,3) == J, 'Incorrect number of rotations supplied')
end
if isempty(opts.seeds)
    rng("shuffle");
    seeds = randi(2^31-1, nBatches, 1);
else
    seeds = opts.seeds(:);
    assert(numel(seeds) >= nBatches, 'opts.seeds must have >= nBatches elements.');
end

s_unw = zeros(J,1);
s_w   = zeros(J,1);

if opts.useParfor
    p = gcp('nocreate'); if isempty(p), parpool; end
    tmp_unw = cell(nBatches,1);
    tmp_w   = cell(nBatches,1);
    parfor b = 1:nBatches
        [tmp_unw{b}, tmp_w{b}] = so3_cv_batch(Bi9, d, w, w_unw, ranges{b}, seeds(b),J,opts);
    end
    for b = 1:nBatches
        s_unw(ranges{b}) = tmp_unw{b};
        s_w(ranges{b})   = tmp_w{b};
    end
else
    for b = 1:nBatches
        [s_unw(ranges{b}), s_w(ranges{b})] = so3_cv_batch(Bi9, d, w, w_unw, ranges{b}, seeds(b),J,opts);
    end
end
cv_unw = std(s_unw) / mean(s_unw);
cv_w   = std(s_w)   / mean(s_w);
end

% ===================== subfunctions (OK for parfor) ======================

function [y_unw, y_w] = so3_cv_batch(Bi9, d, w, w_unw, range, seed,J,opts)
% CV of powder average over rotations of D
% Single batch (fully vectorized). No nested functions => parfor-safe.
persistent U_cache;
rng(seed);

M = numel(range);

% Sample  rotations
if isfield(opts,'rots') %custom
    U = opts.rots(:,:,range);
elseif strcmp(opts.sample, 'random') %Haar
    U = util.randrotations(M);
elseif strcmp(opts.sample, 'uniform') %Hopf-kronecker quasi uniform
    if isempty(U_cache) 
        Q        = util.hopf_kronecker(J);  
        U_cache  = util.quat_to_R(Q);           
    end
    % Reuse cached rotations for this batch
    U = U_cache(:,:,range);
end

% Dj(:,:,m) = sum_k d(k) * u_k u_k^T  (vectorized)
u1 = squeeze(U(:,1,:));  % 3xM
u2 = squeeze(U(:,2,:));
u3 = squeeze(U(:,3,:));
Dj = d(1) * pagemtimes(reshape(u1,3,1,M), reshape(u1,1,3,M)) ...
    + d(2) * pagemtimes(reshape(u2,3,1,M), reshape(u2,1,3,M)) ...
    + d(3) * pagemtimes(reshape(u3,3,1,M), reshape(u3,1,3,M));   % 3x3xM
Dj9 = reshape(Dj, 9, M);                                        % 9 x M

% T = (Bi9)' * Dj9   (N x M);  S = exp(-T)
T = Bi9.' * Dj9;
S = exp(-T);

y_unw = (w_unw.' * S).';   % M x 1
y_w   = (w.'     * S).';
end
