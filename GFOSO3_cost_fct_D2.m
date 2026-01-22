function val = GFOSO3_cost_fct_D2(x, N, Lmax, bandWeight)

%   GFOSO3_COST_FCT_D2  D2-symmetric variant of GFOSO3_COST_FCT via group averaging.
%
%   val = GFOSO3_COST_FCT_D2(x, N, Lmax, bandWeight)
%
%   This function is the dihedral-symmetry (D2) extension of
%   GFOSO3_COST_FCT. It computes the same band-limited SO(3) energy, but
%   enforces *right D2 symmetry* by averaging the cost over the action of
%   the dihedral group D2 on the second argument.
%
%   In representation-theoretic terms, the character
%       chi^ell(h_j h_k^{-1})
%   used in GFOSO3_COST_FCT is replaced by its group-averaged form
%       (1/|D2|) * sum_{K in D2} chi^ell(h_j K h_k^{-1}),
%   which projects the cost onto the D2-invariant subspace.
%
%   INPUTS / OUTPUT
%     Identical to GFOSO3_COST_FCT:
%       x, N, Lmax, bandWeight  -> see GFOSO3_COST_FCT for definitions.
%
%     - D2 invariance:
%         The resulting energy is invariant under right multiplication of
%         all samples by any element of D2, appropriate for objects or
%         orientations with dihedral symmetry.
%   SEE ALSO
%     GFOSO3_COST_FCT

Q = reshape(x, [N,4]);
Q = util.normalize_rows(Q);
nR = size(Q,1);

% D2 lift inside Q8: {1, i, j, k}  
Kset = [ 1 0 0 0;   % 1
         0 1 0 0;   % i
         0 0 1 0;   % j
         0 0 0 1];  % k
nK = size(Kset,1);

val = 0;

for kk = 1:nK
    % Right-multiply second factor by K: q_k -> q_k ⊗ K  (inverse = itself for these K)
    QK = qmul_right(Q, Kset(kk,:));        % nR×4

    % Pairwise cos(beta/2) = |<q_j, q_k ⊗ K>| in [0,1]
    S = Q * QK.';                          % nR×nR
    S = sqrt(S.^2);                        % abs, stable
    S = max(min(S, 1-1e-12), 0);           % clamp to [0,1)

    % Recurrence wants diagonal = 1; we later subtract diag(U) so jk self-terms are ignored
    S(1:nR+1:end) = 1;

    % Chebyshev U_k(S) recurrence
    U_prev = ones(nR);     % U0
    U_curr = 2*S;          % U1

    % advance to U2 (k=2), no contribution (ell=1)
    U_next = 2*S .* U_curr - U_prev;
    U_prev = U_curr; U_curr = U_next;

    valK = 0;
    for k = 3:2*Lmax
        U_next = 2*S .* U_curr - U_prev;
        if mod(k,2)==0
            ell = k/2;   % ell = 2..Lmax
            sum_all = sum(U_next(:)) - sum(diag(U_next));  % drop self-terms
            valK = valK + bandWeight(ell+1) * sum_all;
        end
        U_prev = U_curr; U_curr = U_next;
    end

    val = val + valK;
end
% Average over K in D2
val = val / nK;

end

function Qout = qmul_right(Q, r)
% Qout = Q ⊗ r, with quaternions in [t x y z]
t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);
rt = r(1);  rx = r(2);  ry = r(3);  rz = r(4);

Qout = zeros(size(Q),'like',Q);
Qout(:,1) = t*rt - x*rx - y*ry - z*rz;
Qout(:,2) = t*rx + x*rt + y*rz - z*ry;
Qout(:,3) = t*ry - x*rz + y*rt + z*rx;
Qout(:,4) = t*rz + x*ry - y*rx + z*rt;
end
