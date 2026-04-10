function R = rotVec2Vec(a, b)
% function R = util.rotVec2Vec(a, b)
% Rotation matrix that rotates vector a into vector b.
%
%   R = vec2vec_rotmat(a, b)
%
%   Returns a 3x3 rotation matrix R such that R*a is parallel to b
%   (and in the same direction). Vectors need not be unit length.
%
%   Uses the Rodrigues formula:
%       R = I + [v]x + [v]x^2 * (1 / (1 + c))
%   where v = a x b, c = a . b, and [v]x is the skew-symmetric
%   cross-product matrix of v.
%
%   Edge case: if a and b are antiparallel (c ≈ -1), the rotation axis
%   is undefined. This is handled by finding an arbitrary perpendicular
%   axis and rotating by pi.

a = a(:) / norm(a);
b = b(:) / norm(b);

v = cross(a, b);        % axis direction (unnormalized)
c = dot(a, b);          % cos(angle)

if abs(c + 1) < 1e-10   % antiparallel: rotate pi around any perp axis
    perp = null(a');    % any vector perpendicular to a
    k = perp(:, 1);
    K = skew(k);
    R = eye(3) + 2 * K^2;
    return
end

K = skew(v);
R = eye(3) + K + K^2 / (1 + c);

end


function K = skew(v)
% Skew-symmetric cross-product matrix of v
K = [  0   -v(3)  v(2);
      v(3)   0   -v(1);
     -v(2)  v(1)   0  ];
end