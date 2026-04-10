function u = applyRotMatToVec(v, R)
% function u = util.applyRotMatToVec(v, R)
%
% Proper use of rotation matrices depends on the vector to be rotated.
% For column vectors [3x1] we use u = R*v
% For row vectors [1x3] we use u = v*R' = (R*v')'
%
% Here we ensure the input form by first making it a column vector, but
% then we assume the output is the more conventional nx3 format.

v = v(:);

for i = 1:size(R,3)

    u(i,:) = (R(:,:,i) * v)';

end


