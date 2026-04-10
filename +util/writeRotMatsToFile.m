function writeRotMatsToFile(R, filename)
% function writeRotMatsToFile(R, filename)
%
% Inputs:
%   R        - 3x3xN array of rotation matrices
%   filename - output text file path
%
% Output format: N rows, 9 columns, row-major flattening of each 3x3 matrix
%   Col:  1    2    3    4    5    6    7    8    9
%         R11  R12  R13  R21  R22  R23  R31  R32  R33
%
%   where Rij = R(i,j,:), i = row index, j = column index

[d1, d2, n] = size(R);
assert(d1 == 3 && d2 == 3, 'Input must be 3x3xN.');

% Reshape: permute to Nx3x3, then reshape to Nx9 (row-major per matrix)
R_flat = reshape(permute(R, [3 1 2]), n, 9);

fid = fopen(filename, 'w');
if fid == -1
    error('Could not open file: %s', filename);
end

% Write header
fprintf(fid, '# Rotation matrices: %d entries, row-major flattening of 3x3 -> 1x9\n', n);
fprintf(fid, '# https://github.com/Neurophysics-CFIN/GFO\n');
fprintf(fid, '# Col:  1    2    3    4    5    6    7    8    9\n');
fprintf(fid, '#       R11  R12  R13  R21  R22  R23  R31  R32  R33\n');
fprintf(fid, '# where Rij = R(i,j), i = row index, j = column index\n\n');

for i = 1:n
    fprintf(fid, '%14.8f', R_flat(i, :));
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Wrote %d rotation matrices to %s\n', n, filename);
end