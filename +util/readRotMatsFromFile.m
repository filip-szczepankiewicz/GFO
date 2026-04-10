function R = readRotMatsFromFile(filename)
% function R = readRotMatsFromFile(filename)
%
% Input:
%   filename - text file written by writeRotationMatricesToFile
%
% Output:
%   R        - 3x3xN array of rotation matrices
%
% Expected file format: N rows, 9 columns, row-major flattening of 3x3 matrix
%   Col:  1    2    3    4    5    6    7    8    9
%         R11  R12  R13  R21  R22  R23  R31  R32  R33
%
%   where Rij = R(i,j,:), i = row index, j = column index

fid = fopen(filename, 'r');
if fid == -1
    error('Could not open file: %s', filename);
end

% Skip comment lines starting with '#'
R_flat = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || line(1) == '#'
        continue
    end
    R_flat(end+1, :) = sscanf(line, '%g')'; %#ok<AGROW>
end

fclose(fid);

n = size(R_flat, 1);
assert(size(R_flat, 2) == 9, 'Expected 9 columns per row, got %d.', size(R_flat, 2));

% Inverse of permute(R,[3 1 2]) + reshape: go from Nx9 -> Nx3x3 -> 3x3xN
R = permute(reshape(R_flat, n, 3, 3), [2 3 1]);

end