function R = loadSet(n, mode)
% function R = library.loadSet(n, mode)

basePath = fileparts(mfilename('fullpath'));

fn = [basePath filesep mode '_' sprintf('%04d', n)];

R = util.readRotMatsFromFile(fn);
