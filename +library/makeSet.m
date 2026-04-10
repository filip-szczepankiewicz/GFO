function fnl = makeSet(n, mode)
% function fnl = library.makeSet(n, mode)

basePath = fileparts(mfilename('fullpath'));

for i = 1:numel(n)

    ofn    = [basePath filesep mode '_' sprintf('%04d', n(i))];
    fnl{i} = ofn;

    if exist(ofn, 'file')
        disp(['Skipping: ' ofn])
        continue
    end

    rotMats = GFO_generateSet(n(i), mode);
    util.writeRotMatsToFile(rotMats, ofn);

end

