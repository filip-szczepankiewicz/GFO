clear

n = [6 12 18 24 32 64];

B = diag([0 1/4 3/4]+0.05);

p = fix.subplot();
p.nr = 2;
p.nc = 3;

clf

for i = 1:numel(n)

    R = library.loadSet(n(i), 'GFOD2');

    fix.subplot(p,i)
    for j = 1:size(R,3)
        fix.plot.superQuadratic( R(:,:,j) * B * R(:,:,j)', 'col', abs([1 0 0]*R(:,:,j)')', 'level', 1, ...
            'eps1', 0.2, 'eps2', 0.5, 'alpha', 0.5);
    end

    title ''

    lighting phong;
    lightangle(-45,30);
    camlight headlight;

    axis([-1 1 -1 1 -1 1])

    campos([1 1 1]*6)

    text(0,0,1.4,['N = ' num2str(n(i))], 'HorizontalAlignment','center', 'FontName','times')

    axis off

end


fix.figure.save('rotationSets', '', [6 4], 600)