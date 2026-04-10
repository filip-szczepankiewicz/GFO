function btensor_triangle(colorFcn, opts)
% btensor_triangle  Draw a b-tensor shape triangle colored by a triaxiality metric.
%
%  btensor_triangle()
%  btensor_triangle(@myFcn)
%  btensor_triangle(@myFcn, 'nPts', 100, 'cmap', 'viridis', 'title', 'My metric')
%
%  Corner eigenvalues (trace=1, sorted descending):
%    Linear   (top):          [1,   0,   0  ]
%    Planar   (bottom-left):  [0.5, 0.5, 0  ]
%    Spherical(bottom-right): [1/3, 1/3, 1/3]
%
%  colorFcn receives [λ1 λ2 λ3] sorted descending and returns a scalar.
%  Default: (λ1-λ2)*(λ2-λ3), normalised to [0,1] → peaks at maximally triaxial.

arguments
    colorFcn  (1,1) function_handle = @(e) (e(1)-e(2)) * (e(2)-e(3))
    opts.nPts       (1,1) double  = 100    % vertices per edge
    opts.cmap             double  = jet
    opts.clim             double = []
    opts.title            char   = ''
    opts.showCornerLabels (1,1) logical = true
    opts.showColorbar     (1,1) logical = true
end

% --- Corner positions (2-D) ---
pL = [0.5,  sqrt(3)/2];   % Linear    (top)
pP = [0,    0         ];   % Planar    (bottom-left)
pS = [1,    0         ];   % Spherical (bottom-right)

% --- Corner eigenvalues ---
eL = [1,   0,   0  ];
eP = [0.5, 0.5, 0  ];
eS = [1/3, 1/3, 1/3];

% --- Build triangular mesh via barycentric sampling ---
% Use a regular grid in (wL, wP) with wS = 1-wL-wP >= 0
n   = opts.nPts;
t   = linspace(0, 1, n);
[wL, wP] = ndgrid(t, t);
wS  = 1 - wL - wP;
inside = wS >= -1e-10;   % include boundary

wL = wL(inside);
wP = wP(inside);
wS = max(wS(inside), 0);

% 2-D positions
xy = wL*pL + wP*pP + wS*pS;   % Nx2

% Eigenvalues at each vertex
E  = wL*eL + wP*eP + wS*eS;   % Nx3

% Evaluate color function at every vertex
nv   = size(E,1);
vals = zeros(nv,1);
for k = 1:nv
    ev      = sort(E(k,:), 'descend');
    vals(k) = colorFcn(ev);
end

% --- Build triangle connectivity (Delaunay on the 2-D points) ---
tri = delaunay(xy(:,1), xy(:,2));

% --- Plot ---
% figure('Color','w');
ax = axes('DataAspectRatio',[1 1 1], 'Visible','off');
hold(ax,'on');

% Gouraud-shaded patch: color interpolated smoothly across each triangle
p = patch('Faces', tri, ...
          'Vertices', xy, ...
          'FaceVertexCData', vals, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

% Triangle outline on top
tx = [pL(1), pP(1), pS(1), pL(1)];
ty = [pL(2), pP(2), pS(2), pL(2)];
plot(ax, tx, ty, 'k--', 'LineWidth', 0.5);

% Colormap / colorbar
colormap(ax, opts.cmap);
if ~isempty(opts.clim), clim(ax, opts.clim); end
if opts.showColorbar,   colorbar(ax);         end

% Corner labels
if opts.showCornerLabels
    d = 0.05;
    text(pL(1), pL(2)+d, 'Linear',    'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontSize',11,'FontWeight','bold');
    text(pP(1)-d, pP(2), 'Planar',    'HorizontalAlignment','right', ...
         'VerticalAlignment','middle', 'FontSize',11,'FontWeight','bold');
    text(pS(1)+d, pS(2), 'Spherical', 'HorizontalAlignment','left', ...
         'VerticalAlignment','middle', 'FontSize',11,'FontWeight','bold');
end

if ~isempty(opts.title)
    title(ax, opts.title, 'FontSize',12);
end

hold(ax,'off');
end