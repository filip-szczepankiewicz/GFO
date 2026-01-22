%% Define problem
%hyper pars
s = 8; kappa = 7; 
nDirs    = 32;
Lmax = 8;
NL = @(L)1/6 *(2 + L) *(3 + 10 *L + 4* L^2);
lvals = 0:Lmax;
Vdeg = (2*lvals+1).*util.sobolev(kappa,s,Lmax).^2;
Vdeg = Vdeg.*(mod(lvals,2)==0);

% tensors
B = diag([0 1 2])/2;
D = diag([0.1 0.1 2.8]);
%% Get orientations ("directions")

[t,x,y,z,fval,exitflag,output] = isotropic_SO3(nDirs, "GFOD2", Lmax,Vdeg);

% Optional conversions
Q = cat(2,t,x,y,z); %quaternions

[alpha,beta,gamma] = util.quat_to_euler(t,x,y,z); %Euler angles

rotMats = util.Rzyz(alpha,beta,gamma); %rotation omatrices

% Rotation axis and angle
rotAx = zeros(nDirs, 3); % Initialize rotation axis array
rotAngle = zeros(nDirs, 1); % Initialize rotation angle array
for j = 1:nDirs
    [rotAx(j,:),rotAngle(j)] = util.quat_to_axis_angle(Q(j,:));
    if rotAx(j,3) < 0
        rotAx(j,3) = - rotAx(j,3);
        rotAngle(j) = -rotAngle(j);
    end
end
%% Compute CV

Euler = [alpha,beta,gamma];
useParfor= false;
cvOpts = struct('batch',5e3, 'useParfor', useParfor);

% define rotations over which to rotate D tensor - here a an exact
% quadrature grid up to L = 6
EU = util.euler_so3_grid(6);
U = util.Rzyz(EU(:,1),EU(:,2),EU(:,3));
cvOpts.rots = U;
nRot = size(U,3);
 w    = ones(1,nDirs)/nDirs;
 %CV contains coefficient of variation of powder averages over rotations of
 %D, sPowder each estimate of the powder average (for each "R D R^T").
[CV, ~, sPowder, ~]   = cv_SO3(Euler,w,B,D,nRot,cvOpts); 
sbar   = mean(sPowder); %average powder average