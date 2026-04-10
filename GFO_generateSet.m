function rotMats = GFO_generateSet(nRot, mode, s, kappa, Lmax)
% function rotMats = GFO_generateSet(nRot, mode, s, kappa, Lmax)

if nargin < 2
    mode = 'GFOD2';
end

if nargin < 3
    s = 8;
    kappa = 7;
    Lmax = 8;
end

% If multiple sets are requested, call the function multiple times
if numel(nRot)>1
    for i = 1:numel(nRot)
        rotMats(i).mode = mode;
        rotMats(i).nRot = nRot(i);
        rotMats(i).rotMats = GFO_generateSet(nRot(i), mode, s, kappa, Lmax);
        disp(['Optimized set with n = ' num2str(nRot(i)) ' rotations.'])
    end
    return
end

switch mode
    case {'ESRS2'}
        u = isotropic_S2(nRot);

        rotMats = zeros(3,3,nRot);

        for i = 1:nRot
            rotMats(:,:,i) = util.rotVec2Vec([0 0 1], u(i,:));
        end

    case {'ESRS2f'} % Fast variant of ESRS2
        u = isotropic_S2_fast2(nRot);

        rotMats = zeros(3,3,nRot);

        for i = 1:nRot
            rotMats(:,:,i) = util.rotVec2Vec([0 0 1], u(i,:));
        end

    case {'ESR', 'ESRD2', 'GFO', 'GFOD2'}
        lvals = 0:Lmax;
        Vdeg = (2*lvals+1).*util.sobolev(kappa,s,Lmax).^2;
        Vdeg = Vdeg.*(mod(lvals,2)==0);

        [t,x,y,z] = isotropic_SO3(nRot, mode, Lmax, Vdeg);
        [alpha,beta,gamma] = util.quat_to_euler(t,x,y,z);
        rotMats = util.Rzyz(alpha,beta,gamma);

    otherwise
        error('Mode not recognized!')
end
