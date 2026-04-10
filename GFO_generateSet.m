function rotMats = GFO_generateSet(nRot, mode, Lmax, s, kappa)
% function rotMats = GFO_generateSet(nRot, mode, Lmax, s, kappa)

if nargin < 2
    mode = 'GFOD2';
end

if nargin < 3
    [Lmax, s, kappa] = this_getDefault(mode);
end


switch mode
    case {'ESRS2', 'GFOS2'}
        Vdeg = util.sobolev(kappa,s,Lmax);
        u = isotropic_S2(nRot, mode, Lmax, Vdeg);

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
end


% Function to generate default hyperparameters for optimization
function [Lmax, s, kappa] = this_getDefault(mode)

switch mode
    case {'ESRS2', 'GFOS2'}
        Lmax = 8; s = 5; kappa = 4;

    case {'ESR', 'ESRD2', 'GFO', 'GFOD2'}
        Lmax = 8; s = 8; kappa = 7;

    otherwise
        error('Mode not recognized!')
end

end