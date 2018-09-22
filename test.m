%
% test for FUMOT class.
% 
% The mesh size has to be fine enough to avoid error on 
% inconsistent transport along different lines. A good 
% choice is to set ``minArea`` smaller than 0.001 for the 
% unit square.
%
% Angular space is chosen to resolve the anisotropy. To fully
% resolve g = 0.8, 64 directions would be at least used.
%
% The solver did not use any preconditioner to accelarate. Usually
% a diffusion acceleration will work well in near isotropic case,
% for anisotropic case, more modes should be involved if the scattering
% is not strong enough.

clear;clc;

boundary = @(x,y,v) (1);

sigma = struct(...
    'xs', @(x) (2 + 0.2 * x(1, :)),...
    'xa', @(x) (0.2 + 0.02 * x(2, :)) , ...
    'xf', @(x) (1.0 + 0.5 * x(1, :)));
gamma = struct('x', 0.3, 'f', 0.7);

opt = struct('anisotropy', 0, 'angle', 64, ...
    'nodes', [0 0; 1 0; 1 1; 0 1]', 'minArea', 0.0002,...
    'sigma', sigma, 'boundary', boundary, 'gamma', gamma);

% construct FUMOT instance.
ft = FUMOT(opt);

% initialize fluorescent coefficient at excitation wavelength.
XF =  ones(1, ft.RTE.nPoint); % row vec.

% internal data from excitation stage.

H = ft.ExciteForwardOp(ft.sigma.xf);




for i = 1:100
    [H0, psi0] = ft.ExciteForwardOp(XF);

    disp(norm(H - H0));
    
    XF = H ./(opt.gamma.f * psi0);
    
end

%%
figure(1);
ft.RTE.plot(XF');

figure(2);
ft.RTE.plot(ft.sigma.xf');







