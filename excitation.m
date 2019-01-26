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

boundary = @(x,y,v) (20 * cos(16 * pi * y) * (cos(16 * pi * y) > 0));%(5*(sin(4*pi * y)^2 + sin(4*pi*x)^2));
h_boundary = @(x, y, v)  (1);

sigma = struct(...
    'xs', @(x) (0.2 + 0.2 * x(1, :)),...
    'xa', @(x) (0.2 + 0.2 * x(2, :)) , ...
    'xf', @coefficientXF, ...
    'ms', @(x) (2.0 + 0.2 * x(1,:)), ...
    'ma', @(x) (0.4 + 0.2 * x(2, :)), ...
    'eta',@coefficientETA);
gamma = struct('x', 1, 'f', 1);

% opt must have even number of angles in the angular space.
% used for perfect flipping. For mild anisotrpy, a relative small number
% suffices.
opt = struct('anisotropy', 0.5, 'angle', 24, ...
    'nodes', [0 0; 1 0; 1 1; 0 1]', 'minArea', 4e-5,...
    'sigma', sigma, 'boundary', boundary, 'gamma', gamma, 'h', h_boundary, ...
    'tau', 0.01);

%% construct FUMOT instance.
ft = FUMOT(opt);


%% initialize fluorescent coefficient at excitation wavelength.
XF =  rand(1, ft.RTE.nPoint); % row vec.
% XF(ft.edge) = ft.sigma.xf(ft.edge);
%%
[xf,hist] = ft.ExciteSolve(ft.sigma.xf');

ft.RTE.plot(xf')


