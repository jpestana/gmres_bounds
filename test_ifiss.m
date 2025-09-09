% test_ifiss
%
% Script for getting GMRES bounds for the IFISS test matrix
%
% R. Abu-Labdeh and J. Pestana 5 September 2025
% 
% Requirements: IFISS
% (https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/)
%
% First run IFISS convection-diffusion test problem 4 (recirculating wind)
% with epsilong = 1/300 and all other parameters with default values

% Setup
addpath(genpath("./util"))

% Preconditioner
here = pwd;
cd_solve_amg_only
cd(here);
fname = 'ifiss'

% Coefficient matrix and RHS
A = Asupg; 
b = fsupg;
n = length(b);
AP = A*P;

% GMRES
fprintf('Run GMRES\n\n');
tol = 1e-10;
maxdim = 50;
[~,~,~,~,rvg] = rpgmres(A,b,tol,maxdim,1,inv(P),0*b);


% Convergence curves for all eigs
fprintf('Convergence curves for all eigs\n\n')

% Eigtool levels
levels_all = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_all = [-0.6 1.5 -0.5 0.5];

% Eigenvalues to zoom in on in eigtool plots
l1 = 0.121633378620697;
l2 = 0.414762823284323;
l3 = 1.232466378067176 - 0.043737879310167i;
l4 = 1.232466378067176 + 0.043737879310167i;

% Set of axes for individual eigtool plots
ax_a1 = [0.6 1.3 -.3 .3; l1-2e-3 l1+2e-3 -2e-3 2e-3; l2-2e-3 l2+2e-3 -2e-3 2e-3; real(l3)-1e-2 real(l3)+1e-2 imag(l3)-1e-2 imag(l3)+1e-2; real(l4)-1e-2 real(l4)+1e-2 imag(l4)-1e-2 imag(l4)+1e-2];
ax_a2 = [0.6 1.3 -.3 .3; l1-2e-6 l1+2e-6 -2e-6 2e-6; l2-2e-6 l2+2e-6 -5e-6 5e-6; real(l3)-1e-5 real(l3)+1e-5 imag(l3)-1e-5 imag(l3)+1e-5; real(l4)-1e-5 real(l4)+1e-5 imag(l4)-1e-5 imag(l4)+1e-5];
axes_set_all = {ax_a1 ax_a2};
npts_all = {[300 200 200 500 500] [1000 500 500 500 500]};

info_all = struct('levels',levels_all,'npts',npts_all,'axes',axes_set_all);
[r_all,c_all] = pseudo_bound_multiple_components(AP,fname,'all_eigs',rvg,tol,info_all,ax_all);


% Convergence curves dropping outlier
fprintf('Convergence curves for dropped outlier\n\n')

% Eigtool levels
levels_clust = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_clust = [0.6 1.3 -0.5 0.5];

% Eigenvalues to zoom in on in eigtool plots
l1 = 0.638538852461685;
l3 = 1.232466378067176 - 0.043737879310167i;
l4 = 1.232466378067176 + 0.043737879310167i;

% Set of axes for individual eigtool plots
ax_c1 = [0.65 1.3 -.3 .3; l1-4e-3 l1+4e-3 -4e-3 4e-3; real(l3)-1e-2 real(l3)+1e-2 imag(l3)-1e-2 imag(l3)+1e-2; real(l4)-1e-2 real(l4)+1e-2 imag(l4)-1e-2 imag(l4)+1e-2];
ax_c2 = [0.65 1.3 -.3 .3; l1-5e-6 l1+5e-6 -5e-6 5e-6; real(l3)-1e-5 real(l3)+1e-5 imag(l3)-1e-5 imag(l3)+1e-5; real(l4)-1e-5 real(l4)+1e-5 imag(l4)-1e-5 imag(l4)+1e-5];
axes_set_clust = {ax_c1 ax_c2};
npt_clust = {[300 500 500 500] [1000 700 500 500]};


info_clust = struct('levels',levels_clust,'npts',npts_clust,'axes',axes_set_clust); 
[r_clust,c_clust] = pseudo_bound_multiple_components(AP,fname,'drop_outlier',rvg,tol,info_clust,ax_clust);


% Bounds
fprintf('Bounds\n\n')
levels_clust = -3:-1:-8;
ax_clust = [0.6 1.3 -.3 .3];
[cl,V,W,D] = linear_bound(AP,fname,rvg,tol,levels_clust,ax_clust,r_clust,c_clust,1);

save(fname,"cl","D", "W", "V", "c_all", "r_all", "c_clust", "r_clust","ax_clust","A","P")