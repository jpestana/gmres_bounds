% test_suitesparse_gre_1107
%
% Script for getting GMRES bounds for the matrix GRE_1107
%
% R. Abu-Labdeh and J. Pestana 5 September 2025


% Setup 
fname = 'gre_1107';
close all;

% Get coefficient matrix
fprintf('Get problem: %s\n\n',fname);
cd ('./SuiteSparse/')
load([fname,'.mat']);
cd('../')

A = Problem.A; 
n = size(A,1);

% Random right-hand side
b = randn(n,1);

% Compute ILU preconditioner
fprintf('Get ILU preconditioner\n\n')
options.droptol = 1e-4;
options.type = "ilutp";
[L,U] = ilu(A,options);

% Compute preconditioned matrix
fprintf('Compute preconditioned matrix\n\n')
AP = full(A/U/L);

% GMRES
fprintf('Run GMRES\n\n');
tol = 1e-10;
maxdim = 50;
[~,~,~,~,rvg] = rpgmres(A,b,tol,maxdim,1,L*U,0*b);


% Convergence curves for all eigs
fprintf('Convergence curves for all eigs\n\n')

% levels
levels_all = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_all = [-0.6 1.5 -0.5 0.5];

% Eigenvalues to zoom in on in eigtool plots
l1 = -0.2443;
l2 = 0.0876;
l3 = 0.5142;
l4 = 0.7844;
l5 = 1.3976;

% Set of axes for individual eigtool plots
ax_a1 = [0.9 1.12 -0.1 0.1; -0.35 0.61 -0.2 0.2;  0.74 0.84 -0.05 0.05; 1.29 1.47 -0.1 0.1];
ax_a2 = [0.96 1.05 -0.05 0.05; l1-5e-4 l1+5e-4 -5e-4 5e-4; l2-5e-4 l2+5e-4 -5e-4 5e-4; l3-5e-4 l3+5e-4 -5e-4 5e-4; l4-5e-4 l4+5e-4 -5e-4 5e-4; l5-5e-4 l5+5e-4 -5e-4 5e-4];
axes_set_all = {ax_a1 ax_a2};
npts_all = {300*ones(1,4) [1000 1000 1000 1000 1000 1000]};

info_all = struct('levels',levels_all,'npts',npts_all,'axes',axes_set_all);
[r_all,c_all] = pseudo_bound_multiple_components(AP,fname,'all_eigs',rvg,tol,info_all,ax_all);

% Convergence curves dropping outlier
fprintf('Convergence curves for dropped outlier\n\n')

% levels
levels_clust = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_clust = [0.9 1.12 -0.1 0.1];

% Eigenvalues to zoom in on in eigtool plots
l1 = 1.01202+0.0418584*1i;
l2 = 1.01202-0.0418584*1i;
l3 = 1.04301;

% Set of axes for individual eigtool plots
ax_c2 = [0.982 1.007 -0.011 0.011; real(l1)-1e-4 real(l1)+1e-4 imag(l1)-1e-4 imag(l1)+1e-4; real(l2)-1e-4 real(l2)+1e-4 imag(l2)-1e-4 imag(l2)+1e-4; real(l3)-1e-4 real(l3)+1e-4 imag(l3)-1e-4 imag(l3)+1e-4];
axes_set_clust = {ax_clust ax_c2};
npts_clust = {50 [1000;500;500;500]};

info_clust = struct('levels',levels_clust,'npts',npts_clust,axes_set_clust,'axes');
[r_clust,c_clust] = pseudo_bound_multiple_components(AP,fname,'drop_outlier',rvg,tol,info_clust,ax_c1);

% Bounds
fprintf('Bounds\n\n')
ax_clust = ax_c1;
[cl,V,W,D] = linear_bound(AP,fname,rvg,tol,levels_clust,ax_clust,r_clust,c_clust,1);

save(fname,"cl","options","r_clust","c_clust","rvg","tol","levels_clust","ax_clust","r_clust")