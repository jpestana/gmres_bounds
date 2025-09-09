% test_CD_1D
%
% Script for getting GMRES bounds for the 1D convection-diffusion
% (Toeplitz) matrix
%
% R. Abu-Labdeh and J. Pestana 5 September 2025
%
% Requirements: smt toolbox (https://bugs.unica.it/~gppe/soft/smt/)

% setup
close all;
addpath(genpath('./util'));
addpath(genpath('./smt'))
fname = 'CD_1D_right_pre';

% Compute matrix and right-hand side
N = 2000;
e = 1/500;
h = 1/(N+1);
I=eye(N,N);

aa = -e/h^2-1/h;
bb = 2*e/h^2 +1/h;
cc = -e/h^2;

A = full(gallery('tridiag',N,aa,bb,cc)/e);
b = ones(N,1); b(N) = b(N) + e/h^2;

% Preconditioner
P = smtcprec('Optimal',A);
AP = full(A)*inv(P);

% GMRES
fprintf('Run GMRES\n\n');
tol = 1e-10;
maxdim = 50;
[xg,flg,rrg,itg,rvg] = rpgmres(A,b,tol,maxdim,1,full(P),0*b); 

% Convergence curves for all eigs
fprintf('Convergence curves for all eigs\n\n')

% Eigtool levels
levels_all = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_all = [0 10 -5 5];

% Eigenvalues to zoom in on in eigtool plots
l1 = 9.100434398925564;
l2 = 0.251364989521556 + 0.147579889047979i;
l3 = 0.251364989521556 - 0.147579889047979i;

% Set of axes for individual eigtool plots
ax_a1 = [0.2 1.02 -0.4 0.4; l1-8e-3 l1+8e-3 -8e-3 8e-3];
ax_a2 = [0.3 1.02 -0.4 0.4; l1-8e-6 l1+8e-6 -8e-6 8e-6; real(l2)-3e-5 real(l2)+3e-5 imag(l2)-3e-5 imag(l2)+3e-5; real(l3)-3e-5 real(l3)+3e-5 imag(l3)-3e-5 imag(l3)+3e-5];
axes_set_all = {ax_a1 ax_a2};
npts_all = {[400 200] [500 200 200 200]};

info_all = struct('levels',levels_all,'npts',npts_all,'axes',axes_set_all);
[r_all,c_all] = pseudo_bound_multiple_components(AP,fname,'all_eigs',rvg,tol,info_all,ax_all);


% Convergence curves dropping outlier
fprintf('Convergence curves for dropped outlier\n\n')

% Eigtool levels
levels_clust = {-3:-1:-5 -6:-1:-8};

% Full axis for plotting
ax_clust = [0.85 1.02 -0.24 0.24];

% Eigenvalues to zoom in on in eigtool plots
l1 = 1.000249996368900;
l2 = 0.899208565607912 + 0.232263056073509i;
l3 = 0.899208565607912 - 0.232263056073509i;

% Set of axes for individual eigtool plots
ax_c1 = [0.85 1.02 -0.24 0.24; l1-1e-2 l1+1e-2 -1e-2 1e-2; real(l2)-1e-2 real(l2)+1e-2 imag(l2)-1e-2 imag(l2)+1e-2; real(l3)-1e-2 real(l3)+1e-2 imag(l3)-1e-2 imag(l3)+1e-2];
ax_c2 = [0.85 1.02 -0.24 0.24; l1-1e-5 l1+1e-5 -1e-5 1e-5; real(l2)-1e-5 real(l2)+1e-5 imag(l2)-1e-5 imag(l2)+1e-5; real(l3)-1e-5 real(l3)+1e-5 imag(l3)-1e-5 imag(l3)+1e-5];
axes_set_clust = 'axes'; value3 = {ax_c1  ax_c2};
npts_clust = 'npts'; value2 = {[500 200 300 300] [500 200 300 300]};

info_all = struct('levels',levels_clust,'npts',npts_all,'axes',axes_set_clust);
[r_clust,c_clust] = pseudo_bound_multiple_components(AP,fname,'drop_outlier',rvg,tol,info_clust,ax_clust);

% Bounds
fprintf('Bounds\n\n')

levels_clust = -3:-1:-8;
[cl,V,W,D] = linear_bound(AP,fname,rvg,tol,levels_clust,ax_clust,r_clust,c_clust,1);

save(fname,"cl","N","e","b","A","P","r_clust","c_clust","ax_clust","levels_clust","rvg","r_all","c_all","ax_all");


% Plot of outlying eigenvalues
d = diag(D);
x_ax = [ax_clust(1) ax_clust(1) ax_clust(2) ax_clust(2)];
y_ax = [ax_clust(3) ax_clust(4) ax_clust(4) ax_clust(3)];
in = inpolygon(real(d),imag(d),x_ax,y_ax);
outlier = find(~in);
n_out = length(outlier);

d_out = d(outlier);
d_in = d(in);

figure; 
hold on
plot(real(d_in),imag(d_in),'.k','Markersize',8);
plot(real(d_out),imag(d_out),'o','Markersize',8,"MarkerFaceColor",[0 0.7 0.7]);
contour_settings;
axis([0 1.2 -0.6 0.6])
hold off
saveas(gcf,['./fig/',fname,'_eig_cluster.png'])