% test_jordan
%
% Script for getting GMRES bounds for the two Jordan examples
%
% R. Abu-Labdeh and J. Pestana 5 September 2025

addpath(genpath('./util'));
rng("default")
close all;

% Example settings
ex = 3;
n = 100;
I = eye(n);

% Choose between examples
switch ex
    case 1 
        % Single outlier at 3, Jordan block of dimension 10, other eigenvalues
        % in [0.9,1.1]. Eigenvector matrix has perturbation of 1e-2 affects outlier.

        fname = 'Jordan_1';
        
        % Set up Jordan matrix
        n_out = 1;
        n_jor = 10;
        lambda_0 = 3;
        J = eye(n_jor) + diag(ones(n_jor-1,1),1);
        D =  diag(linspace(0.9,1.1,n-n_jor-n_out));
        
        % Set up eigenvectors
        V = I;
        V(1:2,1) = [1e-2,1];

        % Get coefficient matrix
        A = V*blkdiag(lambda_0,J,D)*inv(V);

        % Pseudospectra information: all eigenalues
        
        % levels 
        levels_all = {-3:-1:-6 -7:-1:-10}; 

        % full axis for plotting
        ax_all = [0.2 lambda_0+.15 -1 1];

        % Set of axes for individual eigtool plots
        ax_a1 = [0.2 1.8 -1 1; lambda_0-1.1e-1 lambda_0+1e-1 -1e-1 1e-1];
        ax_a2 = [0.2 1.8 -1 1; lambda_0-1.1e-5 lambda_0+1e-5 -1e-5 1e-5];
        ax_set_all = {ax_a1 ax_a2};
        npts_all = {[300 1500] [300 1500]};
        
        % Pseudospectra information: cluster
        levels_clust = -3:-1:-10;
        ax_clust = [0.2 1.8 -1 1];
        npts_clust = 300;

    case 3

        % Single outlier at 3, Jordan block of dimension 10, other eigenvalues
        % in [0.9,1.1]. Nonnormality of eigenvectors doesn't affect outlier.

        fname = 'Jordan_3';

        % Set up Jordan matrix
        n_out = 1;
        n_jor = 10;
        lambda_0 = 3;
        J = eye(n_jor) + diag(ones(n_jor-1,1),1);
        D = diag(linspace(0.9,1.1,n-n_jor-n_out)); %diag(0.9+0.2*rand(n-n_jor-n_out,1));

        % Set up eigenvectors
        V = I;
        V(n-1:n,n-1) = [1e-2,1];

        % Get coefficient matrix
        A = V*blkdiag(lambda_0,J,D)*inv(V);

        % Pseudospectra information: all eigenalues
        
        % levels 
        levels_all = {-3:-1:-6 -7:-1:-10}; 

        % full axis for plotting
        ax_all = [0.2 3.1 -1 1];

        ax_a1 = [0.2 1.8 -1 1;lambda_0-1e-3 lambda_0+1e-3 -1e-3 1e-3];
        ax_a2 = [0.2 1.8 -1 1;lambda_0-1e-7 lambda_0+1e-7 -1e-7 1e-7];
        ax_set_all = {ax_a1 ax_a2};
        npts_all = {[300 1500] [300 1500]};
        
        % Pseudospectra information: cluster
        levels_clust = -3:-1:-10;
        ax_clust = [0.2 1.8 -1 1];
        npts_clust = 300;

end

% Right-hand side
b = randn(n,1);

% GMRES
fprintf('Run GMRES\n\n');
tol = 1e-10;
maxdim = 50;
[~,~,~,~,rvg] = rpgmres(A,b,tol,maxdim,1,I,0*b);


% Convergence curves for all eigs
fprintf('Convergence curves for all eigs\n\n')

info_all = struct('levels',levels_all,'npts',npts_all,'axes',ax_set_all);
[r_all,c_all] = pseudo_bound_multiple_components(A,fname,'all_eigs',rvg,tol,info_all,ax_all);


% Convergence curves dropping outlier
fprintf('Convergence curves for dropped outlier\n\n')

info_clust = struct('levels',levels_clust,'npts',npts_clust,'axes',ax_clust);
[r_clust,c_clust] = pseudo_bound_multiple_components(A,fname,'drop_outlier',rvg,tol,info_clust,ax_clust);


% Bounds
fprintf('Bounds\n\n')
cl = linear_bound(A,fname,rvg,tol,levels_clust,ax_clust,r_clust,c_clust,1);
cq = quad_bound(A,fname,rvg,tol,levels_clust,ax_clust,r_clust,c_clust,1);

save(fname,'cl','cq','levels_all','ax_all','r_all','c_all','levels_clust','ax_clust','r_clust','c_clust');