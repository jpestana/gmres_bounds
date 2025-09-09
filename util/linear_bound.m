function [cl,V,W,D] = linear_bound(A,fname,rv,tol,levels,ax,r_clust,c_clust,add_iters)

% [cl,V,W,D] = linear_bound(A,fname,rv,tol,levels,ax,r_clust,c_clust,add_iters)
% 
% Function to compute linear bound with one or more outliers
%
% Inputs:   A:          Matrix for which bounds are computed
%           fname:      Name for saving plot
%           rv:         GMRES residual vector
%           tol:        GMRES tolerance
%           levels:     Pseudospectra contour levels
%           ax:         Cluster axis limits
%           r_clust:    Vector of circle radii for circles enclosing
%                       pseudospectra
%           c_clust:    Vector of circle centres for circles enclosing
%                       pseudospectra
%           add_iters:  Any additional iteration steps to compute bounds
%                       for
%
% Outputs:  cl:         Constant for bound
%           V:          Matrix of right eigenvectors
%           D:          Matrix of eigenvalues
%           W:          Matrix of left eigenvectors
%
% R. Abu-Labdeh and J. Pestana 5 September 2025

% Get eigenvalues and eigenvectors
n = size(A,1);
[V,D,W]=eig(full(A)); %V is the right eigenvec, Wis left eigenvecs. Outlier eigval is D(1,1)
d = diag(D);
I = eye(n);

% Find outliers
x_ax = [ax(1) ax(1) ax(2) ax(2)];
y_ax = [ax(3) ax(4) ax(4) ax(3)];
in = inpolygon(real(d),imag(d),x_ax,y_ax);
outlier = find(~in);
n_out = length(outlier);

% Find projectors
P0 = zeros(n);
for j = 1:n_out
    ind = outlier(j);
    P0 = P0 + 1/(W(:,ind)'*V(:,ind)) * V(:,ind)*W(:,ind)';
end

Pc= I-P0; %cluster spectral projector
Ac=Pc*A;

% Get multiplier
fac = I;
for j = 1:length(outlier)
    lam = d(outlier(j));
    disp(norm(Pc-(Ac/lam)));
    fac = fac*(Pc-(Ac/lam));
end
cl = norm(fac);

% Set up for linear bound
eps = 10.^levels;
epsnum = length(eps);
last_it = length(rv)-1;

% Compute and plot bounds
figure;
hold on
semilogy(0:last_it,rv/rv(1),'-k','Linewidth',2);

for j=1:epsnum
    clear lin_bound
    for k = n_out: last_it+add_iters
        lin_bound(1,k-n_out+1)=cl*(r_clust(1,j)/eps(1,j))*(r_clust(1,j)/c_clust(1,j))^(k-n_out);
    end
    iter_vec=n_out:last_it+add_iters;
    plot(iter_vec,lin_bound,'-','Color',[0.6 0.6 0.6]);
end
yline(tol,'r--');
plot_settings;
saveas(gcf,['./fig/',fname,'_linear_bound.png'])