function [rad,cent] = pseudo_bound_multiple_components(A,fname,type,rv,tol,info,ax_all)
%
% [rad,cent] = pseudo_bound_multiple_components(A,fname,type,rv,tol,info,ax_all)
%
% Function for plotting standard pseudospectral bound
%
% Inputs:   A:          Matrix for which bounds are computed
%           fname:      Name for saving plot
%           type:       String to describe type of computation, e.g.,
%                       'all_eigs'
%           rv:         GMRES residual vector
%           tol:        GMRES tolerance
%           info:       Structure for use with enclosing_circles code
%           ax_all:     Axis limits for final contour plots
%
% Outputs:  rad:        Vector of radii of circles enclosing pseudospectra
%           cent:       Vector of centres of circles enclosing pseudospectra
%
% R. Abu-Labdeh and J. Pestana 5 September 2025


% Get vector of all pseudospecra levels
n_blocks = length(info);
levels = [];
for k = 1:n_blocks
    levels = [levels info(k).levels];
end

% Get pseudospectra enclosing circle information
[rad,cent] = enclosing_circles_multiple_components(A,[fname,'_',type],info,ax_all); 

% Set up for bounds
eps = 10.^levels;
epsnum = length(eps);


% Compute and plot convergence bounds
last_it = length(rv)-1;
figure;
hold on
semilogy(0:length(rv)-1,rv/rv(1),'-k','LineWidth',2)

bound=zeros(1,last_it+1);
for j=1:epsnum
    for k = 0: last_it
        bound(1,k+1) = (rad(1,j)/eps(1,j))*(rad(1,j)/cent(1,j))^k;
    end
    iter_vec=0:last_it;
    hold on
    plot(iter_vec,bound,'-','Color',[0.6,0.6,0.6]);
end
yline(tol,'r--');
plot_settings
saveas(gcf,['./fig/',fname,'_',type,'.png'])