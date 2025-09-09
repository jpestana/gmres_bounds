function [radius_vec,centre_vec] = enclosing_circles_multiple_components(A,fname,info,ax_all)

% [radius_vec,centre_vec] = enclosing_circles_multiple_compoenents(A,levels,fname,info,ax_all)

% Function for computing the smallest circle enclosing pseuspectra contours
%
% Inputs:   A:      Matrix for which pseudospectra are computed
%           fname:  Name for saving plots
%           info:   Structure containing information for eigtool plots. 
%                   For each set of contour levels, there should be an 
%                   array of eigtool plot axes (one four-element vector 
%                   per plot) and an array of integers specifing the 
%                   number of points to be used by eigtool.
%           ax_all: Full axis for final contour plots
%
% Outputs:  radius_vec: 1 x n_lev vector of circle radii for each contour
%                       level
%           centre_vec: 2 x n vector of x- and y- coordinates of circle for
%           each contour level
% Requires: https://github.com/AntonSemechko/Bounding-Spheres-And-Circles
%           https://github.com/eigtool/eigtool
%
%
% R. Abu-Labdeh and J. Pestana 5 September 2025

addpath(genpath('./github_repo')) % add Bounding-Spheres-And-Circles to path
addpath(genpath('./eigtool-master')); % add eigtool to path

% Eigtool options
eigtoolopts.direct = 1;

% Number of sets of contour levels
n_sets = length(info);

% Vector for storing all contour levels
all_levels = [];


%%% Eigtool plots %%%

% Loop over sets of contour levels
for k = 1:n_sets

    % Set contour levels for eigtool plots and add to all_levels
    levels_k = info(k).levels;
    eigtoolopts.levels = levels_k;
    all_levels = [all_levels levels_k];

    % Pull out array of axes and ascertain number of eigtool plots needed
    axes_k = info(k).axes;
    n_comp = size(axes_k,1);

    % Pull out array of number of points
    npts_k = info(k).npts;

    % For each eigtool plot
    for j = 1:n_comp

        % Set eigtool options for axis and number of points
        eigtoolopts.ax = axes_k(j,:);
        eigtoolopts.npts = npts_k(j);

        % Run eigtool and save contour matrices
        [eX{j,k},eY{j,k},eZ{j,k}] = eigtool(A,eigtoolopts);
    end
end

% Get total number of contour levels
n_levels = length(all_levels);

% Get eigenvalues within ax_all
d = eig(A);
d(real(d) < ax_all(1) | real(d) > ax_all(2)) = NaN;
d(imag(d) < ax_all(3) | imag(d) > ax_all(4)) = NaN;


%%% Plot pseudospectral contours and eigenvalues %%%
figure; hold on;

% Plot eigenvalues
plot(real(d),imag(d),'.k','Markersize',8);

% Plot pseudospectra
for k = 1:n_sets
    % Get contour levels
    levels_k = info(k).levels; 
    
    % Get number of plots
    axes_k = info(k).axes;
    n_comp = size(axes_k,1);

    % Plot pseudospectra
    for j = 1:n_comp
        [~,ncc_p] = contour(eX{j,k},eY{j,k},log10(eZ{j,k}),levels_k);
        ncc_p.LineWidth = 2;
    end
end

% Plot settings
axis(ax_all);
contour_settings;
colormap("winter");
colorbar('TickLabelInterpreter','latex');
axis equal;
hold off;
saveas(gcf,['./fig/',fname,'_pseudo.png'])

%%% Get bounding circles %%%

% Vectors for radius and centre of each circle
radius_vec = zeros(1,n_levels);
centre_vec = zeros(2,n_levels);

% Plot pseudospectra and bounding circles
f_handle = figure;
hold on;

% Loop over contour levels
for p = 1:n_levels
    % Find which element of info structure contains the given level
    for k = 1:n_sets
        ind = find(info(k).levels==all_levels(p));
        if ~isempty(ind)
            break
        end
    end

    % Contour value
    v = 10.^all_levels(p);

    % Array to hold points on contour boundary
    C = [];

    % Get number of eigtool plots involving given contour level
    axes_k = info(k).axes;
    n_comp = size(axes_k,1);


    for j = 1:n_comp
        % Add contour to overall plot
        [C_comp,ncc_comp] = contour(eX{j,k},eY{j,k},eZ{j,k},[v v]);
        ncc_comp.LineWidth=2;
        ncc_comp.EdgeColor = [0,0,0];

        % Store contour boundary points
        C = [C C_comp];
    end

    % Get the convex hull of the contour boundaries
    Cix = find(C(1,:)==v);
    Cxy = C; Cxy(:,Cix) = [];
    ind = convhull(Cxy');
    Cxy = Cxy(:,ind);

    % Find the bounding circle and store its radius and centre
    [radius,centre] = ExactMinBoundCircle(Cxy');
    radius_vec(p) = radius;
    centre_vec(:,p) = centre';

    % Plot the circle
    xs=linspace(0,2*pi,100);
    circle = centre(1) + radius*exp(1i*xs);
    figure(f_handle);
    plot(circle,'Linewidth',2,'color',[0,0.6,0.6],'LineStyle','--');
end

% Plot settings
axis(ax_all);
axis equal
contour_settings
saveas(gcf,['./fig/',fname,'_contours.png']);
hold off;
