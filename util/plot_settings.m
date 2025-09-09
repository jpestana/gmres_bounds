% plot_settings
%
% GMRES residual plot settings
%
% R. Abu-Labdeh and J. Pestana 5 September 2025


yl = ylim; 
yl(1) = 1e-15;
ylim(yl);

set(gca,'FontSize',20);
xlabel('Iteration','Interpreter','latex')
ylabel('$\|r_k\|_2/\|r_0\|_2$','Interpreter','latex')
set(gca, 'YScale', 'log')
set(gca,'TickLabelInterpreter','latex')  
