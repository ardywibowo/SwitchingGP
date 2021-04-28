function make_plot_figure_2_middle( exp_root, exp_name, params, Ns )
%MAKE_PLOT_FIGURE_2_MIDDLE Produces a plot for Figure 2 (middle)
% after the respective data are sampled.
%
%
% Copyright 2016, Anastasia Podosinnikova



  N = Ns(2); irun = 1; % choosing a single dataset
  
  data_folder = get_folder( exp_name, exp_root );
  data_file = get_file( N, irun );
  data_name = strcat( data_folder, '/', data_file );
  rr = load( data_name, 'X1' ); 
  X1 = rr.X1; 
  
  ff=figure; hold on

    set(gca, 'Fontname', 'Times New Roman')
    set(gca, 'FontSize', 25)
    xlim([0,1000])
    ylim([0,1000])
    set(gcf, 'Color', 'None') % set transparent background

    lblue = [.5 .5 1];
    plot( X1(1,:), X1(2,:), '-', 'Color', 'w', 'Marker', '+', ...
      'MarkerEdgeColor', lblue,  ...
      'MarkerSize', 7, 'LineWidth', 2)
    
    D1 = 1500 * params.D1;
    plot([0; D1(1)], [0; D1(2)],'-','LineWidth',2, 'Color', 'k')

    F1 = 1000 * params.F1;
    for k=1:params.K1
        plot([0; F1(1,k)],[0; F1(2,k)], ':', 'LineWidth', 2, 'Color', 'k')
    end
    ll = legend( '$X_1$', '$D_1$', '$F_{11}$', '$F_{12}$' );
    set(ll,'interpreter','latex');
    set(ll,'Fontname', 'TimesNewRoman')
    set(ff, 'Position', [300, 300, 800, 600]);
    box on
    
    
    % save the plot
    plots_folder = strcat( exp_root, 'plots/' ); 
    if ~(exist( plots_folder, 'dir') == 7)
      mkdir( plots_folder );
    end
    plotname2save = [ plots_folder, 'Figure_2_middle' ];
    saveas(ff,[plotname2save,'.eps'],'epsc')
    disp(['The plot was saved as ', plotname2save, '.eps'])
    
  hold off
  
  
end
