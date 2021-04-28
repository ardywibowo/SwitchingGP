function make_plots_figure_3_left_and_middle(...
  exp_name, exp_root, Ns, nruns, opts, algs )
%MAKE_PLOT_FIGURE_3_LEFT_AND_MIDDLE Produces a plot for Figure 3 (left and
% middle, from Appendix/supplementary material) after the respective data 
% are sampled.
%
%
% Copyright 2016, Anastasia Podosinnikova

  data_folder = get_folder( exp_name, exp_root );

  len = length(Ns);
  
  errors_ncca = zeros(len,nruns); 
  errors_cca  = zeros(len,nruns);
  errors_jade = zeros(len,nruns);
  errors_spec = zeros(len,nruns);
  
  for in = 1:len
    N = Ns(in);
    for irun = 1:nruns
      data_file = get_file(N,irun);

      if sum( strcmp( algs, 'cca' ) ) > 0
        rr = load( strcat(...
          data_folder, '/CCA_', data_file ), 'err1' );
        errors_cca(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'jade' ) ) > 0
        Kjade = opts.Kjade;
        rr = load( strcat( ...
          data_folder, '/JADE_Kjade', num2str(Kjade), '_', data_file ), 'err1' );
        errors_jade(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'ncca' ) ) > 0
        rr = load( strcat( ...
          data_folder, '/gNCCA_', data_file ), 'err1' );
        errors_ncca(in,irun) = rr.err1;
      end
      
      if sum( strcmp( algs, 'ncca_spec' ) ) > 0
        rr = load( strcat( ...
          data_folder, '/gnSpec_', data_file ), 'err1' );
        errors_spec(in,irun) = rr.err1;
      end

    end
  end
  
  % compute l1errors
  nalgs = length(algs);
  l1errors = cell(len,1);
  for i = 1:len
    temp1 = zeros(nruns, nalgs);
    ind = 1;
    if sum( strcmp( algs, 'cca' ) ) > 0
      temp1(:,ind) = errors_cca(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'jade' ) ) > 0
      temp1(:,ind) = errors_jade(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'ncca' ) ) > 0
      temp1(:,ind) = errors_ncca(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'ncca_spec' ) ) > 0
      temp1(:,ind) = errors_spec(i,:);
    end
    l1errors{i} = temp1; 
  end
  save( strcat( data_folder, '/l1errors' ), 'l1errors' )
  
  
  
  % Set plot parameters
  xs = Ns;
  len_x = length(xs);
  
  % make lines/algorithms options
  algopts.nalg = numel( algs );
  algopts.algspecs = {'-', '-','-',':'};
  
  pink = [1 0 .5];
  lpink = [1 .5 .5];
  lblue = [.5 .5 1];
  blue = [0 .2  .6];
  dcian = [0 .5 .5];
  dviolet = [.5 0 .5];
  green =  [.1 .4 .1];
  lgreen = [.2 .5 .2];
  
  algopts.algcolors = { dcian, dviolet, blue,  green };
  algopts.markcols  = { lpink, pink,    lblue, lgreen };
  algopts.markers = {'^', '*', 'd', 'o'};
  algopts.marksizes = {17, 17, 17, 17};
  algopts.algnames = { 'CCA', 'JADE', 'gNCCA', 'Spec' }; 
  nalg = algopts.nalg;
  
  plotopts.position = [0 0 1280 800];
  plotopts.fontname = 'Times New Roman';
  plotopts.fontsize = 25;
  plotopts.islatex = 1;
  plotopts.istitle = 0;
  plotopts.titlestr = '';
  
  plotopts.yname = '$\ell_1$-error';
  plotopts.isylim = 1;
  plotopts.ylims = [0 1];
  plotopts.isyticks = 0;
  plotopts.yticks = 0;
  
  plotopts.linewidth = 3;
  
  plotopts.islegend = 0;
  plotopts.legendsize = 32;
  
  plotopts.ylim = [0 .5];
  plotopts.xname = 'Number of samples';
  plotopts.isxlim = 1;
  delta = xs(end)/20;
  plotopts.xlims = [xs(1)-delta xs(end)+delta];
  plotopts.isxticks = 0;
  plotopts.xticks = [xs(1)-1 xs xs(end)+1];
  
  plotopts.islegend = 1;
  plotopts.legendsize = 22;
  plotopts.legposition = [.62 .67 .3 .1];
  
  
  
  % make the plot
  plots_folder = strcat( exp_root, 'plots/' ); 
  if ~(exist( plots_folder, 'dir') == 7)
    mkdir( plots_folder );
  end
  temp = '_unknown';
  if opts.Kjade == 1, temp = '_left'; end
  if opts.Kjade == 10; temp = '_middle'; end
  plotopts.plotname2save = [ plots_folder, 'Figure_3', temp ];
  [ys, ys_L, ys_U] = convert_to_errorbar_format( l1errors, len_x, nalg );
  make_single_plot(xs, ys, ys_L, ys_U, algopts, plotopts)
  
  disp(['The plot was saved as ', plotopts.plotname2save, '.eps'])
    
end




