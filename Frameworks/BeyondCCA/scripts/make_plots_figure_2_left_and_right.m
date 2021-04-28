function make_plots_figure_2_left_and_right( ...
  exp_name, exp_root, Ns, nruns, opts, algs )
%MAKE_PLOTS_FIGURE_2_LEFT_AND_RIGHT Produces plots for Figure 2
% (left and right) after the respective data are sampled and respective
% algorithms are evaluated.
%
%
% Copyright 2016, Anastasia Podosinnikova

  data_folder = get_folder( exp_name, exp_root );

  len = length(Ns);
  
  errors_dcca  = zeros(len, nruns); 
  errors_dica  = zeros(len, nruns);
  errors_gdcca = zeros(len, nruns);
  errors_nmf   = zeros(len, nruns);
  errors_nmf0  = zeros(len, nruns);
  errors_gdica = zeros(len, nruns);
  
  % load the results of the experiment
  for in = 1:len
    N = Ns(in);
    for irun = 1:nruns
      data_file = get_file(N, irun);
 
      if sum( strcmp( algs, 'dcca' ) ) > 0
        rr = load( strcat( ...
          data_folder, '/DCCA_', data_file ), ...
          'err1' );
        errors_dcca(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'gdcca' ) ) > 0
        rr = load( strcat( ...
          data_folder, '/gDCCA_', data_file ), ...
          'err1' );
        errors_gdcca(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'dica' ) ) > 0
        Kdica = opts.Kdica;
        rr = load( strcat( ...
          data_folder, '/DICA_Kdica', num2str(Kdica), '_', data_file ), ...
          'err1' );
        errors_dica(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'gdica' ) ) > 0
        Kdica = opts.Kdica;
        rr = load( strcat( ...
          data_folder, '/gDICA_Kdica', num2str(Kdica), '_', data_file ), ...
          'err1' );
        errors_gdica(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'nmf' ) ) > 0
        rr = load( strcat( ...
          data_folder, '/NMF_', data_file ), ...
          'err1' );
        errors_nmf(in,irun) = rr.err1;
      end

      if sum( strcmp( algs, 'nmf0' ) ) > 0
        rr = load( ...
          strcat( data_folder, '/NMF0_', data_file ), ...
          'err1' );
        errors_nmf0(in,irun) = rr.err1;
      end

    end
  end
  
  % compute l1errors
  nalgs = length(algs);
  l1errors = cell(len,1);
  for i = 1:len
    temp = zeros(nruns, nalgs);
    ind = 1;
    if sum( strcmp( algs, 'nmf' ) ) > 0
      temp(:,ind) = errors_nmf(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'nmf0' ) ) > 0
      temp(:,ind) = errors_nmf0(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'dica' ) ) > 0
      temp(:,ind) = errors_dica(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'gdica' ) ) > 0
      temp(:,ind) = errors_gdica(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'dcca' ) ) > 0
      temp(:,ind) = errors_dcca(i,:);
      ind = ind + 1;
    end
    if sum( strcmp( algs, 'gdcca' ) ) > 0
      temp(:,ind) = errors_gdcca(i,:);
    end
    l1errors{i} = temp; 
  end
  save( strcat( data_folder, '/l1errors' ), 'l1errors' )
  
  
  
  
  % setting options for plots
  xs = Ns;
  len_x = length(xs);
  
  algopts.nalg = numel( l1errors );
  algopts.algspecs = {'-', ':','-','-',':'};
  
  orange = [.8 .3 0];
  lorange = [.9  .5 0];
  pink = [1 0 .5];
  lpink = [1 .5 .5];
  lblue = [.5 .5 1];
  blue = [0 .2  .6];
  dcian = [0 .5 .5];
  dviolet = [.5 0 .5];
  green =  [.1 .4 .1];
  lgreen = [.2 .5 .2];
  
  algopts.algcolors = { dcian,  orange, dviolet, green,  blue };
  algopts.markcols  = { lpink, lorange, pink, lgreen,   lblue };
  algopts.markers = {'^', 'v', '*', 'o', 'd'};
  algopts.marksizes = {17, 17, 17, 17, 17};
  algopts.algnames = { 'NMF', 'NMF\circ', 'DICA', 'DCCA', 'DCCAg' }; 
  nalg = algopts.nalg;
  
  plotopts.position = [0 0 1280 800];
  plotopts.fontname = 'Times New Roman';
  plotopts.fontsize = 25;
  plotopts.islatex = 1;
  plotopts.istitle = 0;
  plotopts.titlestr = '';
  
  plotopts.yname = '$\ell_1$-error';
  plotopts.isylim = 1;
  plotopts.ylims = [0 .5];
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
  temp = '_undefined';
  if opts.Kdica == 4,  temp = '_left'; end
  if opts.Kdica == 40, temp = '_right'; end
  plotopts.plotname2save = [ plots_folder, 'Figure_2', temp ];
  [ys, ys_L, ys_U] = convert_to_errorbar_format( l1errors, len_x, nalg );
  make_single_plot(xs, ys, ys_L, ys_U, algopts, plotopts)
  
  disp(['The plot was saved as ', plotopts.plotname2save, '.eps' ])
  
end
