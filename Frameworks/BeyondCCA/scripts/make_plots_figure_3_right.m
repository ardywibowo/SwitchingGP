function make_plots_figure_3_right( ...
  exp_name, exp_root, Ns, nruns, deltas )
%MAKE_PLOTS_FIGURE_3_RIGHT Produces plots for Figure 3 (right, from the
% supplementary material/Appendix) after the respective data are sampled 
% and respective algorithms are evaluated.
%
%
% Copyright 2016, Anastasia Podosinnikova


  data_folder = get_folder( exp_name, exp_root );

  len = length(Ns);
  
  errors_1   = zeros(len,nruns); 
  errors_2   = zeros(len,nruns);
  errors_3   = zeros(len,nruns);
  errors_4   = zeros(len,nruns);
  errors_5   = zeros(len,nruns);
  errors_def = zeros(len,nruns);
  
  for in = 1:len
    N = Ns(in);
    for irun = 1:nruns
      data_file = get_file(N,irun);

      rr = load( strcat( ...
        data_folder, '/gDCCA_delta_1_', data_file ), 'err1' );
      errors_1(in,irun) = rr.err1;
  
      rr = load( strcat( ...
        data_folder, '/gDCCA_delta_2_', data_file ), 'err1' );
      errors_2(in,irun) = rr.err1;
  
      rr = load( strcat( ...
        data_folder, '/gDCCA_delta_3_', data_file ), 'err1' );
      errors_3(in,irun) = rr.err1;
  
      rr = load( strcat( ...
        data_folder, '/gDCCA_delta_4_', data_file ), 'err1' );
      errors_4(in,irun) = rr.err1;
  
      rr = load( strcat( ...
        data_folder, '/gDCCA_delta_5_', data_file ), 'err1' );
      errors_5(in,irun) = rr.err1;
  
      rr = load( strcat( ...
        data_folder, '/gDCCA_def_delta_', data_file ), 'err1', 'delta' );
      errors_def(in,irun) = rr.err1;
    end
  end
  
  % compute l1errors
  nalgs = 6;
  l1errors = cell(len,1);
  for i = 1:len
    temp1 = zeros(nruns, nalgs);
    temp1(:,1) = errors_1(i,:);
    temp1(:,2) = errors_2(i,:);
    temp1(:,3) = errors_3(i,:);
    temp1(:,4) = errors_4(i,:);
    temp1(:,5) = errors_5(i,:);
    temp1(:,6) = errors_def(i,:);
    l1errors{i} = temp1; 
  end
  save( strcat( data_folder, '/l1errors' ), 'l1errors' )
  
  
  % Set plot options
  xs = Ns;
  len_x = length(xs);
  
  algopts.nalg = nalgs;
  algopts.algspecs = {'-', '--',':','-','-',':'};
  
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
  
  algopts.algcolors = { blue,  green, dviolet, orange,  dcian, 'r' };
  algopts.markcols  = { lblue, lgreen, pink, lorange,   lpink, 'r' };
  algopts.markers = { '*', '*', 'o', 'v', '^', 'd' };
  algopts.marksizes = {17, 17, 17, 17, 17, 17};
  algopts.algnames = { num2str(deltas(1)), num2str(deltas(2)), num2str(deltas(3)), num2str(deltas(4)), num2str(deltas(5)), 'def' }; 
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
  plotopts.plotname2save = [ plots_folder, 'Figure_3_right' ];
  [ys, ys_L, ys_U] = convert_to_errorbar_format( l1errors, len_x, nalg );
  make_single_plot(xs, ys, ys_L, ys_U, algopts, plotopts)
  
  disp(['The plot was saved as ', plotopts.plotname2save, '.eps'])
    
end



