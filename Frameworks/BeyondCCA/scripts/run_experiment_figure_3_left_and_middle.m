function run_experiment_figure_3_left_and_middle( ...
  exp_name, exp_root, N, irun, opts, algs )
%RUN_EXPERIMENT_FIGURE_3_LEFT_AND_MIDDLE

  disp(['RUN FOR N = ', num2str(N), ' and irun = ', num2str(irun)])

  data_folder = get_folder( exp_name, exp_root );

  rr = load( strcat( data_folder, '/params' ), 'params' );
  params = rr.params;
  D1 = params.D1; D2 = params.D2;

  data_file = get_file( N, irun );
  data_name = strcat( data_folder, '/', data_file );
  rr = load( data_name, 'X1', 'X2' ); 
  X1 = rr.X1; X2 = rr.X2;
  
  % (classical) cca by hotelling
  if sum( strcmp( algs, 'cca' ) ) > 0
    warning('off','stats:canoncorr:NotFullRank')
    tt=tic; [A,B] = canoncorr( full(X1)', full(X2)'); time=toc(tt); %#ok
    err1 = l1_error_cont( [A; B], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/CCA_', data_file ), ...
      'A', 'time', 'err1' )
  end
  
  % (classical) ica with JADE (by Cardoso and Souloumiac)
  if sum( strcmp( algs, 'jade' ) ) > 0
    Kjade = opts.Kjade; M1 = params.M1; M2 = params.M2;
    tt=tic; B = jadeR( [X1; X2], Kjade ); time=toc(tt); %#ok
    B = B';
    [D1est, D2est] = get_dcca_matrices_oracle( B, M1, M2 );
    err1 = l1_error_cont( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/JADE_Kjade', num2str(Kjade), '_', data_file ), ...
      'B', 'time', 'err1' )
  end
  
  % ncca
  if sum( strcmp( algs, 'ncca' ) ) > 0
    K = opts.K; scale = opts.scale;
    tt=tic; [D1est, D2est] = ncca( X1, X2, K, scale ); time=toc(tt); %#ok
    err1 = l1_error_cont( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/gNCCA_', data_file ), ...
      'D1est', 'D2est', 'time', 'err1' )
  end
  
  % spectral for ncca
  if sum( strcmp( algs, 'ncca_spec' ) ) > 0
    K = opts.K; scale = opts.scale;
    tt=tic; [D1est, D2est] = ncca_spectral( X1, X2, K, scale ); time=toc(tt); %#ok
    err1 = l1_error_cont( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/gnSpec_', data_file ), ...
      'D1est', 'D2est', 'time', 'err1' )
  end
end
