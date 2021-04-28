function run_experiment_figure_2( data_folder, N, irun, opts, algs )
%RUN_EXPERIMENT_FIGURE_2
  
  disp(['RUN FOR N = ', num2str(N), ' and irun = ', num2str(irun)])

  rr = load( strcat( data_folder, '/params' ), 'params' );
  params = rr.params;
  D1 = params.D1; D2 = params.D2;

  data_file = get_file( N, irun );
  data_name = strcat( data_folder, '/', data_file );
  rr = load( data_name, 'X1', 'X2' ); 
  X1 = rr.X1; X2 = rr.X2; X = [X1; X2];

  if sum( strcmp( algs, 'dica' ) ) > 0
    Kdica = opts.Kdica; M1 = params.M1; M2 = params.M2; K = opts.K;
    tt=tic; Dest = dica(X, Kdica); time=toc(tt); %#ok
    [D1est, D2est] = get_dcca_matrices_oracle( Dest, M1, M2 );
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/DICA_Kdica', num2str(Kdica), '_', data_file ), ...
      'Dest', 'time', 'err1' )
  end

  if sum( strcmp( algs, 'dcca' ) ) > 0
    K = opts.K;
    tt=tic; [D1est, D2est] = dcca(X1, X2, K); time=toc(tt); %#ok
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/DCCA_', data_file ), ...
      'D1est', 'D2est', 'time', 'err1' )
  end

  if sum( strcmp( algs, 'gdcca' ) ) > 0
    K = opts.K; scale = opts.scale;
    tt=tic; [D1est, D2est] = dcca_gencov(X1, X2, K, scale); time=toc(tt); %#ok
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/gDCCA_', data_file ), ...
      'D1est', 'D2est', 'time', 'err1' )
  end

  if sum( strcmp( algs, 'nmf' ) ) > 0
    M1 = params.M1; M2 = params.M2; 
    M = M1 + M2;
    K1 = params.K1; K2 = params.K2; K = opts.K;
    Knmf = K + K1 + K2;
    % Initialization for NMF
    D0 = rand(M, Knmf);
    H0 = rand(Knmf, N);
    tt = tic; Dest = nmf( full(X), D0, H0 ); time=toc(tt); %#ok
    [D1est, D2est] = get_dcca_matrices_oracle( Dest, M1, M2 );
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/NMF_', data_file ), ...
      'Dest', 'time', 'err1' )
  end

  if sum( strcmp( algs, 'nmf0' ) ) > 0
    M1 = params.M1; M2 = params.M2; 
    M = M1 + M2;
    K1 = params.K1; K2 = params.K2; K = opts.K;
    Knmf = K + K1 + K2;
    % Initialization for NMF
    D = [rand(M,K) [ rand(M1,K1);zeros(M2,K1)] [zeros(M1,K2);rand(M2,K2)]];
    D0 = D; H0 = rand(Knmf, N);
    tt = tic; Dest = nmf( full(X), D0, H0 ); time=toc(tt); %#ok
    [D1est, D2est] = get_dcca_matrices_oracle( Dest, M1, M2 );
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/NMF0_', data_file ), ...
      'Dest', 'time', 'err1' )
  end

end
