function run_experiment_figure_3_right( ...
  exp_name, exp_root, N, irun, K, deltas )
%RUN_EXPERIMENT_FIGURE_3_RIGHT

  data_folder = get_folder( exp_name, exp_root );

  rr = load( strcat( data_folder, '/params' ), 'params' );
  params = rr.params;
  D1 = params.D1; D2 = params.D2;

  data_file = get_file( N, irun );
  data_name = strcat( data_folder, '/', data_file );
  rr = load( data_name, 'X1', 'X2' ); 
  X1 = rr.X1; X2 = rr.X2;
  
  len = length(deltas);
  
  for rr = 1:len
    delta = deltas(rr);
    tt = tic; [D1est,D2est] = dcca_gencov( X1, X2, K, delta ); time=toc(tt); %#ok
    err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
    save( strcat( ...
      data_folder, '/gDCCA_delta_',num2str(rr),'_', data_file ), ...
      'D1est', 'D2est', 'time', 'err1', 'delta' )
  end
  
  delta = 0.1 / ( mean( mean(X1) + mean(X2) ) );
  % THIS IS THE DEFALT VALUE - THERE IS A TYPO IN THE PAPER IN THE CAPTION
  % OF FIGURE 3
  tt = tic; [D1est,D2est] = dcca_gencov( X1, X2, K, delta ); time=toc(tt); %#ok
  err1 = l1_error( [D1est; D2est], [D1; D2] ); %#ok
  save( strcat( ...
    data_folder, '/gDCCA_def_delta_', data_file ), ...
    'D1est', 'D2est', 'time', 'err1', 'delta' )
  
  
end
