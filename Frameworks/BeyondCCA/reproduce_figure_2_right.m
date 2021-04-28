%REPRODUCE_FIGURE_2_RIGHT Reproduces Figure 2 (right) from the paper:
%   A. Podosinnikova, F. Bach, S. Lacoste-Julien. Beyond CCA: 
%   Moment matching for multi-view models.
%   In Proc. ICML, 2016.
%
% This is a synthetic experiment. The sampled data are "discrete CCA" type,
% i.e.non-negative count data for each view. See the description of
% the experiments in Section 5 (Experiments) of the paper for more details.
%
% Unless exp_root is changed, the data is saved to the experiments folder.
%
%
% Copyright 2016, Anastasia Podosinnikova

clear

exp_name = 'experiment_figure_2_right';
exp_root = 'experiments/';
exp_name = create_exp_folder(exp_name, exp_root);
disp(['Data directory: "',pwd,'/',get_folder(exp_name, exp_root),'"'])


% Initialize the parameters for sampling data
params = sample_parameters_for_figure_2_right;
save( strcat( exp_root, exp_name, '/params' ), 'params' );



% Sample data
Ns = [500 1000 2000 5000 10000];
nruns = 5;
for N = Ns
  for irun = 1:nruns
    
    [X1, X2] = sample_from_linear_dcca(params, N);
    data_folder = get_folder(exp_name, exp_root);
    data_file = get_file(N, irun);
    save( strcat(data_folder,'/',data_file), 'X1','X2' );

  end
end


% Run algoirthms
scale = 0.001;
K = params.K; K1 = params.K1; K2 = params.K2;
M1 = params.M1; M2 = params.M2;
opts = struct('K', K, 'Kdica', min(M1+M2,K+K1+K2), 'Knmf', K, 'scale', scale);
algs = { 'dcca', 'gdcca', 'dica', 'nmf', 'nmf0' };

% NOTE: you may want to parallizer this loop
for N = Ns
  for irun = 1:nruns
    run_experiment_figure_2( data_folder, N, irun, opts, algs )
  end
end


% Make plots
make_plots_figure_2_left_and_right( ...
  exp_name, exp_root, Ns, nruns, opts, algs )


disp(['Data was saved to: "',pwd,'/',get_folder(exp_name, exp_root),'"'])

