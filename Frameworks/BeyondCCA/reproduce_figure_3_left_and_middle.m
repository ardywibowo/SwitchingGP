%REPRODUCE_FIGURE_2_LEFT_AND_MIDDLE Reproduces Figure 3 (left and middle) 
% from the Appendix (Supplementary material) of the paper:
%   A. Podosinnikova, F. Bach, S. Lacoste-Julien. Beyond CCA: 
%   Moment matching for multi-view models.
%   In Proc. ICML, 2016.
%
% This is a synthetic experiment. The sampled data are "non-Guassian CCA" 
% type. See Appendix G (Supplementary experiments) section from the 
% supplementary material of the paper for more details.
%
% Unless exp_root is changed, the data is saved to the experiments folder.
%
%
% Copyright 2016, Anastasia Podosinnikova

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEFT
clear

exp_name = 'experiment_figure_3_left';
exp_root = 'experiments/';
exp_name = create_exp_folder(exp_name, exp_root);
disp(['Data directory: "',pwd,'/',get_folder(exp_name, exp_root),'"'])

% Initialize the parameters for sampling data
params = sample_parameters_for_figure_3_left;
save( strcat( exp_root, exp_name, '/params' ), 'params' );



% Sample data
Ns = [100, 300, 600, 1000, 2000];
nruns = 5;
for N = Ns
  for irun = 1:nruns
    [X1, X2, ~,~,~] = sample_from_linear_ncca(params, N);
    data_folder = get_folder(exp_name, exp_root);
    data_file = get_file(N, irun);
    save( strcat(data_folder,'/',data_file), 'X1','X2' );
  end
end


% Run algoirthms
scale = 0.001;
K = params.K;
opts = struct( 'K', K, 'Kjade', K, 'scale', scale );
algs = { 'cca', 'jade', 'ncca', 'ncca_spec' };
for N = Ns
  for irun = 1:nruns
    run_experiment_figure_3_left_and_middle(...
      exp_name, exp_root, N, irun, opts, algs )
  end
end

make_plots_figure_3_left_and_middle(...
  exp_name, exp_root, Ns, nruns, opts, algs )


disp(['Data was saved to: "',pwd,'/',get_folder(exp_name, exp_root),'"'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIDDLE
clear

exp_name = 'experiment_figure_3_middle';
exp_root = 'experiments/';
exp_name = create_exp_folder(exp_name, exp_root);
disp(['Data directory: "',pwd,'/',get_folder(exp_name, exp_root),'"'])

% Initialize the parameters for sampling data
params = sample_parameters_for_figure_3_middle;
save( strcat( exp_root, exp_name, '/params' ), 'params' );



% Sample data
Ns = [100, 300, 600, 1000, 2000];
nruns = 5;
for N = Ns
  for irun = 1:nruns
    [X1, X2, ~,~,~] = sample_from_linear_ncca(params, N);
    data_folder = get_folder(exp_name, exp_root);
    data_file = get_file(N, irun);
    save( strcat(data_folder,'/',data_file), 'X1','X2' );
  end
end


% Run algoirthms
scale = 0.001;
K = params.K;
opts = struct( 'K', K, 'Kjade', K, 'scale', scale );
algs = { 'cca', 'jade', 'ncca', 'ncca_spec' };
for N = Ns
  for irun = 1:nruns
    run_experiment_figure_3_left_and_middle(...
      exp_name, exp_root, N, irun, opts, algs )
  end
end

make_plots_figure_3_left_and_middle(...
  exp_name, exp_root, Ns, nruns, opts, algs )


disp(['Data was saved to: "',pwd,'/',get_folder(exp_name, exp_root),'"'])


