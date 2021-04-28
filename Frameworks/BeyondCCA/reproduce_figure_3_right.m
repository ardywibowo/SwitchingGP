%REPRODUCE_FIGURE_3_RIGHT Reproduces Figure 3 (right) from the Appendix 
% (Supplementary material) of the paper:
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

clear

exp_name = 'experiment_figure_3_right';
exp_root = 'experiments/';
exp_name = create_exp_folder(exp_name, exp_root);
disp(['Data directory: "',pwd,'/',get_folder(exp_name, exp_root),'"'])

params = sample_parameters_for_figure_3_right;
save( strcat( exp_root, exp_name, '/params' ), 'params' );



% Sample data
Ns = [ 100 500 1000 2000 3000 ];
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
K = params.K;
deltas = [ 0.001 0.01 0.1 0.5 1 ]; % scales of the processing points
for N = Ns
  for irun = 1:nruns
    run_experiment_figure_3_right( exp_name, exp_root, N, irun, K, deltas )
  end
end


make_plots_figure_3_right( exp_name, exp_root, Ns, nruns, deltas )


    
    