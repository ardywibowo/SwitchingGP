function [logtheta_all, nl] = learn_mtgpnew(logtheta_all, deriv_range, data, niter)
%LEARN_MTGP Learns hyperparameters of mtgp model using minimize
%
% INPUT:
% - logtheta_all : All initial hyper-parameter values 
% - deriv_range  : Indices of hyper-parameters to optimize for
% - data         : cell data in the order 
%                  [covfunc_x, xtrain, ytrain, M, irank, nx, ind_kf_train, ind_kx_train]
%
% OUTPUT:
% - logtheta_all : all learned hyperparameters
% - nl           : Final negative marginal likelihood
%
% Edwin V. Bonilla

%% This can be changed: See minimize function for details

% ************* Learning Hyperparameters here *******************************
logtheta0 = logtheta_all(deriv_range); % selects hyper-parameters to optimize

[cov_func, time, signals, num_features, rank_approx, ~, task_index, time_index] = deal(data{:});

[logtheta, nl] = minimize(logtheta0,'nmargl_mtgpnew',niter, logtheta_all, ...
			 cov_func, time, signals, num_features, rank_approx, task_index, time_index, deriv_range);

%% Update whole vector of parameters with learned ones
logtheta_all(deriv_range) = logtheta;

end