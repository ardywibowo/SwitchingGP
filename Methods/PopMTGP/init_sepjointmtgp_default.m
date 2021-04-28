function [log_theta, deriv_range] = init_sepjointmtgp_default(signals_train, cov_func, num_features, rank_approx)
% Initializes parameters of mtgp by default
% You should pay careful attention to modifying this as it may be a 
% bad initialization for your problem!
% 
% INPUT:
% - xtrain: Input training data
% - covfunc_x: Input covariance function
% - M: Number of tasks
% - irank: Rank required for Kf
%
% OUTPUT:
% - logtheta_all: Vector of all hyper-paramters
% - deriv_range: Indices of hyper-parameters to optimize for
%
% Edwin V. BOnilla

num_contexts = length(signals_train);

D = 1;
% Baseline
L_baseline = eval(feval(cov_func{1}));
theta_kx_baseline = log(ones(L_baseline, 1) + 0.01 * randn(L_baseline, 1));

% Contexts
L_context = eval(feval(cov_func{2}));
theta_kx_context = log(ones(L_context * num_contexts, 1) + 0.01 * randn(L_context * num_contexts, 1));

theta_lf0 = [];
for j = 1 : num_contexts
	theta_lf0 = [theta_lf0; init_Kf(num_features, rank_approx)]; %#ok
end

theta_sigma0 = init_sigma(num_features);

log_theta = [theta_lf0; theta_kx_baseline; theta_kx_context; theta_sigma0];
deriv_range = 1 : length(log_theta);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_kx0 = init_kx(xtrain, covfunc_x) %#ok
D         = size(xtrain,2);
L         = eval(feval(covfunc_x{:}));
theta_kx0 = log(ones(L,1) + 0.01 * randn(L, 1)); 
 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_lf0 = init_Kf(M,irank)
Kf0 = eye(M);                 % Init to diagonal matrix (No task correlations)
Lf0 = chol(Kf0)';
theta_lf0 = lowtri2vec_inchol(Lf0,M,irank);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_sigma0 = init_sigma(M) % noise variances
theta_sigma0 =  (1e-7)*rand(M,1);  
 
return;



