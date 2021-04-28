function negLogEmission = emissionMTGP(signal, time, feature_model, gp_model, sub_num)
%EMISSIONMTGP Summary of this function goes here
%   Detailed explanation goes here

% Constants
num_context = length(gp_model.coeffs);
rank_approx = gp_model.rank_approx;
num_features = rank_approx;

num_obs = ones(num_features, 1);

% Feature model
means = feature_model.means;
pca_coeffs = feature_model.pca_coeffs;

% GP model
gp_coeffs = gp_model.coeffs;
cov_func = gp_model.cov_func;
rank_approx = gp_model.rank_approx;

time = {time(:)};

negLogEmission = zeros(num_context, 1);
for j = 1 : num_context
	% Feature Extraction
	mean_j = means{j};
	signal_features = (signal - mean_j) * pca_coeffs{j};
% 	signal_features = signal;
	
	% Convert data format to fit with GP function
	task_index = kron(1 : size(signal_features, 2), ones(1, size(signal_features, 1)))';
	time_index = repmat(1 : size(signal_features, 1), 1, size(signal_features, 2))';
	
	task_index = {task_index};
	time_index = {time_index};
	
	signal_features = {signal_features(:)};
	
	log_theta = gp_coeffs{j};
	% End GP data conversion
	
	[negLogEmission(j), ~] = nmargl_popmtgp([], log_theta, cov_func, time, signal_features, num_features, rank_approx, num_obs, task_index, time_index, []);
	% Remark: Standing ~ Laying ~ Sitting
end

end

