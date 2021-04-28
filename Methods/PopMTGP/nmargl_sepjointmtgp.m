function [nl, gradnl] = nmargl_sepjointmtgp(logtheta, logtheta_all, cov_func, time, signals,...
				      num_features, rank_approx, num_obs, task_index, time_index, deriv_range)
%NMARGL_POPMTGP Summary of this function goes here
% Marginal likelihood and its gradients for multi-task Gaussian Processes
% 
% [nl, gradnl] = nmargl_mtgp(logtheta, logtheta_all, covfunc_x, x, y,...
%                	      m, irank, nx, ind_kf, ind_kx, deriv_range)
%
% To be used in conjunction with Carl Rasmussen's minimize function
% and the gpml package http://www.gaussianprocess.org/gpml/code/matlab/doc/
%
% nl = nmargl_mtgp(logtheta, ...) Returns the marginal negative log-likelihood
% [nl gradnl] =  nmargl_mtgp(logtheta, ...) Returns the gradients wrt logtheta
%
% logtheta    : Column vector of initial values of parameters to be optimized
% logtheta_all: Vector of all parameters: [theta_lf; theta_x; sigma_l]
%                - theta_lf: the parameter vector of the
%                   cholesky decomposition of k_f
%                - theta_x: the parameters of K^x
%                - sigma_l: The log of the noise std deviations for each task
% covfunc_x   : Name of covariance function on input space x
% x           : Unique input points on all tasks 
% y           : Vector of target values
% m           : The number of tasks
% irank       : The rank of K^f 
% nx          : number of times each element of y has been observed 
%                usually nx(i)=1 unless the corresponding y is an average
% ind_kx      : Vector containing the indexes of the data-points in x
%                which each observation y corresponds to
% ind_kf      : Vector containing the indexes of the task to which
%                each observation y corresponds to
% deriv_range : The indices of the parameters in logtheta_all
%                to which each element in logtheta corresponds to
%

% *** General settings here ****
config = get_mtgp_config();
MIN_NOISE = config.MIN_NOISE;
% ******************************

D = size(time{1}{1}, 2); %#ok
num_contexts = length(signals);

logtheta_all(deriv_range) = logtheta;

% Get parameter lengths
length_baseline = eval(feval(cov_func{1}));
length_contexts = eval(feval(cov_func{2}));
length_Lf = num_contexts * rank_approx * (2*num_features - rank_approx +1) / 2; % number of parameters for Lf

length_Lf_curr = length_Lf / num_contexts;

% Get Parameters
theta_Lf = logtheta_all(1 : length_Lf);

theta_baseline = logtheta_all(length_Lf+1 : length_Lf+length_baseline);
theta_contexts = logtheta_all(length_Lf+length_baseline+1 : length_Lf+length_baseline + num_contexts * length_contexts);

sigma2n = exp(2 * logtheta_all(length_Lf+length_baseline + num_contexts * length_contexts + 1 : end));
Sigma2n = diag(sigma2n);

W_saved = cell(num_contexts, 1);
Kx_saved = cell(num_contexts, 1);
Kf_saved = cell(num_contexts, 1);
Lf_saved = cell(num_contexts, 1);

nl = 0;
for j = 1 : num_contexts
	num_realizations = length(time{j});
	W_saved{j} = cell(num_realizations, 1);
	Kx_saved{j} = cell(num_realizations, 1);
	
	current_theta_Lf = theta_Lf((j-1)*length_Lf_curr +1 : j*length_Lf_curr);
	Lf = vec2lowtri_inchol(current_theta_Lf, num_features, rank_approx);
	Kf = Lf*Lf';
	
	Kf_saved{j} = Kf;
	Lf_saved{j} = Lf;
	
	for i = 1 : num_realizations
		n = length(signals{j}{i});
		
		% Spatial correlation
		Kx_baseline = feval(cov_func{1}, theta_baseline, time{j}{i});
		
		current_theta_context = theta_contexts((j-1)*length_contexts + 1 : j*length_contexts);
		Kx_contexts = feval(cov_func{2}, current_theta_context, time{j}{i});
		
		Kx = Kx_baseline + Kx_contexts;
		Kx_saved{j}{i} = Kx;
		
		% Multivariate correlation
		K = kron(Kf, Kx);
		
		% Likelihood
		K = K + diag(diag(Sigma2n(task_index{j}{i}, task_index{j}{i}))) / num_obs(1); 
		Sigma_noise = MIN_NOISE * eye(n);
		K = K + Sigma_noise;

		L = chol(K)'; % Cholesky factorization of the covariance
		alpha = solve_chol(L',signals{j}{i});
		
		W = L' \ (L\speye(n)) - alpha*alpha'; % precompute for convenience
		W_saved{j}{i} = W;
		
		% Negative log-likelihood
		nl = nl + 0.5*signals{j}{i}'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi);
	end
end

if (nargout == 2) % If requested, its partial derivatives
  gradnl = zeros(size(logtheta)); % set the size of the derivative vector
	
	count = 1;
	for zz = 1 : length(deriv_range)
		z = deriv_range(zz);
		if (z <= length_Lf) % Gradient wrt  Kf			
			j = floor((z-1)/length_Lf_curr) + 1;
			num_realizations = length(time{j});
			z_index = mod(z-1, length_Lf_curr) + 1;
			
			Lf = Lf_saved{j};
			for i = 1 : num_realizations					
				Kx = Kx_saved{j}{i};
				W = W_saved{j}{i};

				% Gradient calculations
				[o, p] = pos2ind_tri_inchol(z_index, num_features, rank_approx); % determines row and column
				J = zeros(num_features, num_features);
				J(o,p) = 1;
				Val = J*Lf' + Lf*J';
				dK_baseline = Val(task_index{j}{i}, task_index{j}{i}) .* Kx(time_index{j}{i}, time_index{j}{i});

				gradnl(count) = gradnl(count) + sum(sum(W.*dK_baseline, 2), 1) / 2;
			end
			count = count + 1;
		elseif (z <= length_Lf + length_baseline) % Gradient wrt baseline parameters
			for j = 1 : num_contexts
				num_realizations = length(time{j});
				Kf = Kf_saved{j};
				for i = 1 : num_realizations
					W = W_saved{j}{i};

					% Gradient calculations
					z_x =  z - length_Lf;
					
					Kx_baseline = feval(cov_func{1}, theta_baseline, time{j}{i}); %#ok
					dKx = feval(cov_func{1}, theta_baseline, time{j}{i}, z_x);
					dK_baseline = Kf(task_index{j}{i}, task_index{j}{i}) .* dKx(time_index{j}{i}, time_index{j}{i});
					
					gradnl(count) = gradnl(count) + sum(sum(W.*dK_baseline, 2), 1) / 2;
				end
			end
			count = count + 1;
		elseif (z <= length_Lf + length_baseline + num_contexts * length_contexts) % Gradient wrt context parameters
			% Pick the subjects from specific context
			z_index = z - (length_Lf + length_baseline);
			j = floor((z_index-1)/length_contexts) + 1;
			
			num_realizations = length(time{j});
			Kf = Kf_saved{j};
			for i = 1 : num_realizations
				W = W_saved{j}{i};

				% Gradient calculations
				z_x =  z_index - (j-1) * length_contexts;
				
				current_theta_context = theta_contexts((j-1)*length_contexts + 1 : j*length_contexts);
				Kx_context = feval(cov_func{2}, current_theta_context, time{j}{i}); %#ok
				dKx = feval(cov_func{2}, current_theta_context, time{j}{i}, z_x);
				dK_context = Kf(task_index{j}{i}, task_index{j}{i}) .* dKx(time_index{j}{i}, time_index{j}{i});

				gradnl(count) = gradnl(count) + sum(sum(W.*dK_context, 2), 1) / 2;
			end
			count = count + 1;
			
		elseif (z >= length_Lf + length_baseline + num_contexts * length_contexts + 1) % Gradient wrt Noise variances
			for j = 1 : num_contexts
				num_realizations = length(time{j});
				for i = 1 : num_realizations
					W = W_saved{j}{i};
					
					% Gradient calculations
					Val = zeros(num_features, num_features);
					kk = z - length_Lf - length_baseline - num_contexts * length_contexts;
					Val(kk,kk) = 2 * Sigma2n(kk,kk);
					dK_noise = diag(diag(Val(task_index{j}{i}, task_index{j}{i})));
					
					gradnl(count) = gradnl(count) + sum(sum(W.*dK_noise, 2), 1) / 2;
				end
			end
			count = count + 1;
		end
	end
end

end
