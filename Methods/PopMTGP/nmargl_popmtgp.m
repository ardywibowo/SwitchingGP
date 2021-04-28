function [nl, gradnl] = nmargl_popmtgp(logtheta, logtheta_all, cov_func, time, signal,...
				      m, rank_approx, num_obs, task_index, time_index, deriv_range)
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

D = size(time{1}, 2); %#ok

if ischar(cov_func), cov_func = cellstr(cov_func); end % convert to cell if needed

logtheta_all(deriv_range) = logtheta;
ltheta_x = eval(feval(cov_func{:}));     % number of parameters for input covariance

nlf = rank_approx * (2*m - rank_approx +1) / 2;        % number of parameters for Lf
vlf = logtheta_all(1:nlf);             % parameters for Lf

theta_lf = vlf; 
Lf = vec2lowtri_inchol(theta_lf, m, rank_approx);

theta_x = logtheta_all(nlf+1:nlf+ltheta_x);                     % cov_x parameters
sigma2n = exp(2*logtheta_all(nlf+ltheta_x+1:end));              % Noise parameters
Sigma2n = diag(sigma2n);                                % Noise Matrix

num_realizations = length(time);
nl = 0;
for i = 1 : num_realizations
	n = length(signal{i}); 
	
	Kx = feval(cov_func{:}, theta_x, time{i});
	Kf = Lf*Lf';
	K = kron(Kf, Kx);
	K = K + diag(diag(Sigma2n(task_index{i}, task_index{i}))) / num_obs(1); 
	Sigma_noise = MIN_NOISE * eye(n);
	K = K + Sigma_noise;

	L = chol(K)';                        % cholesky factorization of the covariance
	alpha = solve_chol(L',signal{i});

	% negative log-likelihood
	nl = nl + 0.5*signal{i}'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi);
end

if (nargout == 2)                      % If requested, its partial derivatives
  gradnl = zeros(size(logtheta));        % set the size of the derivative vector
	
	for i = 1 : num_realizations
		n = length(signal{i});
		
		Kx = feval(cov_func{:}, theta_x, time{i});
		Kf = Lf*Lf';
		K = kron(Kf, Kx);
		K = K + diag(diag(Sigma2n(task_index{i}, task_index{i}))) / num_obs(1); 
		Sigma_noise = MIN_NOISE * speye(n);
		K = K + Sigma_noise;

		L = chol(K)'; % cholesky factorization of the covariance
		alpha = solve_chol(L', signal{i});
		
		W = L' \ (L\speye(n)) - alpha*alpha'; % precompute for convenience

		count = 1;
		for zz = 1 : length(deriv_range) 
			 z = deriv_range(zz);

			if ( z <= nlf ) % Gradient wrt  Kf
				[o, p] = pos2ind_tri_inchol(z, m, rank_approx); % determines row and column
				J = zeros(m,m);
				J(o,p) = 1;
				Val = J*Lf' + Lf*J';
				dK = Val(task_index{i}, task_index{i}) .* Kx(time_index{i}, time_index{i});
				
			elseif ( z <= (nlf+ltheta_x) ) % Gradient wrt parameters of Kx
				z_x =  z - nlf;
				dKx = feval(cov_func{:}, theta_x, time{i}, z_x);      
				dK = Kf(task_index{i}, task_index{i}) .* dKx(time_index{i}, time_index{i});
				
			elseif ( z >= (nlf + ltheta_x + 1) ) % Gradient wrt Noise variances
				Val = zeros(m,m);
				kk = z - nlf - ltheta_x;
				Val(kk,kk) = 2 * Sigma2n(kk,kk);
				dK = diag(diag(Val(task_index{i}, task_index{i})));

			end % endif z

			gradnl(count) = gradnl(count) + sum(sum(W.*dK, 2), 1) / 2;
			count = count + 1;
		end % end for derivarives
	end
  
end % end if nargout ==2

% disp('Marginal Iterated');

end

