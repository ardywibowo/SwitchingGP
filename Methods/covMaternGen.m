function [A, B] = covMaternGen(logtheta, x, z)
%COVMATERNGEN Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0, A = '2'; return; end          % report number of parameters

persistent K;
range = exp(logtheta(1));
nu = exp(logtheta(2));

dist_x = cell(num_contexts, 1);
for j = 1 : num_contexts
	dist_x{j} = {};
	for i = 1 : length(x{j})
		dist_x{j}{i} = squareform(pdist(x{j}{i}'));
	end
end

dist_z = cell(num_contexts, 1);
for j = 1 : num_contexts
	dist_z{j} = {};
	for i = 1 : length(z{j})
		dist_z{j}{i} = squareform(pdist(z{j}{i}'));
	end
end


num_samples = length(dist_x);
cov_m = cell(num_samples, 1);
for i = 1 : num_samples
	cov_m{i} = matern(dist_x{i}, range, nu);
end
cov_mat = spdiags(blkdiag(cov_m{:}));

if nargin == 2
	K = 
	A = K;
elseif nargout == 2                              % compute test set covariances
	A = sf2*ones(size(z,1),1);
	B = sf2 * exp( -sq_dist( diag(1./ell)*x', diag(1./ell)*z') / 2);
else                                                % compute derivative matrix
	if z <= D                                           % length scale parameters
		A = K.*sq_dist(x(:,z)'/ell(z));  
	else                                                    % magnitude parameter
		A = 2*K;
		clear K;
	end
end

end