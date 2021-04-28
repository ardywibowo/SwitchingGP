function [logAlpha, logLikelihood] = HMMforwardMTGP(X, Ujs, pSGivenSm, pS1, A, Bjs, C, sigma2, Tskip)
%HMMFORWARDSAR Switching Autoregressive HMM with switches updated only every Tskip timesteps
% [logAlpha, logLikelihood] = HMMforwardInputSAR(X, Ujs, pSGivenSm, pSt1, A, Bjs, C, sigma2, Tskip)
% 
% Inputs:
% X             : State observations
% Ujs           : Input observations
% pSGivenSm     : State (switch) transition matrix p(s(t)|s(t-1))
% pS1           : Prior state distribution
% A             : Matrix of AR coefficients for X. Column A(:,i) are the AR coeffs for switch state i.
% Bjs           : Cell vector of input AR coefficients. Each cell is one input variable. 
%                 Each column Bj(:,i) are the AR coeffs for switch input i.
% C             : Column vector of constant coefficients
% The AR coefficients are in reverse order X(t-L) : X(t-1)
% Each column represents a different switching state
% sigma2        : The innovation noise
% Tskip         : The number of timesteps to skip before a switch update is allowed
%
% Outputs:
% logAlpha      : Log forward messages
% logLikelihood : Sequence log likelihood log p(X(1:T))

T = length(X);
S = size(A, 2);

% Logalpha recursion:
% Mean zero initial
X1 = X(1);

logAlpha(:,1) = -0.5 * repmat(X1.^2, S, 1) ./ sigma2(:) - 0.5 * log(2 * pi * sigma2(:)) + log(pS1);
for t = 2 : T
	[m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
	
	XMinusMean = repmat(X(t), S, 1) - m;
	pHatXGivenS = exp(-0.5 * XMinusMean.^2 ./ sigma2(:)) ./ sqrt(2 * pi * sigma2(:)) + eps;
	
	if mod(t, Tskip) == 0 || Tskip == 0 % Only make a transition every Tskip timesteps
		pStGivenStm = pSGivenSm;
	else
		pStGivenStm = eye(S);
	end
	
	logAlpha(:, t) = brml.logsumexp(repmat(logAlpha(:, t-1), 1, S), repmat(pHatXGivenS', S, 1) .* pStGivenStm');
end

logLikelihood = brml.logsumexp(logAlpha(:, T), ones(S, 1)); % Log likelihood

end