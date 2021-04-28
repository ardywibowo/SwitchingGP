function [pStGivenX1T, pStStpGivenX1T] = HMMsmoothMTGP(logAlpha, logBeta, A, Bjs, C, sigma2, pSGivenSm, X, Ujs, Tskip)
%HMMSMOOTHSAR Switching Autoregressive HMM smoothing
% [pStGivenX1T, pStStpGivenX1T] = HMMsmoothInputSAR(logAlpha, logBeta, A, Bjs, C, sigma2, pSGivenSm, X, Ujs, Tskip)
% Returns the smoothed pointwise posterior p(s(t)|x(1:T)) and pairwise smoothed posterior p(s(t),s(t+1)|x(1:T)). 
%
% Inputs:
% logAlpha  : Log alpha messages (see HMMforwardSAR.m)
% logBeta   : Log beta messages (see HMMbackwardSAR.m)
% A         : Matrix of AR coefficients for x. Column a(:,i) are the AR coeffs for switch state i. 
% Bjs       : Cell vector of input AR coefficients. Each cell is one input variable. 
%             Each column Bj(:,i) are the AR coeffs for switch input i.
% C         : Column vector of constant coefficients
% The AR coefficients are in reverse order X(t-L) : X(t-1)
% Each column represents a different switching state
% sigma2    : The innovation noise
% pSGivenSm : State (switch) transition matrix p(s(t)|s(t-1))
% X         : Observations
% Ujs       : Inputs. Each j corresponding to a different type.
% Tskip     : The number of timesteps to skip before a switch update is allowed
%
% Outputs:
% pStGivenX1T     : Smoothed posterior p(s(t)|x(1:T))
% pStStpGivenX1T  : Smoothed posterior p(s(t),s(t+1)|x(1:T))
% See also HMMforwardSAR.m, HMMbackwardSAR.m, demoSARlearn.m

import brml.*

T = length(X);
S = size(A, 2);

% Smoothed posteriors: pointwise marginals:
for t = 1 : T
	logPStGivenX1T(:, t) = logAlpha(:, t) + logBeta(:, t); % alpha-beta approach
	pStGivenX1T(:, t) = condexp(logPStGivenX1T(:, t));
end

% Smoothed posteriors: pairwise marginals p(h(t),h(t+1)|v(1:T)):
for t=2:T
	aTemp = condexp(logAlpha(:, t-1));
	bTemp = condexp(logBeta(:, t));
	
	[m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
	XMinusMean = repmat(X(t), S, 1) - m;
	
	logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));
	pHatXGivenS = condexp(logPHatXGivenS);
	if mod(t, Tskip) == 0 || Tskip == 0
		pSGivenSmt = pSGivenSm;
	else
		pSGivenSmt=eye(S);
	end
	cTemp = repmat(aTemp, 1, S) .* pSGivenSmt' .* repmat(pHatXGivenS' .* bTemp', S, 1) + eps; % Two timestep potential
	pStStpGivenX1T(:, :, t-1) = cTemp ./ sum(sum(cTemp));
end

end