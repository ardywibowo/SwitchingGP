function [SHat, logProbability] = HMMinputViterbi(X, Ujs, THat, A, Bjs, C, pS1, pSGivenSm, sigma2, Tskip)
%HMMINPUTVITERBI Viterbi most likely joint hidden state of a HMM
%   [maxState logProbability] = HMMinputViterbi(X, Ujs, pSGivenSm, pS1, pXGivenS)
%
% Inputs:
% X              : Row vector of observations
% Ujs            : Cell vector. Each cell containing an input row vector.
% A              : AR coefficients for x. Shared between each
%                  subject.
% Bjs            : Cell vector. Each cell contains learned AR coefficients for
%                  different input types Uis, shared between each subject.
% C              : Constant Coefficient
% pS1            : Prior state distribution 
% pSGivenSm      : State (switch) transition matrix p(s(t)|s(t-1))
% sigma2         : The innovation noise
% Tskip          : The number of timesteps to skip before a switch update is allowed
%
% Outputs:
% maxState       : Most likely joint state sequence
% logProbability : Associated log probability of the most likely hidden sequence

import brml.*
S = size(pSGivenSm, 1);
SHat = size(pSGivenSm, 1);

mu(:, THat) = ones(SHat, 1);
for t = THat : -1 : 2
	[m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
	XMinusMean = repmat(X(t), SHat, 1) - m;
	logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));
	tmp = repmat(logPHatXGivenS .* mu(:, t), 1, S) .* pSGivenSm;
	
	% Normalize to avoid underflow
	mu(:,t-1)= condp(max(tmp)'); 
end

% Backtrack
[m, ~, ~] = sarMean(X, Ujs, 1, A, Bjs, C);
XMinusMean = repmat(X(t), S, 1) - m;
logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));

[~, SHat(1)] = max(pS1 .* logPHatXGivenS .* mu(:, 1));
for t = 2 : THat
	if mod(t, Tskip) == 0 || Tskip == 0	
		[m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
		XMinusMean = repmat(X(t), S, 1) - m;
		logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));

		tmp = logPHatXGivenS .* pSGivenSm(:, SHat(t-1));
		[~, SHat(t)] = max(tmp .* mu(:,t));
	else
		SHat(t) = SHat(t-1);
	end
end

[m, ~, ~] = sarMean(X, Ujs, 1, A, Bjs, C);
XMinusMean = repmat(X(t), S, 1) - m;
logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));
logProbability = log(pS1(SHat(1))) + logPHatXGivenS;
for t = 2 : THat
	[m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
	XMinusMean = repmat(X(t), S, 1) - m;
	logPHatXGivenS = -0.5 * XMinusMean.^2 ./ sigma2(:) -0.5 * log(2 * pi * sigma2(:));
	logProbability = logProbability + log(pSGivenSm(SHat(t), SHat(t-1))) + logPHatXGivenS(SHat(t));
end

end