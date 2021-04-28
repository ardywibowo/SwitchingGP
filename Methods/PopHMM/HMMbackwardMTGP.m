function logBeta = HMMbackwardMTGP(X, Ujs, pSGivenSm, A, Bjs, C, sigma2, Tskip)
%HMMBACKWARDSAR Backward Pass (beta method) for the Switching Autoregressive HMM
% logBeta = HMMbackwardInputSAR(X, Ujs, pSGivenSm, A, Bjs, C, sigma2, Tskip)
%
% Inputs:
% X         : State observations
% Ujs       : Input observations
% pSGivenSm : State (switch) transition matrix p(s(t)|s(t-1))
% A         : Matrix of AR coefficients for state x. Column a(:,i) are the AR coeffs for switch state i. 
% Bjs       : Cell vector of input AR coefficients. Each cell is one input variable. 
%             Each column Bj(:,i) are the AR coeffs for switch input i.
% C         : Column vector of constant coefficients
% The AR coefficients are in reverse order X(t-L) : X(t-1)
% Each column represents a different switching state
% sigma2    : The innovation noise
% Tskip     : The number of timesteps to skip before a switch update is allowed
%
% Outputs:
% logBeta   : Log backward messages log p(x(t+1:T)|s(t),x(t-L+1:t))

T = length(X);
S = size(A, 2);

% logBeta recursion
logBeta(:, T) = zeros(S, 1);
for t = T : -1 : 2
    [m, ~, ~] = sarMean(X, Ujs, t, A, Bjs, C);
    XMinusMean = repmat(X(t), S, 1) - m;
	
    pHatXGivenS = exp(-0.5 * XMinusMean.^2 ./ sigma2(:)) ./ sqrt(2 * pi * sigma2(:)) + eps;
    if mod(t, Tskip) == 0 || Tskip == 0
        pSgSmt = pSGivenSm;
    else
        pSgSmt = eye(S);
    end
    logBeta(:, t-1) = brml.logsumexp(repmat(logBeta(:, t), 1, S), repmat(pHatXGivenS, 1, S) .* pSgSmt);
end

end

