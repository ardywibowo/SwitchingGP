function [alpha,loglik]=HMMforward(v,phghm,ph1,pvgh)
%HMMFORWARD HMM Forward Pass
% [alpha,loglik]=HMMforward(v,phghm,ph1,pvgh)
%
% Inputs:
% v : visible (observation) sequence being a vector v=[2 1 3 3 1 ...]
% phghm : homogeneous transition distribution phghm(i,j)=p(h(t)=i|h(t-1)=j)
% ph1 : initial distribution
% pvgh : homogeneous emission distribution pvgh(i,j)=p(v(t)=i|h(t)=j)
% 
% Outputs:
% alpha : alpha messages: p(h(t)|v(1:t))
% loglik : sequence log likelihood log p(v(1:T))
% See also HMMbackward.m, HMMviterbi.m,  HMMsmooth.m, demoHMMinference.m
import brml.*
T=length(v); H=length(ph1);
if 1==0
% logalpha recursion
logalpha(:,1) = log(pvgh(v(1),:)'.*ph1);
for t=2:T
	logalpha(:,t)=logsumexp(repmat(logalpha(:,t-1),1,H),repmat(pvgh(v(t),:),H,1).*phghm');
end
loglik = logsumexp(logalpha(:,T),ones(H,1)); % log likelihood
end

% alpha recursion (with normalisation to avoid numerical underflow)
    alpha=zeros(H,T); 
    z=zeros(1,T); % local normalisation factors
	alpha(:,1) = pvgh(v(1),:)'.*ph1;
    z(1)=sum(alpha(:,1));
    alpha(:,1)=alpha(:,1)/z(1);
	for t=2:T
		alpha(:,t)=pvgh(v(t),:)'.*(phghm*alpha(:,t-1));
        z(t)=sum(alpha(:,t));
        alpha(:,t)=alpha(:,t)/z(t);
	end
	loglik = sum(log(z(:))); % log likelihood
end