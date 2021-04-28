function [alpha,loglik]=HMMforwardTimeDependent(v,phghm,ph1,pvgh)
%HMMFORWARD HMM Forward Pass
% [alpha,loglik]=HMMforward(v,phghm,ph1,pvgh)
%
% Inputs:
% v : visible (observation) sequence being a vector v=[2 1 3 3 1 ...]
% phghm : transition distribution phghm(i,j,t)=p(h(t)=i|h(t-1)=j,t)
% ph1 : initial distribution
% pvgh : emission distribution pvgh(i,j,t)=p(v(t)=i|h(t)=j,t)
%
% Outputs:
% alpha : alpha messages: p(h(t)|v(1:t))
% loglik : sequence log likelihood log p(v(1:T))
% See also HMMbackward.m, HMMviterbi.m,  HMMsmooth.m, demoHMMinference.m
import brml.*
T=length(v); H=length(ph1);

% alpha recursion (with normalisation to avoid numerical underflow)
alpha=zeros(H,T);
z=zeros(1,T); % local normalisation factors
alpha(:,1) = pvgh(v(1),:,1)'.*ph1;
z(1)=sum(alpha(:,1));
alpha(:,1)=alpha(:,1)/z(1);
for t=2:T
    alpha(:,t)=pvgh(v(t),:,t)'.*(phghm(:,:,t)*alpha(:,t-1));
    z(t)=sum(alpha(:,t));
    alpha(:,t)=alpha(:,t)/z(t);
end
loglik = sum(log(z(:))); % log likelihood