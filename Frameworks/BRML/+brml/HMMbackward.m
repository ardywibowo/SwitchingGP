function beta=HMMbackward(v,phghm,pvgh)
%HMMBACKWARD HMM Backward Pass
% beta=HMMbackward(v,phghm,pvgh)
%
% Inputs:
% v : visible (observation) sequence being a vector v=[2 1 3 3 1 ...]
% phghm : homogeneous transition distribution phghm(i,j)=p(h(t)=i|h(t-1)=j)
% ph1 : initial distribution
% pvgh : homogeneous emission disrtribution pvgh(i,j)=p(v(t)=i|h(t)=j)
% 
% Outputs:
% beta: beta messages: p(v(t+1:T)|h(t))/(sum_{h(t)} p(v(t+1:T)|h(t)))
% See also HMMbackward.m, HMMviterbi.m, demoHMMinference.m
import brml.*
T=length(v); H=size(phghm,1);
if 1==0
% logbeta recursion
logbeta(:,T)=zeros(H,1);
for t=T:-1:2
	logbeta(:,t-1)=logsumexp(repmat(logbeta(:,t),1,H),repmat(pvgh(v(t),:)',1,H).*phghm);
end
end
if 1==1 % beta recursion (not recommended due to numerical underflow)
	beta(:,T)=ones(H,1);
	for t=T:-1:2
		beta(:,t-1)=phghm'*(beta(:,t).*pvgh(v(t),:)');
        beta(:,t-1)=beta(:,t-1)/sum(beta(:,t-1));
	end
end