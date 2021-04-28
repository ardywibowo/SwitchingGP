function [maxstate logprob]=HMMviterbiGeneral(phghm,ph1,phiht)
%HMMVITERBI Viterbi most likely joint hidden state of a HMM for an
%inhomogenous emission distribution
%
% [maxstate logprob]=HMMviterbiGeneral(phghm,ph1,phiht)
%
% Inputs:
% v : visible (obervation) sequence being a vector v=[2 1 3 3 1 ...]
% phghm : homogeneous transition distribution phghm(i,j)=p(h(t)=i|h(t-1)=j)
% ph1 : initial distribution
% phiht : p(v,t)\propto phi(h,t) (for a fixed observed sequence, the
% emission distribution depends only on h and potentially t as well for an
% inhomogenous emission distribution)
%
% Outputs:
% maxstate : most likely joint hidden (latent) state sequence
% logprob : associated log probability of the most likely hidden sequence
% See also demoHMMinference.m
import brml.*
T=size(phiht,2); H=size(phghm,1);
mu(:,T)=ones(H,1);
for t=T:-1:2
	tmp = repmat(phiht(:,t).*mu(:,t),1,H).*phghm;
	mu(:,t-1)= condp(max(tmp)'); % normalise to avoid underflow
end
% backtrack
[val hs(1)]=max(ph1.*phiht(:,1).*mu(:,1));
for t=2:T
	tmp = phiht(:,t).*phghm(:,hs(t-1));
	[val hs(t)]=max(tmp.*mu(:,t));
end
maxstate=hs;
logprob=log(ph1(hs(1)))+log(phiht(hs(1,t)));
for t=2:T
	logprob=logprob+log(phghm(hs(t),hs(t-1)))+log(phiht(hs(t),t));
end