function [ph pvgh loglik phgv]=MIXprodCategorical(v,H,opts)
%MIXPRODCATEGORICAL EM training of a Mixture of a product of Categorical distributions
%[ph pvgh loglik phgv]=MIXprodCategorical(v,H,opts)
%
% Inputs:
% v : data matrix : each column is a datapoint
% H : number of mixture components
% opts.maxit
% opts.plotprogress
% opts.meanint (set to true to initalise to the mean of each catgory)
%
% Outputs:
% ph : p(h)
% pvgh : p(v|h)
% loglik : log likelihood of the set of data
% phgv : p(h|v) posterior assignment of datapoints to mixture components
% See also demoMixprodCategorical.m
import brml.*
[D N]=size(v);
ph = condp(rand(H,1)); pvgh=condp(rand(D,H,opts.C),3); % random initialisations
if isfield(opts,'meaninit')
    if opts.meaninit
        ph=brml.condp(ones(H,1));
        for d=1:D
            for c=1:opts.C
                pvgh(d,:,c)=sum(v(d,:)==c);
            end
        end
        pvgh=brml.condp(pvgh+0.001*rand(size(pvgh)),3);
    end
end
for emloop = 1:opts.maxit
    %E-step:
    htot=zeros(H,1);
    vhtot=zeros(D,H,opts.C);
    loglik=0;
    for n = 1:N
        logpold = log(ph);
        for c=1:opts.C
            st{c}=find(v(:,n)==c);
            if ~isempty(st{c})
                logpold = logpold+sum(log(pvgh(st{c},:,c)),1)';
            end
        end
        poldhgvn=condexp(logpold); phgv(:,n)=poldhgvn;
        % get the stats for the M-step:
        htot=htot+poldhgvn;
        for c=1:opts.C
            vhtot(st{c},:,c)=vhtot(st{c},:,c)+repmat(poldhgvn',length(st{c}),1);
        end
        loglik=loglik+logsumexp(logpold,ones(H,1));
    end
    lik(emloop)=loglik;
    if opts.plotprogress; plot(lik,'-o'); title('log likelihood'); xlabel('iteration'); drawnow; end
    ph=condp(htot); 	pvgh=condp(vhtot,3); % M-step
end