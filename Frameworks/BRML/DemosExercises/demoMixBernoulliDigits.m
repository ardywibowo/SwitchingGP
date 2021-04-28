function demoMixBernoulliDigits(H)
%DEMOMIXBERNOULLIDIGITS demo of training a mixture of Bernoulli distributions
import brml.*
v=[];
for d=0:9
    load(['digit',num2str(d),'.mat'])
    xx=x';
    v=[v xx(:,1:50)];
end
v=v>0;

figure; disp('plotting subset of training data')
r=randperm(size(v,2));
v=v(:,r);
for h=1:40*5
    subplot(5,40,h,'align');
    imagesc(reshape(v(:,h),28,28)'); axis off; colormap bone
end


% EM training:
%H=10;
opts.plotprogress=1; opts.maxit=10;
opts.meaninit=0;
% Perform multiple runs (due to local maxima issues:)
loglik=-inf; nruns=3;

for runs=1:nruns
    figure(2); subplot(nruns,1,runs);
    [phr pvghr thisloglik phgvr]=MIXprodBern(v,H,opts);
    thisloglik
    if thisloglik>loglik
        ph=phr; pvgh=pvghr; phgv=phgvr; loglik=thisloglik;
    end
    figure(3); 
    for h=1:H
        subplot(nruns,H,h+(runs-1)*H,'align');
        imagesc(reshape(pvgh(:,h),28,28)'); axis off;% colormap bone
    end
end
disp('plotted Mixture components')