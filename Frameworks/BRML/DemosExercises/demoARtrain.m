function demoARtrain
close all
import brml.*
T=100; t=1:T;
v = 0.02*(t-2).^2-0.2.*t; v = v + 7*randn(1,T); % noisy data
plot(t,v,'.');

L=3; % AR order
[a res]=ARtrain(v,L); % fit the AR model
% prediction:
mv(1:T)=v(1:T);TT=T+40;
for t=T+1:TT
    vhat = mv(t-L:t-1);
    mv(t) = vhat*a;
end
hold on; plot(T+1:TT,mv(T+1:TT),'-m')
s=sqrt(mean(res.^2));%plot(T+1:TT,mv(T+1:TT)-s,':m');plot(T+1:TT,mv(T+1:TT)+s,':m')

% prediction:
Nsamp=1000;
vv(:,1:T)=repmat(v(1:T),Nsamp,1);TT=T+40;
for samp=1:Nsamp
    for t=T+1:TT
        vhat = vv(samp,t-L:t-1);
        vv(samp,t) = vhat*a+s*randn;
    end
end
sd=std(vv);
plot(T+1:TT,mv(T+1:TT)-sd(T+1:TT),':m');plot(T+1:TT,mv(T+1:TT)+sd(T+1:TT),':m')
