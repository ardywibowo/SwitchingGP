function y=logsoftmax(X,C,W)
% LOGSOFTMAX: log of the softmax probability
y=0;
for n=1:length(C)
    y=y+W(C(n),:)*X(:,n);
end
y=y-sum(log(sum(exp(W*X))));