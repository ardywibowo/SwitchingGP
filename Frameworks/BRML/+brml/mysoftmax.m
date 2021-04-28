function y=mysoftmax(x,W)
% MYSOFTMAX(X,W)
% softmax probability p(x=c|W)\propto exp(W(c,:)*x)
% X : each column contains a D dimensional datapoint
% W : CxD matrix of softmax weights
y=exp(W*x)./repmat(sum(exp(W*x)),size(W,1),1);