function Wopt=trainsoftmax(W,X,C,options)
% TRAINSOFTMAX: softmax classification
% For the n-th datapoint:
% p(class(n)=c|X(:,n))\propto exp(W(c,:)*x(:,n)
% W : inital setting of the weight matrix
% X : dataset (each column contains a datapoint)
% C : vector of class labels (each label is an integer 1,2,...)
% options: minFunc options
%
% See demoSoftMaxRawDigits.m
sW=size(W);
objfcn = @(par) brml.softmax_loss_grad(par,X,C,sW);
[allpar, funval, exitflag, output] = minFunc(objfcn, W(:), options);
Wopt=reshape(allpar,size(W));