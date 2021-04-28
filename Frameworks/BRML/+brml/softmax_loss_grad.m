function [E,g]=softmax_loss_grad(par,X,C,sW)
% SOFTMAX_LOSS_GRAD : wrapper for minFunc
E=brml.softmaxloss(par,X,C,sW);
g=brml.gradsoftmaxloss(par,X,C,sW);