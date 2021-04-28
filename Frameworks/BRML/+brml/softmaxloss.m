function L=softmaxloss(allW,X,C,sW)
W=reshape(allW,sW);
L=-brml.logsoftmax(X,C,W);