function gW=gradsoftmaxloss(allW,X,C,sW)
W=reshape(allW,sW);
gW=brml.gradlogsoftmax(X,C,W);
gW=-gW(:);