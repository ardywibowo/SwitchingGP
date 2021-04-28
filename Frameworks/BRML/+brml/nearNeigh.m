function [y, tie] = nearNeigh(xtrain, xtest, trainlabels, K)
%NEARNEIGH Nearest Neighbour classification
% y=nearNeigh(xtrain, xtest, trainlabels,K)
% calculate the nearest  neighbour classification (use squared distance to measure dissimilarity)
% If there is a tie, the single nearest neighbour class is returned
% xtrain : matrix with each column a training vector
% xtest : matrix with each column a test vector
% trainlabels : vector of length size(xtrain,2) of training labels
ntrain = size(xtrain,2); % number of training points
ntest = size(xtest,2); % number of test points
[vals, ind] = sort(brml.sqdist(xtrain,xtest));
tie=zeros(1,size(xtest,2));
if K==1
    y =trainlabels(ind(1,:));
else
    if size(xtest,2)>1
        [y,tie] = brml.majority(trainlabels(ind(1:K,:)),1);
    else
        [y,tie] = brml.majority(trainlabels(ind(1:K,:))',1);
    end
end
if sum(tie)>0
    t=trainlabels(ind(1,:));
    y(tie)=t(tie);
end