function gW=gradlogsoftmax(X,C,W)
% GRADSOFTMAX: gradient of the log softmax probability
gW=zeros(size(W));
S=brml.mysoftmax(X,W);
for n=1:length(C)
    gW=gW-S(:,n)*X(:,n)';
    gW(C(n),:)=gW(C(n),:)+X(:,n)';
end