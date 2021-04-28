function demoSoftmaxRawDigits(varargin)
% DEMO of softmax classification for MNIST digita
% (Note that this is not a great method for the digits problem)
% demoSoftmaxRawDigits(<minFuncPath>)
import brml.*
if nargin==0
    minFuncPath='../minFunc';
else
    minFuncPath=varargin{1};
end

addpath(genpath(minFuncPath))

options.Method = 'lbfgs';
options.progTol=10e-25;
options.optTol=10e-25;
options.DerivativeCheck='off';
options.MaxIter = 50;
options.MaxFunEvals=50;

if ~exist('mnist-original.mat')
    f='http://mldata.org/repository/data/download/matlab/mnist-original/';
    fprintf(1,['\nAttempting to download the MNIST data from:\n',f])
    urlwrite(f,'mnist-original.mat')
end
load mnist-original
data=double(data); data=data./max(data(:))-0.5;

N=60000; xtrain=data(:,1:N);
xtest=data(:,N+1:end);

W=0.001*randn(10,size(xtrain,1));
Wopt=trainsoftmax(W,xtrain,1+label(1:N),options);

probtrain=mysoftmax(xtrain,Wopt);
[val ind]=max(probtrain);
ctrain=ind-1;
fprintf(1,'\nTrain set accuracy = %f',mean(ctrain==label(1:N)))

probtest=mysoftmax(xtest,Wopt);
[val ind]=max(probtest);
ctest=ind-1;
fprintf(1,'\nTest set accuracy  = %f\n',mean(ctest==label(N+1:end)))