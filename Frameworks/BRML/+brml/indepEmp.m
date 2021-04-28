function [ind,logBayesFactor]=indepEmp(dataX,dataY,X,Y,thresh,varargin)
%INDEPEMP Compute the empirical log Bayes Factor and MI for independence/dependence
% [ind,logBayesFactor]=condindepEmp(dataX,dataY,X,Y,thresh,<opts>)
% here dataX,dataY are data matrices where each row contains the sample states
% with the number of their states in X,Y
% opts.Ux, opts.Uy, are the (scalar) hyperparameters
import brml.*
if nargin==5
    alpha=0.1*size(dataX,2);
    Ux=alpha/prod(X); Uy=alpha/prod(Y); Uxy=alpha/(prod(X)*prod(Y));
else
    opts=varargin{1};
    Ux=opts.Ux; Uy=opts.Uy;
end

% model of p(x)
cx=count(dataX,X);
logZux = logZdirichlet(Ux*ones(prod(X),1));
logpx = logZdirichlet(cx+Ux)-logZux;

% model of p(y)
cy=count(dataY,Y);
logZuy = logZdirichlet(Uy*ones(prod(Y),1));
logpy = logZdirichlet(cy+Uy)-logZuy;

logpindep = logpx+logpy;

% model of p(xy)
cxy=count(vertcat(dataX,dataY),[X Y]);
logZuxy = logZdirichlet(Uxy*ones(prod([X Y]),1));
logpdep = logZdirichlet(cxy+Uxy)-logZuxy;

logBayesFactor=logpindep-logpdep;

ind=logBayesFactor>thresh; % use the Bayes Factor to decide