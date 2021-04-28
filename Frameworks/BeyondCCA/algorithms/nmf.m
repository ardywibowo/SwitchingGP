function [D, H] = nmf(X, D0, H0)
%NMF Nonnegative matrix factorization
%
% [D, H] = nmf(X, D0, H0)
%
% Input:
%   X  : the M-by-N nonnegative matrix to be factorized as D*H 
%        with D and H nonnegative
%   D0 : the M-by-K (nonnegative) initialization for D
%   H0 : the K-by-N (nonnegative) initialization for H
%
% Output:
%   D    : a M-by-K nonnegative matrix
%   H    : a K-by-N nonnegative matrix
% 
% This is an implementation of the algorithm proposed by:
%   D.D. Lee, H.S. Seung. Algorithms for non-negativ matrix factorization. 
%   In Adv. NIPS, 2000.
%
%
% Copyright 2016, Anastasia Podosinnikova

  if sum(sum( D0 < 0 ))>1  ||  sum(sum( H0 < 0 ))>1 || sum(sum( X < 0 ))>1
    error('Wrong input')
  end
  
  maxiter = 1e3;
  epsilon = 1e-8;
  D = D0; H = H0;
  [M, N] = size(X);
  
  for iter = 1:maxiter
    % Multiplicative update rule
    D = D.*(( X./(D*H+eps))*H')./repmat(sum(H,2)',M,1);
    H = H.*(D'*(X./(D*H+eps)))./repmat(sum(D,1)',1,N);
    
    % Renormalize each term to avoid extremely large/small values
    H = H.*repmat(sum(D,1)',1,N);
    D = D.*repmat(1./sum(D,1),M,1);
    
    % Stopping criterion
    if norm(X - D*H)/norm(X) <= epsilon,
      disp(['NMF terminated after ', num2str(iter), ' iteration'])
      break;
    end
  end
  
end
