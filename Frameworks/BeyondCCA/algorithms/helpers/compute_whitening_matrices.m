function [W1, W2] = compute_whitening_matrices(S12, K)
%COMPUTE_WHITENING_MATRICES Compute whitening matrices W1 and W2 of S12
%
%
% Copyright 2016, Anastasia Podosinnikova

  [U, D, V] = svds(S12,K);
  tmp = pinv( sqrt(D) );
  W1 = tmp * U';
  W2 = tmp * V';
  
end