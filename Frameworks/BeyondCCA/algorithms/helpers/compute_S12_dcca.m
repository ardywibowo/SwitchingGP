function S12 = compute_S12_dcca(X1,X2)
%COMPUTE_S12_DCCA Compute the cross-covariance matrix S12 of the DCCA model
%
%
% Copyright 2016, Anastasia Podosinnikova

  Ex1 = mean(X1,2); Ex2 = mean(X2,2);
  N = size(X1,2);
  S12 = (1/(N-1)) * (X1*X2' - N* Ex1*Ex2');
end