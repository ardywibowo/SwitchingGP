function [D1, D2] = dcca(X1, X2, K)
%DCCA Discrete canonical correlation analysis; the estimation is based 
% on the method of moments of the cumulant-based tensor T with  
% non-orthogonal joint diagonalization by similarity
%
% [D1, D2] = dcca(X1, X2, K)
%
% Input: 
%   X1 : a first sparse M1-by-N data matrix with observations in columns
%   X2 : a second sparse M2-by-N data matrix with observations in columns
%   K  : the number of latent variables (topics)
%
% Output:
%   D1 : the estimation of the parameter D1 (M1-by-K matrix)
%   D2 : the estimation of the second parameter D2 (M2-by-K matrix)
%
% The model and algorithm are proposed by:
%   A. Podosinnikova, F. Bach, S. Lacoste-Julien. Beyond CCA: 
%   Moment matching for multi-view models.
%   In Proc. ICML, 2016.
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  if size(X2,2)~=size(X1,2), error('wrong input'), end
  if ~(issparse(X1)||issparse(X2)), X1 = sparse(X1); X2 = sparse(X2); end
  
  % Compute whitening matrices
  S12 = compute_S12_dcca(X1,X2);
  [W1, W2] = compute_whitening_matrices(S12, K);
  
  % Initialize the projection vectors
  v1s = cell(K,1); v2s = cell(K,1);
  E = eye(K);
  for k=1:K, v1s{k} = W1'*E(:,k); v2s{k} = W2'*E(:,k); end
  
  % Construct the target matrices to be jointly diagonalized
  B = compute_wtw_dcca_multiple(X1, X2, W1, W2, v1s, v2s);
  B = [ eye(K) B ];
  
  % Non-orthogonal joint diagonalization
  [Q, Qinv] = nojd_in_cpp( B );
  D1 = normalize_D( pinv(W1)*Q );
  D2 = normalize_D( pinv(W2)*Qinv' );
  
end
