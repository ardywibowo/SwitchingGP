function [Dest, as, W, vecs] = dica(X, K)
%DICA Discrete independent component analysis; the evaluation is based 
% on the method of moments and the orthogonal joint diagonalization 
% algorithm with pre-whitening
%
% [Dest, as, W, vecs] = dica(SX, K)
%
% Input:
%   X : a sparse M-by-N data matrix with observations in columns
%   K : the number of latent variables (topics)
%
% Output:
%   Dest : an estimate of the parameter D
%   W    : the K-by-M whitening matrix
%   vecs : the projection vectors for the cumulant-based tensor T
%
% The model and algorithm are proposed by:
%   A. Podosinnikova, F. Bach, S. Lacoste-Julien. Rethinking LDA: 
%   Moment matching for discrete ICA.
%   In Adv. NIPS, 2015.
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  if ~issparse(X), X = sparse(X); end
  
  % Compute a whitening matrix
  [W, M1] = compute_S_and_W_dica(X, K); 
  
  % Intialize the projection vectors
  E = eye(K);
  vecs = cell(K,1);
  for k = 1:K, vecs{k} = W'*E(:,k); end
  
  % Compute target matrices to be jointly diagonalized,
  % i.e. (W*T(v)*W', for different v)
  WTWs = compute_multiple_wtw_dica(X, W, K, vecs, M1);
  
  % Joint diagonalization
  [V, as] = jd_in_cpp([eye(K) WTWs]); % Note that eye(K) = W*S*W'
  
  
  % Computing the parameter D
  Dest = pinv(V'*W);
  % Problems:
  %   - the pseudo-inverse introduces negative values 
  %   - each column of Dest is estimated up to the multiplication by a scalar
  %         => checking wheter columns have correct signs
  Dest = flip_column_signs(Dest);
  % Truncating all negative values
  Dest = max(0,Dest);
  
  % WARNING: the normalization is done outside of this function. This is
  % necessary to correctly separate D into D1 and D2. Uncomment the last
  % line to allow the normalization inside of this function.
  % Dest = Dest ./ repmat(sum(Dest), M, 1);
  
end
