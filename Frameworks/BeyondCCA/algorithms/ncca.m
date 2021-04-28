function [D1, D2] = ncca( X1, X2, K, scale )
%NCCA Non-Gaussian canonical correlation analysis; the estimation 
% is based on non-orthogonal joint diagonalization by similarity of 
% generalized covariance matrices
%
% [D1, D2] = ncca( X1, X2, K, scale )
%
% Input: 
%   X1    : a first sparse M1-by-N data matrix with observations in columns
%   X2    : a second sparse M2-by-N data matrix with observations in columns
%   K     : the number of latent variables (sources)
%   scale : scale factor for the processing points of generalized covariance matrices
%           NOTE: the algroithm is sensitive to the choice of this scaling
%                 parameter (update of the code is needed, contact me)
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
  
  type = 'ncca';
  if size(X1,2)~=size(X2,2), error('wrong input'), end
  N = size(X1,2);
  
  % Compute whitening matrices
  S12 = compute_S12_ncca(X1,X2); % <=> S12(t) for t=0
  [W1, W2] = compute_whitening_matrices(S12, K);
  
  % Construct target matrices based on generalized covariance matrices
  E = eye(K);
  B = zeros(K, K*K);
  for k=1:K
    t1 = scale * (W1'*E(:,k)); 
    t2 = scale * (W2'*E(:,k));
    wsw12 = compute_gen_wsw12_cca_all( X1, X2, W1, W2, t1, t2, N, type );
    B(:,K*(k-1)+1:K*k) = wsw12;
  end
  B = [W1*S12*W2' B];
  
  % Non-orthogonal joint diagonalization
  [Q, Qinv] = nojd_in_cpp( B );
  D1 = normalize_D_n(pinv(W1)*Q);
  D2 = normalize_D_n(pinv(W2)*Qinv');

end
