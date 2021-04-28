function WTWs = compute_multiple_wtw_dica(SX, W, K, vecs, M1)
%COMPUTE_MULTIPLE_WTW_DICA Computes P K-by-K matrices W*T(v_1)*W', 
%   W*T(v_2)*W', ..., W*T(v_P)*W' for the DICA T-cumulant and P vectors 
%   v_1, v_2, ..., v_P
%
% WTWs = compute_multiple_wtw_dica(SX, W, K, vecs, M1)
%
% Input: 
%   SX   : sparse M-by-N matrix of word counts X with docs in columns
%   W    : K-by-M whitening matrix of the DICA S-cumulant
%   K    : number of topics
%   vecs : cell array {v_1, v_2, ..., v_P}
%   M1   : expectation of document's counts
%
% Output:
%   WTWs : K-by-(P*K) matrix [W*T(v_1)*W' W*T(v_2)*W' ... W*T(v_P)*W']
%
% Comment: M is the number of words in the dictionary, N is the number
%   of documents in the corpus, and P is the number of projections, i.e.
%   number of vectors, v_1, v_2, ..., v_P. Note that in most cases it is
%   beneficial to set each vector v_p to W'*u_p, where u_p is a K-vector. 
%   It is often sufficient to set P = K and vectors u_p to the canonical 
%   basis of R^K, i.e. the colums of the K-identity matrix.
%
%
% Copyright 2016, Anastasia Podosinnikova

  P = length(vecs);
  W = sparse(W);
  precomputed = precompute_for_wtw(SX, W, 'dica', M1);
  
  WTWs = zeros(K,K*P);
  for p = 1:P
    v = vecs{p};
    WTWs(:,K*(p-1)+1:K*p) = compute_wtw_single_dica(SX,W,v,M1,precomputed);
  end
  
end