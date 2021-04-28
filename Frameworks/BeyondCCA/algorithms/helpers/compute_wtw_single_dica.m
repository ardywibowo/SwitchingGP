function WTW = compute_wtw_single_dica(SX, W, v, M1, precomputed)
% COMPUTE_WTW_SINGLE_DICA Computes W*T(v)*W' for the DICA T-cumulant
%
% WTW = compute_wtw_single_dica(SX, W, v, M1, precomputed)
%
% INPUT:
%   SX : sparse M-by-N matrix of word counts X with docs in columns
%   W  : K-by-M whitening matrix of the DICA S-cumulant
%   v  : M-vector, in practice, W'*u where u is a K-vector
%   M1 : expectation of document's counts
%   precomputed : some precomputed values (see precompute_for_wtw.m)
%
% OUTPUT:
%   WTW : K-by-K matrix W*T(v)*W'
%
% COMMENTS: This function computes W*T(v)*W', where T is the DICA
%   T-cumulant, W is a whitening matrix of the DICA S-cumulant, and v is a
%   vector. When computing several K-by-K matrices W*T(v_1)*W', 
%   W*T(v_2)*W', ..., W*T(v_P)*W', it is sufficient to perform some 
%   computations only ones instead of P times. The result of these 
%   computations is contained in the variable called precomputed. M is the 
%   number of words in the dictionary, N is the number of documents in the 
%   corpus, and K is the number of topics.
%
%
% Copyright 2016, Anastasia Podosinnikova

  [M, N] = size(SX);
  
  WX   = precomputed.WX;
  WM1  = precomputed.WM1;
  temp = precomputed.temp;
  Xt   = precomputed.Xt;
  XWX  = precomputed.XWX;
  
  Xv  = Xt*v;
  M1v = M1'*v;
  Wv  = W*sparse(1:M,1:M,v);

  temp0 = Wv*XWX;
  temp1 = WM1 * (WX*Xv)';
  temp2 = N*WM1*(Wv*M1)';

  WTW = (N/((N-1)*(N-2))) * ( ...
                    WX*sparse(1:N,1:N,Xv)*WX' ...
                  + M1v * temp ...
                  - (temp1 + temp1')) ...
        + 2*W*sparse(1:M,1:M,v.*M1)*W' ...
        + (1/(N-1)) * ( ...
                  - temp0 - temp0' ...
                  + temp2 + temp2' ...
                  + W*sparse(1:M,1:M,N*M1v*M1 - SX*Xv)*W') ...
                ;

  
end