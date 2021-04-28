function precomputed = precompute_for_wtw(SX, W, momtype, M1, M2)
%PRECOMPUTE_FOR_WTW Performs some precompuations for further computation of
%   W*T(v)*W' where T is either the DICA T-cumulant or LDA T-moment
%
% precomputed = precompute_for_wtw(SX, W, momtype, M1, M2)
%
% INPUT:
%   SX      : sparse M-by-N matrix of word counts X with docs in columns
%   W       : K-by-M whitening matrix of the DICA S-cumulant
%   momtype : either 'lda' for LDA moments or 'dica' for DICA cumulants
%   M1      : as it is defined for either LDA moments or DICA cumulants
%   M2      : the second LDA moment, only for LDA moments (momtype='lda')
%
% OUTPUT:
%   precomputed : struct with precomputed values
%
% COMMENTS: When computing P K-by-K matrices W*T(v_1)*W', W*T(v_2)*W', 
%   ..., W*T(v_P)*W', it is sufficient to perform some computations only 
%   ones instead of P times. The result of these computations is contained 
%   in the variable called precomputed. M is the number of words in the 
%   dictionary, N is the number of documents in the corpus, and K is the 
%   number of topics.
%
%
% Copyright 2016, Anastasia Podosinnikova

  if ~( strcmp(momtype,'dica') || strcmp(momtype,'lda') )
    error('Wrong momtype')
  end
  
  if strcmp(momtype,'dica')
    
    N    = size(SX,2);
    W    = sparse(W); 
    WX   = W*SX;
    WM1  = W*M1;
    temp = (2*N*WM1)*WM1' - WX*WX';
    Xt   = SX';
    XWX  = SX*WX';
    
    precomputed = struct('WX',WX, ...
                         'WM1',WM1, ...
                         'temp',temp, ...
                         'Xt',Xt, ...
                         'XWX',XWX);
  end
  
  if strcmp(momtype,'lda')
     
    N      = size(SX,2);
    Ls     = sum(SX)';
    delta3 = 1./( Ls.*(Ls-1).*(Ls-2) );
    W      = sparse(W); 
    WX     = W*SX;
    WXd3t  = (WX*sparse(1:N,1:N,delta3))';
    Xd3    = SX*delta3;
    WM2W   = W*M2*W';
    WM1    = W*M1;
    WM1WM1 = WM1*WM1';
    
    precomputed = struct('delta3',delta3, ...
                         'WX',WX, ...
                         'WXd3t',WXd3t, ...
                         'Xd3',Xd3, ...
                         'WM2W',WM2W, ...
                         'WM1',WM1, ...
                         'WM1WM1',WM1WM1);
  end
  
end