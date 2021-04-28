function B = compute_wtw_dcca_multiple(X1, X2, W1, W2, t1s, t2s)
%COMPUTE_WTW_DCCA_MULTIPLE Compute target matrices W1*T1(t1)*W2' and
% W1*T2(v2)*W2', speed up is due to some precomputation of terms common
% for all these matrices
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  if min(size(t1s))~=1 || min(size(t2s))~=1, error('wrong input'), end
  if max(size(t1s))~=max(size(t2s)), error('wrong input'), end
  P = max(size(t1s));
  
  W1 = sparse(W1); W2 = sparse(W2);
  precomp = precompute_for_wtw(X1,X2,W1,W2);
  
  B = [];
  for p = 1:P
    t1 = t1s{p}; t2 = t2s{p};
    [wtw1, wtw2] = compute_wtw_cca_single(X1,X2,W1,W2,t1,t2,precomp);
    B = [B wtw1 wtw2]; %#ok
  end
  
end

function [wtw121, wtw122] = compute_wtw_cca_single(X1,X2,W1,W2,t1,t2,precomp)
  
  N = size(X1,2);
  EX1 = precomp.EX1; EX2 = precomp.EX2;
  WX1 = precomp.WX1; WX2 = precomp.WX2;
  WEX1 = precomp.WEX1; WEX2 = precomp.WEX2;
  temp = precomp.temp; % 2*N*((WEX1)*(WEX2)') - (WX1)*(WX2)' 
  
  M1 = size(X1,1); M2 = size(X2,1);
  Wv1 = W1*sparse(1:M1,1:M1,t1); %diag(v1)
  Wv2 = W2*sparse(1:M2,1:M2,t2); %diag(v2)
  Xv1 = precomp.X1t*t1;
  Xv2 = precomp.X2t*t2;
  
  A1 = (WX1)*sparse(1:N,1:N,Xv1)*(WX2)';
  A2 = ((WX1)*(Xv1))*(WEX2)';
  A3 = (WEX1)*((WX2)*(Xv1))';
  A4 = (t1'*EX1)*temp;
  A5 = Wv1*precomp.X1X2W2;
  A6 = N*(Wv1*EX1)*(WEX2)';
  wtw121 = (N/((N-1)*(N-2)))*(A1-A2-A3+A4) - (1/(N-1))*(A5-A6);

  B1 = (WX1)*sparse(1:N,1:N,Xv2)*(WX2)';
  B2 = ((WX1)*(Xv2))*(WEX2)';
  B3 = (WEX1)*((WX2)*(Xv2))';
  B4 = (t2'*EX2)*temp;
  B5 = precomp.W1X1X2*Wv2';
  B6 = N*(WEX1)*(Wv2*EX2)';
  wtw122 = (N/((N-1)*(N-2)))*(B1-B2-B3+B4) - (1/(N-1))*(B5 - B6);

end

function precomp = precompute_for_wtw(X1,X2,W1,W2)
  EX1 = mean(X1,2);
  EX2 = mean(X2,2);
  WX1 = W1*X1;
  WX2 = W2*X2;
  W1X1X2 = WX1*X2';
  X1X2W2 = X1*WX2';
  WEX1 = W1*EX1;
  WEX2 = W2*EX2;
  N = size(X1,2);
  temp = 2*N*((WEX1)*(WEX2)') - (WX1)*(WX2)';
  X1t = X1'; X2t = X2';
  precomp = struct('WX1',WX1,'WX2',WX2,'EX1',EX1,'EX2',EX2, ...
                   'WEX1',WEX1,'WEX2',WEX2,'temp',temp, ...
                   'X1t',X1t,'X2t',X2t,...
                   'W1X1X2',W1X1X2,'X1X2W2',X1X2W2);
end
