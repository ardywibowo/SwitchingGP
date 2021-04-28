function S12 = compute_S12_ncca(X1,X2)
  Ex1 = mean(X1,2); Ex2 = mean(X2,2);
  N = size(X1,2);
  S12 = (1/(N-1)) * (X1*X2' - N* Ex1*Ex2');
end