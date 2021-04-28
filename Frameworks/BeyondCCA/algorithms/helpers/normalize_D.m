function D = normalize_D(Dest)
% normalize Dest to have each column in the probability simplex
  M = size(Dest,1);
  Dest = flip_column_signs(Dest);
  Dest = max(0,Dest);
  D = Dest./repmat(sum(Dest),M,1);
end
