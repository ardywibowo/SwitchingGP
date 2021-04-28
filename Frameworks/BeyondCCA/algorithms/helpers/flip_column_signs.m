function [Dest, signs] = flip_column_signs(Dest)
%FLIP_COLUMN_SIGNS
%
%
% Copyright 2016, Anastasia Podosinnikova

  K = size(Dest,2);
  signs = ones(1,K);
  % if in the column more negative values than positive => switch the sign
  % can be implemented in many different ways
  for k = 1:K
    val1 = sum(min(0,Dest(:,k)).^2);
    val2 = sum(max(0,Dest(:,k)).^2);
    signs(k) = sign(val2 - val1);
  end
  Dest = Dest*diag(signs);
  
end
