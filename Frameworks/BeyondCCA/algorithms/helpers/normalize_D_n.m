function D = normalize_D_n(Dest)
% normalize each column of Dest to have unit l1-norm
  K = size( Dest, 2 );
  D = zeros( size( Dest ) );
  for k = 1:K
    D(:,k) = Dest(:,k) / norm( Dest(:,k), 1 );
  end
end