function [D1est, D2est] = get_dcca_matrices_oracle( Dest, M1, M2 )
%GET_DCCA_MATRICES_ORACLE Separates single matrix D into D1 and D2 
% row-wise; oracle means that we know that Dest = [D1est; D2est], i.e. that
% the first M1 rows correspond to D1 and the consequent M2 rows correspond 
% to D2.
%
% Copyright 2016, Anastasia Podosinnikova
  
  D1est = Dest( 1:M1, : );
  D1est = renormalize( D1est );

  D2est = Dest( M1+1:M1+M2, : );
  D2est = renormalize( D2est );

end

function Dest = renormalize(Dest)
  M = size(Dest,1);
  Dest( :, logical( sum(Dest) < M*1e-14 ) ) = 1;
  Dest = Dest ./ repmat( sum(Dest), M, 1 );
end
