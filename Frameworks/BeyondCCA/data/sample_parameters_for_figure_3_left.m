function params = sample_parameters_for_figure_3_left
%SAMPLE_PARAMETERS_FOR_FIGURE_3_LEFT Samples parameters for the experiment 
% from Figure 3 (left) of the supplementary material (Appendix) of the 
% paper.
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  % see the paper for the notation
  M1 = 20; 
  M2 = 20; 
  K1 = 1; 
  K2 = 1; 
  K = 1;
  c0 = 0.1; 
  c1 = 0.1; 
  c2 = 0.1;
  Ls = 1000; 
  Ln = 1000;

  D1 = 2*rand( M1, K ) - 2;  D1 = normalize_D_n(D1);
  D2 = 2*rand( M2, K ) - 2;  D2 = normalize_D_n(D2);
  F1 = 2*rand( M1, K1 ) - 2; F1 = normalize_D_n(F1);
  F2 = 2*rand( M2, K2 ) - 2; F2 = normalize_D_n(F2);
  
  D = [ D1 F1 zeros(M1,K2) ; D2 zeros(M2,K1) F2 ];
  
  b0 = (K*c0) / Ls; b1 = (K1*c1) / Ln; b2 = (K2*c2) / Ln;
  c0 = c0*ones(K,1);  b0 = b0*ones(K,1);
  c1 = c1*ones(K1,1); b1 = b1*ones(K1,1);
  c2 = c2*ones(K2,1); b2 = b2*ones(K2,1);
  
  params = struct('D1',D1,'D2',D2,'F1',F1,'F2',F2,'D',D,...
                  'c0',c0,'b0',b0,'c1',c1,'b1',b1,'c2',c2,'b2',b2,...
                  'K1',K1,'K2',K2,'K',K,'M1',M1,'M2',M2,'type','linear');
end

% normalize every column of D to unit l1-norm
function D = normalize_D_n(Dest)
  K = size( Dest, 2 );
  D = zeros( size( Dest ) );
  for k = 1:K
    D(:,k) = Dest(:,k) / norm( Dest(:,k), 1 );
  end
end
