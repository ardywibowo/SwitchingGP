function params = sample_parameters_for_figure_2_right
%SAMPLE_PARAMETERS_FOR_FIGURE_2_LEFT Samples parameters for the experiment 
% from Figure 2 (right) of the paper.
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  % see the paper for the notation
  M1 = 20;
  M2 = 20;
  K1 = 20; 
  K2 = 20; 
  K = 10;
  ad = 0.5; % each column of D1(2) is sampled as Dirichlet(ad)
  af = 0.5; % each column of F1(2) is sampled as Dirichlet(af)
  c0 = 0.3; 
  c1 = 0.1; 
  c2 = 0.1;
  Ls = 1000; % expected document length (due to signal)
  Ln = 1000; % expected document legnth (due to noise)

  iswellcond = 0;
  D1 = sample_toy_topic_matrix(M1, K, ad, iswellcond);
  D2 = sample_toy_topic_matrix(M2, K, ad, iswellcond);
  F1 = [];
  if K1 > 0, F1 = sample_toy_topic_matrix(M1, K1 ,af, iswellcond); end
  F2 = [];
  if K2 > 0, F2 = sample_toy_topic_matrix(M2, K2, af, iswellcond); end
  
  D = [ D1 F1 zeros(M1,K2) ; D2 zeros(M2,K1) F2 ];
  
  b0 = (K*c0) / Ls; b1 = (K1*c1) / Ln; b2 = (K2*c2) / Ln;
  c0 = c0*ones(K,1);  b0 = b0*ones(K,1);
  c1 = c1*ones(K1,1); b1 = b1*ones(K1,1);
  c2 = c2*ones(K2,1); b2 = b2*ones(K2,1);
  
  params = struct('D1',D1,'D2',D2,'F1',F1,'F2',F2,'D',D,...
                  'c0',c0,'b0',b0,'c1',c1,'b1',b1,'c2',c2,'b2',b2,...
                  'K1',K1,'K2',K2,'K',K,'M1',M1,'M2',M2,'type','linear');
end
