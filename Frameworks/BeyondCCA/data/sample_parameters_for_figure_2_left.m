function params = sample_parameters_for_figure_2_left
%SAMPLE_PARAMETERS_FOR_FIGURE_2_LEFT Returns (deterministic) parameters 
% for the experiment from Figure 2 (left and middle) of the paper.
%
%
% Copyright 2016, Anastasia Podosinnikova

  c0 = .1; 
  c1 = .1; 
  c2 = .1;
  Ls = 100; 
  Ln = 100;

  K = 1;
  K1 = 2;
  K2 = 2;
  M1 = 2;
  M2 = 2;

  D1 = [0.5; 0.5];
  F1 = [0.9 0.1; 0.1 0.9];
  D2 = [0.5; 0.5];
  F2 = [0.9 0.1; 0.1 0.9];

  D = [ D1 F1 zeros(M1,K2) ; D2 zeros(M2,K1) F2 ];

  b0 = (c0*K)  / Ls;
  b1 = (c1*K1) / Ln;
  b2 = (c2*K2) / Ln;
  c0 = c0*ones(K,1);  b0 = b0*ones(K,1);
  c1 = c1*ones(K1,1); b1 = b1*ones(K1,1);
  c2 = c2*ones(K2,1); b2 = b2*ones(K2,1);

  params = struct('D1',D1,'D2',D2,'F1',F1,'F2',F2,'D',D,...
                  'c0',c0,'b0',b0,'c1',c1,'b1',b1,'c2',c2,'b2',b2,...
                  'K1',K1,'K2',K2,'K',K,'M1',M1,'M2',M2,'type','linear');
end
