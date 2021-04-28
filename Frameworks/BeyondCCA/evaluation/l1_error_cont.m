function [err1, perm] = l1_error_cont(Dest, D)
%L1_ERROR_CONT Computes l1-error between two mixing matrices Dest and D,
% no constraints (such as non-negativity of the entries) are made
%
%
% Copyright 2016, Anastasia Podosinnikova

  [Dest, D] = verify_correctness_of_the_input_and_normalize(Dest, D);
  
  Kest = size(Dest,2); K = size(D,2);
  Perf = zeros(Kest,K);
  for i=1:Kest, 
    for j=1:K, 
      Perf(i,j) = norm( Dest(:,i)-D(:,j), 1 ); % l1-norm!
    end
  end
  
  if sum( sum( isnan( Perf ) ) ) > 0
    error( ' L1 error: Perf has NaNs ' )
  end

  % as the columns of the topic matrix are estimated up to permutations
  [matching, cost] = HungarianBipartiteMatching(Perf);
  
  [perm, ~] = find( sparse(matching) );
  
  err1 = cost/K;
  err1 = err1/2; % to scale in [0,1]
  
end

function [Dest, D] = verify_correctness_of_the_input_and_normalize(Dest, D)

  if sum( sum( isnan( Dest ) ) ) > 0 || sum( sum( isnan( D ) ) ) > 0
    error('l1error : Dest or D has NaNs')
  end
  
  if (size(Dest,1) ~= size(D,1))
    error('Wrong input size')
  end
  
  Dest = normalize_D_n(Dest);
  Dest = [Dest -Dest];
  D = normalize_D_n(D);
  
end
