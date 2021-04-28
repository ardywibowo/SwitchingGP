function [err1, perm] = l1_error(Dest, D)
%L1_ERROR Computes the L1-error between two topic matrices Dest and D of
%   the same size; the matrices Dest and D are assumed to have 
%   non-negative entries and to be normalized to unit l1-norm for each
%   column
%
% [err1, rp] = l1_error(Dest, D)
%
% Input:
%   Dest : M-by-K topic matrix estimated by an algorithm
%   D    : M-by-K ground truth topic matrix
%
% Output:
%   err1 : the estimated error scaled to the [0,1] interval
%   perm : permutation of the columns of Dest corresponding to the
%          estimated error err1 (Dest(:,perm) \approx D) when properly
%          renormalized
%
%
% Copyright 2016, Anastasia Podosinnikova

  [Dest, D] = verify_correctness_of_the_input_and_normalize(Dest, D);
  
  Kest = size(Dest,2); K = size(D,2);
  Perf = zeros(Kest,K);
  for i=1:Kest, 
    for j=1:K, 
      Perf(i,j) = norm(Dest(:,i)-D(:,j), 1); % l1-norm!
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

  [M, ~] = size(D);
  
  if sum( sum( Dest < 0 ) ) > 0 || sum( sum( D < 0 ) ) > 0
    % check non-negativity
    error('Some values of D or Dest are negative')
  end
  
  if abs( sum(sum(Dest)) - size(Dest,2) ) > 1e-14*M || abs( sum(sum(D)) - size(D,2) ) > 1e-14*M
    % check the simplex constraint
    Dest = Dest ./ repmat( sum(Dest), M, 1 );
    D = D ./ repmat( sum(D), M, 1 );
  end

end
