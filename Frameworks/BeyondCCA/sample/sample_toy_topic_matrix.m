function D = sample_toy_topic_matrix(M, K, a, iswellcond)
%SAMPLE_TOY_TOPIC_MATRIX Sampels an M-by-K matrix D such that each
% column of this matrix is from Dirichlet(a). If iswllcond is 1, the 
% condition number (largest/smallest singular value) does not exceed 15.
% When a is 1, each column is sampled uniformly in the simplex.
%
%
% Copyright 2016, Anastasia Podosinnikova

  basemeasure = ones(M, 1);
  param = a*basemeasure;
  encore = 1; iter = 1;
  while encore
    temp = gamrnd(repmat(param,1,K), 1, M, K);
    r = sum(temp, 1);
    r(logical(r == 0)) = 1;
    D = temp./repmat(r, M, 1);
    sss = svd(D);
    % ensures that the sampled matrix is well conditioned
    if sss(1)/sss(K) <= 15, encore = 0; end
    iter = iter + 1;
    if iswellcond == 1
      if iter > 100, 
        disp(['Condition number: ', num2str( sss(1)/sss(K) )])
        error('Can not sample a well-conditioned matrix'), 
      end
    else
      break;
    end
  end

end
