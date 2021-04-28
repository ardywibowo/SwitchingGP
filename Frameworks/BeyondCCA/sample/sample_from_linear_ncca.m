function [X1, X2, Alph, Beta1, Beta2] = sample_from_linear_ncca(params,N)
% Samples from model
%  aa, bb1, bb2 are i.i.d. Rademacher r.v. (1 or -1 with probability 0.5)
%  beta1 ~ bb1*Gamma( c1,b1 )
%  beta2 ~ bb2*Gamma( c2,b2 )
%  alpha ~ aa*Gamma( c0, b0 ) 
%     x1 = D1*alpha + F1*beta1 
%     x2 = D2*alpha + F2*beta2 
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  [X1,X2,Alph,Beta1,Beta2] = initialize(params,N);
  
  % sample in batches
  n = 1000;
  times = floor(N/n);
  rest = mod(N,n);
  
  for i=1:times
    inds = (i-1)*n + 1:i*n;
    [x1,x2,alph,beta1,beta2] = sample_batch(params,n);
    X1(:,inds) = x1;
    X2(:,inds) = x2;
    Alph(:,inds) = alph;
    Beta1(:,inds) = beta1;
    Beta2(:,inds) = beta2;
  end
  
  inds = times*n + 1 : times*n + rest;
  [x1,x2,alph,beta1,beta2] = sample_batch(params,rest);
  X1(:,inds) = x1;
  X2(:,inds) = x2;
  Alph(:,inds) = alph;
  Beta1(:,inds) = beta1;
  Beta2(:,inds) = beta2;
    
end

function [x1,x2,alph,beta1,beta2] = sample_batch(params,n)
  K1 = params.K1; K2 = params.K2; K = params.K;
  aa  = binornd( 1, 0.5, K, n );  aa( logical( aa==0 ) ) = -1;
  bb1 = binornd( 1, 0.5, K1, n ); bb1( logical( bb1==0 ) ) = -1;
  bb2 = binornd( 1, 0.5, K2, n ); bb2( logical( bb2==0 ) ) = -1;
  alph  = aa.*gamrnd( repmat(params.c0,1,n), repmat(1./params.b0,1,n) );
  beta1 = bb1.*gamrnd( repmat(params.c1,1,n), repmat(1./params.b1,1,n) );
  beta2 = bb2.*gamrnd( repmat(params.c2,1,n), repmat(1./params.b2,1,n) );
  if params.K1 > 0
    x1 = sparse( (params.D1*alph + params.F1*beta1) );
  else
    x1 = sparse( (params.D1*alph ) );
  end
  if params.K2 > 0
    x2 = sparse( (params.D2*alph + params.F2*beta2) );
  else
    x2 = sparse( (params.D2*alph) );
  end
end

function [X1,X2,Alph,Beta1,Beta2] = initialize(params,N)

  D1 = params.D1; D2 = params.D2;
  F1 = params.F1; F2 = params.F2;
  
  M1 = size(D1,1);
  M2 = size(D2,1);
  K0 = size(D1,2);
  K1 = size(F1,2);
  K2 = size(F2,2);
  
  X1 = sparse(M1,N);
  X2 = sparse(M2,N);
  Alph  = zeros(K0,N);
  Beta1 = zeros(K1,N);
  Beta2 = zeros(K2,N);
  
end
