function wsw12 = compute_gen_wsw12_cca_all(...
  X1, X2, W1, W2, t1, t2, N, type )

  if ~( strcmp( type, 'dcca' ) || strcmp( type, 'ncca' ) || strcmp( type, 'mcca' ) ), error('Wrong type'), end
  if sum( t1==0 ) > 0 || sum( t2==0 ) > 0, error('t1 or t2 must not have zeros'), end

  W1 = sparse(W1);
  W2 = sparse(W2);
  
  M1 = size(X1,1); M2 = size(X2,1);
  
  if strcmp( type, 'dcca' )
    XW1 = W1*sparse( 1:M1, 1:M1, exp(t1).^(-1) )*X1;
    XW2 = W2*sparse( 1:M2, 1:M2, exp(t2).^(-1) )*X2;
  end
  
  if strcmp( type, 'mcca' )
    XW1 = W1*X1;
    XW2 = W2*sparse( 1:M2, 1:M2, exp(t2).^(-1) )*X2;
  end
  
  if strcmp( type, 'ncca' )
    XW1 = W1*X1;
    XW2 = W2*X2;
  end
  
  proj  = X1'*t1 + X2'*t2;
  eproj = exp(proj);
  
  wsw12 = ( XW1*sparse(1:N,1:N,eproj)*XW2' ) / sum(eproj) ...
        - ( (XW1*eproj) / sum(eproj) )*( (XW2*eproj) / sum(eproj) )';
  
end