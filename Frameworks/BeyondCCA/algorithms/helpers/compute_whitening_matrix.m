function [W, VK, DK] = compute_whitening_matrix(S, K)
%COMPUTE_WHITENING_MATRIX Computes a whitening matrix W of S
%
% [W, VK, DK] = compute_whitening_matrix(S, K)
%
% INPUT:
%   S  : symmetric M-by-M real matrix (DICA S-cumulant or LDA S-moment)
%   K  : number of topics
%
% OUTPUT:
%   W  : K-by-M whitening matrix of S
%   VK : first K largest eigenvalues of S
%   DK : eigenvectors that correspond to VK
%
% COMMENTS: Computes a K-by-M whitening matrix W of a symmetric M-by-M
%   real matrix S, i.e. W*S*W'=I where I is the K-identity. M is the number 
%   of words in the dictionary. 
%
%
% Copyright 2016, Anastasia Podosinnikova

  verify_correctness_of_the_input(S, K)
  
  M = size(S,1);
  
  % S = 0.5*(S + S') to avoid complex values
  if M < 500,
    % matrix S is small
    [V, D] = eig(0.5*(S + S'));
  else
    % matrix S is larger
    [V, D] = eigs(0.5*(S + S'), K + 100);
  end
  
  [d,inds] = sort(diag(D),'descend'); 
  DK  = d(1:K); 
  if sum(DK<0) > 0, error('At whitening: negative eigenvalues'), end
  VK  = V(:, inds(1:K));
  W   = diag(DK.^(-1/2))*VK';
  
end

function verify_correctness_of_the_input(S, K)

  [n, m] = size(S);
  if n ~= m, error('S have to be a square matrix'), end
  if K > n, error('K can not exceed M'), end
  
  if n > 15000,
    warning(strcat('This function is not designed for large M. It i', ...
                   's highly recomended to terminate execution. If c', ...
                   'ontinued, it will run eigs of M-by-M matrix.'))
    STR = input('Do you want to terminate execution? (y/n)','s');
    if strcat(STR,'y')
      error('Execution of the program was terminated.')
    end
  end
  
end
