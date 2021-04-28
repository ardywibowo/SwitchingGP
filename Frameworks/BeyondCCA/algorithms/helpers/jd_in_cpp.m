function [V, diags] = jd_in_cpp(B)
%JD_IN_CPP Orthogonal joint diagonalization, calls the respective 
% c++ function
%
% [V, diags] = jd_in_cpp(B)
% 
% Input: 
%     B      : m-by-(n*m) matrix of n m-by-m matrices to be jointly
%              diagonalized
%
% Output: 
%     V     : orthogonal matrix, the output of the joint diagonalization
%     diags : cell array of (almost) diagonal matrices
%
% Comment: This function for calling the C++/MEX-Matlab function
%   "jd_in.cpp". See its comments for details.
%
%
% Copyright 2016, Anastasia Podosinnikova

  [m, nm] = size(B);
  n = nm/m; if m*n~=nm, error('Wrong input size'), end

  % set parameters
  eps = 1e-8;
  V0 = eye(m); % initialization is only with the idenity for now..
  
  % run jd_in.cpp
  [out1, out2] = jd_in( B(:), m, n, eps, V0(:)); 
  
  [V, diags] = reshape_the_output_of_jd_in_cpp(out1, out2, m, n);
  
end

function [V, diags] = reshape_the_output_of_jd_in_cpp(out1, out2, m, n)
% reshape the output of jd_in.cpp
  V = reshape(out1, m, m);
  diags_temp = reshape(out2, m, m*n);
  diags = cell(1,n);
  for i = 1:n, 
    diags{i} = diags_temp(:, (i-1)*m+1:i*m); 
  end
end
