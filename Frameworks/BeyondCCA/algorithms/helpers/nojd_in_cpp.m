function [Q, Qinv] = nojd_in_cpp(B)
%NOJD_IN_CPP Non-orthogonal joint diagonalization by similarity, calls 
% the respective c++ function
%
% [Q, Qinv] = nojd_in_cpp(B);
% 
% Input:
%   B : the K-by-K*N matrix with N target matrices of size K-by-K
%
% Output:
%   Q    : the K-by-K diagonalizer matrix
%   Qinv : the inverce of the matrix Q
%
% Comment: This function is calling the C++/MEX-Matlab function
%   "nojd_in.cpp". See its comments for details.
%
%
% Copyright 2016, Anastasia Podosinnikova
  
  kmax = 1000; % max number of sweeps
  epsilon = 1E-4;
  [P, Pinv] = nojd_in(B, kmax, epsilon);
  Q = P';
  Qinv = Pinv';
  
end