function [ D ] = designgrad1D( n )
%DESIGNGRAD1D Generate nxn 2nd order accurate finite difference matrix D so
%that fy = D*f and fx = f*D.', where f can be matrix
% 
% SYNOPSIS: [ D ] = designgrad1D( n )
%
% INPUT n: size of vector or matrix in the direction the derivative 
%          is to be taken
%
% OUTPUT D: finite difference design matrix for taking first derivative
%
% See also:
% DESIGNGRAD
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

% interior
D = spdiags(repmat([-1,0,1]/2,[n,1]),[-1, 0, 1],n,n);

% 2nd order forward on boundaries
D(1,1:3) = [-3,4,-1]/2;
D(end,end-2:end) = -[-1,4,-3]/2;

% 1st order forward on boundaries
% D(1,1:2) = [-1,1];
% D(end,end-1:end) = [-1,1];

end

