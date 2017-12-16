function A = imstretchlim( A )
%IMSTRETCHLIM normalize images so that values are between 0 and 1 (supports
%stacked images).
%
% SYNOPSIS: A = imstretchlim( A )
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

A = double(A);

mx = max(max(A,[],1),[],2);
mn = min(min(A,[],1),[],2);
delta = mx-mn;

A = bsxfun(@minus,A,mn);
A = bsxfun(@rdivide,A,delta);

end

