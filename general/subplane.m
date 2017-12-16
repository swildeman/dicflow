function h = subplane( h, wgt )
%SUBPLANE fit and subtract a 2d plane from a matrix representing elevation
%values
% 
% SYNOPSIS: h = subplane( h, wgt )
%
% INPUT h: original height field
%       wgt: (optional) logical (0,1) weight matrix of size(h)
%
% OUTPUT h: height field with fitted plane subtracted from it
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

if nargin < 2
    wgt = true(size(h));
end

[rows, cols] = size(h);
[r, c] = ndgrid(1:rows, 1:cols);

A = [ones(rows*cols, 1), r(:), c(:)];

fitc = A(wgt(:),:) \ h(wgt);

h(:) = h(:) - A*fitc;

end

