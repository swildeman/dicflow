function f = invgrad2( fx, fy, roi )
%INVGRAD2(fx, fy, roi) Integrate a gradient field (up to a constant) by 
%inverting a finite difference system in LSQ sense, with foward/backward 
%difference at ROI edges
% 
% SYNOPSIS: f = invgrad2( fx, fy, roi )
%
% INPUT fx: gradient accross columns (left to right)
%       fy: gradient accross rows (up to down)
%       roi: logical array of size(fx) with 1's at valid positions
%
% OUTPUT f: integrated gradient field such that [fx,fy] = grad(f) in a
%           least squares sense
%
% REMARK uses faster <a href="matlab:help invgrad2_rect">INVGRAD2_RECT</a> for rectangular domains
%
% See also:
% INVGRAD2_RECT
% FFTINVGRAD
% DESIGNGRAD
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

if nargin < 3
    roi = true(size(fx));
end

if sum(roi(:)) == numel(fx) % domain is whole image
    f = invgrad2_rect(fx, fy);
else % any other domain shape
    [Dx, Dy] = designgrad(roi);
    f = zeros(size(fx));

    f(roi)= [0;[Dx(:,2:end);Dy(:,2:end)] \ [fx(roi);fy(roi)]];
end

end