function [ u, v, vals ] = corr2disp( xc, mincorr, subpixel )
%CORR2DISP Find peaks in a set of correlation maps and for each map return 
%the displacement of the peak from the center
%
% SYNOPSIS: [ u, v, vals ] = corr2disp( xc, mincorr, subpixel )
%
% INPUT xc: set of N correlation maps (as m*m*N matrix)
%       mincorr: minimum correlation peak value to consider valid
%       subpixel: (default: true) subpixel gaussian peak fitting?
%
% OUTPUT u,v: N displacement vectors (in x (left-right) and y (top-bottom) 
%             direction
%        vals: corresponding peak values
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

if nargin < 3
   subpixel = true; 
end

[xrows, xcols, xpages] = size(xc);
cnt_row = ceil(xrows/2);
cnt_col = ceil(xcols/2);

[vals , inds] = max( reshape(xc, [], xpages) );
[r, c] = ind2sub([xrows,xcols], inds);

if(subpixel)
    % make sure there are no negative values, as we are going to take a
    % logaritm
    e = eps('double');
%     xc = bsxfun(@minus, xc, min(min(xc-e)));
    xc(xc<=0) = e;
    
    % exlude points on border in subpixel peakfinding
    on_lr = ((c == xcols) | (c == 1));
    on_tb = ((r == xrows) | (r == 1));
    
    % get neighbor indices
    p = 1:xpages;
    
    szx = size(xc);
    c_ind = sub2ind(szx, r, c, p);  % center
    l_ind = sub2ind(szx, r(~on_lr), c(~on_lr)-1, p(~on_lr)); % left
    r_ind = sub2ind(szx, r(~on_lr), c(~on_lr)+1, p(~on_lr)); % right
    t_ind = sub2ind(szx, r(~on_tb)-1, c(~on_tb), p(~on_tb)); % top
    b_ind = sub2ind(szx, r(~on_tb)+1, c(~on_tb), p(~on_tb)); % bottom
    
    % get neighbor values
    cv = xc(c_ind);
    lv = xc(l_ind);
    rv = xc(r_ind);
    tv = xc(t_ind);
    bv = xc(b_ind);
      
    % fit a gaussian
    dc = zeros(size(c), 'like', c);
    dr = zeros(size(r), 'like', r);
    dc(~on_lr) = subpxpeak(cv(~on_lr), lv, rv);
    dr(~on_tb) = subpxpeak(cv(~on_tb), tv, bv);
    
    % just in case the "subpixel" peak finding has failed...
    dc(abs(dc) > 1) = 0;
    dr(abs(dr) > 1) = 0;
    
    c = c + dc;
    r = r + dr;
end

u = c(:)-cnt_col;
v = r(:)-cnt_row;

% make sure we do not include bad correlations
notvalid = vals < mincorr;
u(notvalid) = NaN;
v(notvalid) = NaN;

function di = subpxpeak(iv, lv, rv)
    % assume gaussian peak
    di = - .5*( (reallog(rv)-reallog(lv))./(reallog(lv)+reallog(rv)-2*reallog(iv)) );
end

end