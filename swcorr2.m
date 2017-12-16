function xc = swcorr2( t, sw )
%SWCORR2 fft-based, page-wise, cross corelation between a set of
%templates and search windows, with size(sw) > size(t)
%
% SYNOPSIS: xc = swcorr2( t, sw )
%
% INPUT t: set of N square templates of size m*m (as (m*m*N) matrix)
%       sw: set of N square search windows of size k*k (k > m)
%           (as (k*k*N) matrix)
%
% OUTPUT xc: cross correlations between templates and search windows
%            as (k-m)*(k-m)*N matrix (only valid points are returned)
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

[tRows, tCols, ~] = size(t);
[swRows, swCols, ~] = size(sw);

maxRowDisp = swRows - tRows;
maxColDisp = swCols - tCols;

if any(mod([maxRowDisp, maxColDisp],2)) || any([maxRowDisp, maxColDisp]<0)
   error('Search windows are expected to have symmetric (positive) padding around the templates in both dimensions'); 
else
   maxRowDisp = maxRowDisp/2;
   maxColDisp = maxColDisp/2;
end

f_sw = fft2(sw);
f_t = fft2(padarray(t,[maxRowDisp,maxColDisp]));
xc = fftshift(fftshift( ifft2( bsxfun(@times, f_sw, conj(f_t)) ), 1), 2);

% return only valid points
clpRows = (tRows-1)/2;
clpCols = (tCols-1)/2;
validRows = 1+ceil(clpRows) : swRows-floor(clpRows);
validCols = 1+ceil(clpCols) : swCols-floor(clpCols);

xc = real(xc(validRows, validCols,:));

end

