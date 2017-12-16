function nxc = normswcorr2( t, sw, minEnergy, t_wgt, sw_wgt )
%NORMSWCORR2 Normalized (masked) cross correlation between set of square
%templates t and search windows sw, with size(sw) > size(t)
%
% SYNOPSIS: nxc = normswcorr2( t, sw, minEnergy, t_wgt, sw_wgt )
%
% INPUT t: set of N square templates of size m x m (as (m x m x N) matrix)
%       sw: set of N square search windows of size k x k (k > m)
%           (as (k x k x N) matrix)
%       minEnergy: minimum correlation energy considered as valid
%       t_wgt: (optional) Template weights to take into account for masked 
%              registration, (m x m x N) matrix
%       sw_wgt: (optional) Search window weights to take into account for 
%               masked registration, (k x k x N) matrix
%
% OUTPUT nxc: normalized cross correlation between t and sw as a 
%             (k-m) x (k-m) x N matrix (only valid points are returned)
%
% REMARKS algorithm from Padfield 2010, Masked FFT Registration, 
%                        IEEE CVPR Conference
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

[trows,tcols,~] = size(t);

if nargin < 4
   t_wgt = ones(size(t(:,:,1)));
   sw_wgt = ones(size(sw(:,:,1)));
end

if nargin < 3 || isempty(minEnergy)
   maxEnergy = .25*trows*tcols;
   minEnergy = 1e-2*maxEnergy;
end

% for numerator
t_s         = swcorr2(t, sw); % corr 1
twgt_s      = swcorr2(t_wgt, sw); % corr 2
t_swgt      = swcorr2(t, sw_wgt); % corr 3
twgt_swgt   = swcorr2(t_wgt, sw_wgt); % corr 4

% for denominator 1
twgt_s2     = swcorr2(t_wgt, sw.^2); % corr 5
twgt_s_2    = twgt_s.^2;
% twgt_swgt, already calculated above

% for denominator 2
t2_swgt     = swcorr2(t.^2, sw_wgt); % corr 6
t_swgt_2    = t_swgt.^2;

% prevent division by zero below
invalid = squeeze( ~all(all(twgt_swgt)) );
if any(invalid)
    twgt_swgt(:,:,invalid) = 1;
end

% note: use bsxfun for support of older matlab versions
num =   t_s - bsxfun(@rdivide, (twgt_s.*t_swgt),twgt_swgt );
den1 =  twgt_s2 - bsxfun(@rdivide, twgt_s_2, twgt_swgt);
den2 =  t2_swgt - bsxfun(@rdivide, t_swgt_2, twgt_swgt);
den = den1.*den2;
den(den<eps('double')) = Inf; % prevent complex sqrt and make sure we don't have a zero division later on;
den = realsqrt(den);

% check for shifts where there is no overlap of weights
if any(invalid)
    den(:,:,invalid) = Inf;
end

% check for areas where no good SNR is expected
toolow = squeeze(min(min(den))) < minEnergy;
den(:,:,toolow) = Inf;
% den(den < minEnergy) = Inf;

nxc = num ./ den;

end

