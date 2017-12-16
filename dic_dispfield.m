function [ u, v, br, bc, corr_qual ] = dic_dispfield( Iref, Idef, blockspacing, borderoverlap, maxdisp, refinedisp, minqual, normcorr, wgt )
%DIC_DISPFIELD Use digitial image correlation (DIC) to find the displacement
%field that warps Iref into Idef.
%
% SYNOPSIS: [u,v,br,bc,corr_qual] = dic_dispfield(Iref,Idef,blockspacing,borderoverlap,maxdisp,...)
%
% INPUT Iref: Reference image, to be sliced up in template blocks
%       Idef: Distorted image, to be sliced up in search window blocks
%       blockspacing: Grid spacing between consecutive points for which 
%                     displacement vectors are found
%       borderoverlap: Number of pixels each template block extends (on all 
%                      sides) beyond the grid spacing (block overlap) 
%                      (template_width = blockspacing + 2*borderoverlap)
%       maxdisp: maximum vert/hor displacement within which to search for 
%                a match between each template and search window
%                (search window width = template_width + 2*maxdisp)
%       refinedisp: (optional, for warping schemes) replace maxdisp by 
%                   refinedisp < maxdisp, but keep the block coordinates 
%                   (offset from top and left image border) as they were, 
%                   set to [] if not needed.
%       minqual: (default: 0.01) vectors with a correlation quality lower 
%                than minqual will be disregarded (set to NaN). 
%       normcorr: (default: true) use normalized cross correlation 
%                 (if set to false, basic cross correlation is used)
%       wgt: (optional) weight matrix of size(Iref) to take into account in
%            the normalized cross corelation (e.g. for masked regions)
%            (Padfield 2010, Masked FFT Registration, IEEE CVPR Conference)
%
% OUTPUT u: displacement in x (along columns, left to right) direction
%        v: displacement in y (along rows, top to bottom) direction
%        br: y (row) coord. of (u,v) vector (centered on template blocks)
%        bc: x (col) coord. of (u,v) vector (centered on template blocks)
%        corr_qual: correlation quality for each vector
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

if nargin < 9
    weighted = false;
else
    weighted = true;
    wgt = double(wgt);
end

if nargin < 8
    normcorr = true;
end

if nargin < 7
   minqual = 0.01; 
end

if nargin < 6 || isempty(refinedisp) || refinedisp > maxdisp || refinedisp < 1
    refinedisp = maxdisp;
end

% in case of simple correlation, just apply the mask directly to the images
if weighted && ~normcorr
    Iref = Iref.*wgt;
    Idef = Idef.*wgt;
end

% make sure images are doubles between 0 and 1 and increase the contrast so
% that normswcorr2 detects enough correlation energy.
Iref = imstretchlim(Iref);
Idef = imstretchlim(Idef);

[rows, cols] = size(Iref);

% number of blocks in y and x direction
bRows = floor((rows-2*(maxdisp+borderoverlap)) / blockspacing);
bCols = floor((cols-2*(maxdisp+borderoverlap)) / blockspacing);

% slice up Iref into template blocks
[t, br, bc] = im2blocks(Iref, bRows, bCols, blockspacing, borderoverlap, maxdisp+borderoverlap);

% slice up Idef into search window blocks
sw = im2blocks(Idef, bRows, bCols, blockspacing, refinedisp+borderoverlap, maxdisp+borderoverlap);

if normcorr
    % one massive normalized cross correlation using fft2
    if weighted
        t_wgt = im2blocks(wgt, bRows, bCols, blockspacing, borderoverlap, maxdisp+borderoverlap);
        sw_wgt = im2blocks(wgt, bRows, bCols, blockspacing, refinedisp+borderoverlap, maxdisp+borderoverlap);
        xc = normswcorr2(t, sw, [], t_wgt, sw_wgt);
    else
        xc = normswcorr2(t, sw);
    end
else
   % subtract mean from each block
   sw = bsxfun(@minus, sw, mean(mean(sw)));
   t = bsxfun(@minus, t, mean(mean(t)));
   
   % one massive cross correlation using fft2
   xc = swcorr2(t, sw);
   
   % get an idea of the quality of the correlation
   maxcorr = .25*(blockspacing+2*maxdisp).^2;
   xc = xc / maxcorr;
end

% find subpixel displacements from the position of the correlation peaks
[u, v, corr_qual] = corr2disp(xc, minqual);

% get the arrays back in shape :)
u = reshape(u, [bRows, bCols]);
v = reshape(v, [bRows, bCols]);
br = reshape(br, [bRows, bCols]);
bc = reshape(bc, [bRows, bCols]);
corr_qual = reshape(corr_qual, [bRows, bCols]);

end

