function [I, bCntR, bCntC] = im2blocks( I, bRows, bCols, bSize, pad, bOffset )
%IM2BLOCKS Slice up an image into bRows x bCols square sub images
%(starting from the top-left corner)
%
% SYNOPSIS: [I, bCntR, bCntC] = im2blocks( I, bRows, bCols, bSize, pad, bOffset )
%
% INPUT I: image to slice up
%       bRows: number of blocks in vertical direction
%       bCols: number of blocks in horizontal direction
%       bSize: base block size of sub-image (width,height) = (bSize,bSize)
%       pad: (optional) additional overlap of sub image borders
%            i.e. (w,h) = (bSize+2*pad,bSize+2*pad)
%       bOffset: (optional) additional overall offset from top left image
%                corner
%
% OUTPUT I: sliced up version of I of dimension 
%           (m,n,p) = (bSize+2*pad, bSize+2*pad, bRows*bCols)
%        bCntR: y (row) coordinates of sub image centers
%        bCntC: x (col) coordinates of sub image centers
%
% REMARKS: Blocks extending outside the borders of image will be padded with zeros
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

if nargin < 5
    pad = 0;
    bOffset = 0;
elseif nargin < 6
    bOffset = 0;
end

% make sure padded blocks are centered around blocks
if bOffset < pad
    prepad = pad-bOffset;
    I = padarray(I, [prepad, prepad], 'pre');
else
    prepad = 0;
end

% make sure the requred nummer of blocks fit in the image
reqIsize = [bRows, bCols]*bSize + 2*pad + bOffset;
I = padarray(I, max([0,0], reqIsize-size(I)), 'post');

% top left row and column index of the blocks to be padded
[bTLRows, bTLCols] = ndgrid(1+prepad+bOffset+(0:bRows-1)*bSize, ...
                            1+prepad+bOffset+(0:bCols-1)*bSize);

% relative row and column indices
[pbRelRows, pbRelCols] = ndgrid(-pad:bSize+pad-1);

% convert to indices
[Irows, ~] = size(I);
bTLInd = bTLRows + Irows*(bTLCols-1);
pbRelInd = pbRelRows + Irows*pbRelCols;

pbInd = bsxfun(@plus, pbRelInd, shiftdim(bTLInd(:),-2));

I = I(pbInd);

bCntR = bTLRows(:) + (bSize-1)/2;
bCntC = bTLCols(:) + (bSize-1)/2;

end

