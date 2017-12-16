function I = particleimage( w, h, blocksize, dev, partsize )
%PARTICLEIMAGE Generate an irregular dot-pattern for DIC/OF applications
%
% SYNOPSIS: I = particleimage( w, h, blocksize, dev, partsize ). Generate
% an image of size w*h with a dot placed at random in sub blocks of size
% blocksize.
%
% INPUT w,h: width and height of the generated particle image
%       blocksize: size of the sub blocks in which to put a particle
%       dev: maximum deviation of particle placement (from center of
%            blocks)
%       partsize: size of particles (standard deviation of gaussian profile)
%
% OUTPUT I: generated particle image of size w*h
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

I = zeros(h,w);
cntrow = round(blocksize/2);
cntcol = cntrow;

I = blockproc(I,[blocksize,blocksize], @putrandomdot,'PadPartialBlocks', true);

I = imgaussfilt(I,partsize,'FilterSize', 4*ceil(2*partsize)+1);

I = 1 - I;

I = imadjust(I,stretchlim(I,0.005));

I = I(1:h, 1:w);

    function block = putrandomdot(blockstruct)
        row = cntrow + randi(2*dev+1)-dev-1;
        col = cntcol + randi(2*dev+1)-dev-1;
        block = blockstruct.data;
        block(row,col) = 1;
    end

end

