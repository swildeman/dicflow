function [ Dx, Dy ] = designgrad( roi )
%DESIGNGRAD Generate matrix for taking central difference in x,y-direction
% 
% SYNOPSIS: [ Dx, Dy ] = designgrad( roi )
%
% INPUT roi: size of the target image, or ROI matrix for target image 
%            (where 1's correspond to the 'valid' region), be carefull that
%            the roi does not leave any <3 px gaps...
%
% OUTPUT Dx, Dy: Design matrices for taking (x,y)-derivative of image I
%                unrolled as I(roi). It uses central differences in interior 
%                and forward/backward difference at (roi) boundaries.
%   
% EXAMPLES   For an image I, with ROI R one can use the following to
%            compute the gradients of the image:
%                 [Dx, Dy] = designGrad(R);
%                 IDx = zeros(size(I)); IDy = zeros(size(I));
%                 IDx(R) = Dx*I(R); IDy(R) = Dy*I(R);
%                 imshowpair(IDx,IDy,'montage');
%
%             One can also use Dx and Dy to inverse the gradient and
%             recover I (up to a constant, taken as I_11 = 0):
%                 IInt = zeros(size(IDx));
%                 IInt(R)= [0;[Dx(:,2:end);Dy(:,2:end)]\[IDx(R);IDy(R)]];
%                 imshowpair(I,IInt,'montage');
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

% parse input
if numel(roi) < 2
   error('Expected the mask variable to have at least 2 elements.');
end

if numel(roi) == 2
    roi = true(roi); % assume mask represents the size of the image
else
    roi = logical(roi);
end

if isa(roi, 'gpuArray')
    create_cls = 'gpuArray';
else
    create_cls = 'double';
end

% treat border same as masked areas
roi = padarray(roi,[1,1],false);

[IRows,ICols] = size(roi);

% total number of valid pixels at which to calculate derivatives
N = sum(roi(:));

% linear index in masked image (skipping mask positions in numbering)
linInd = zeros(IRows,ICols,create_cls);
linInd(roi) = 1:N;

% init vectors for holding the row and column positions of design matrix entries
% (2 positions for each pixel)
rowDM = zeros(N,3,create_cls);
colDM = zeros(N,3,create_cls);
valDM = zeros(N,3,create_cls);

% keep track of possition in the entry lists we are filling
cur = 1;

% compute horizontal gradient to detect edges in roi
roiD = diff(roi,1,2);

% find left sided edges (0 0 0 1 1 1)
pxind = find([zeros(IRows,1), roiD] == 1);
npx = length(pxind);
% fill entries with 2nd order foward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1, :) = [linInd(pxind), linInd(pxind+IRows),linInd(pxind+2*IRows)];    % edge pixel locations in unrolled (masked) image
valDM(cur:cur+npx-1,:) = repmat([-3,4,-1]/2,npx,1);
cur = cur+npx;

% find right sided edges
pxind = find([roiD, zeros(IRows,1)] == -1);
npx = length(pxind);
% fill entries with backward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-2*IRows),linInd(pxind-IRows), linInd(pxind)];    % edge pixel locations in unrolled (masked) image
valDM(cur:cur+npx-1,:) = repmat(-[-1,4,-3]/2,npx,1);
cur = cur+npx;

% fill interior with central difference
pxind = find(~imdilate(~roi, [1, 1, 1]));
npx = length(pxind);
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind),linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-IRows),linInd(pxind), linInd(pxind+IRows)];
valDM(cur:cur+npx-1,:) = repmat([-.5 0 .5],npx,1);

Dx = sparse(rowDM,colDM,valDM);

% Do the same for vertical gradient, with substitution (and similar):
% [zeros(IRows,1), maskD] -> [zeros(1,ICols); maskD]
% linInd(pxind+IRows) -> linInd(pxind+1) 

% reset pivot (reuse same matrices for sake of memory efficiency)
cur = 1;

% compute vertical gradient to detect edges in mask
roiD = diff(roi,1,1);

% find top sided edges (0 0 0 1 1 1...)
pxind = find([zeros(1,ICols); roiD] == 1);
npx = length(pxind);
% fill entries with foward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1, :) = [linInd(pxind), linInd(pxind+1), linInd(pxind+2)];    % edge pixel locations in unrolled (masked) image
valDM(cur:cur+npx-1,:) = repmat([-3,4,-1]/2,npx,1);
cur = cur+npx;

% find bottom sided edges
pxind = find([roiD; zeros(1,ICols)] == -1);
npx = length(pxind);
% fill entries with backward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-2), linInd(pxind-1), linInd(pxind)];    % edge pixel locations in unrolled (masked) image
valDM(cur:cur+npx-1,:) = repmat(-[-1,4,-3]/2,npx,1);
cur = cur+npx;

% fill interior with central difference
pxind = find(~imdilate(~roi, [1; 1; 1]));
npx = length(pxind);
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-1), linInd(pxind), linInd(pxind+1)];
valDM(cur:cur+npx-1,:) = repmat([-.5 0 .5],npx,1);

Dy = sparse(rowDM,colDM,valDM);

end

