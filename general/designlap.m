function [ Dxx, Dyy ] = designlap( roi, varargin )
%DESIGNLAP Generate matrix for taking the laplacian of an image
%
% SYNOPSIS: [ Dxx, Dyy ] = designlap( roi, bc )
%
% INPUT roi: Size of the target image or ROI matrix for target image 
%            (where 1's correspond to the 'valid' region), be carefull that
%            the ROI does not leave any <3 px gaps...
%       bc: Boundary condition: 'neumann' or 'forward'
% 
% OUTPUT Dxx, Dyy: Design matrices for taking laplacian of image I unrolled 
%        as I(:) (or I(roi)).
% 
% EXAMPLE   For an image I, with roi R one can use the following to
%           compute the gradients of the image:
%               [Dxx, Dyy] = designLap(R);
%               IDxx = zeros(size(I)); IDyy = zeros(size(I));
%               IDxx(R) = Dxx*I(R); IDyy(R) = Dyy*I(R);
%               imshowpair(IDxx,IDyy,'montage');
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

allowedBC = {'neumann','forward'};
defaultBC = 'neumann';
p = inputParser;
p.addRequired('roi',@(x) numel(x) > 1);
p.addParameter('bc',defaultBC,...
    @(x) any(validatestring(x,allowedBC)));

p.parse(roi,varargin{:});
in = p.Results;

bc = validatestring(in.bc, allowedBC);

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

% treat outer edges as a masked border so that same method applies
roi = padarray(roi,[1,1],false);

[IRows,ICols] = size(roi);

% total number of valid pixels at which to calculate derivatives
N = sum(roi(:));

% linear index in masked image (skipping mask positions in numbering)
linInd = zeros(IRows,ICols,create_cls);
linInd(roi) = 1:N;

% init vectors for holding the row and column positions of design matrix entries
% (3 positions for each pixel)
rowDM = zeros(N,3,create_cls);
colDM = zeros(N,3,create_cls);
valDM = zeros(N,3,create_cls);

% keep track of possition in the entry lists we are filling
cur = 1;

% compute horizontal gradient to detect edges in mask
roiD = diff(roi,1,2);

% find left sided edges (0 0 0 1 1 1)
pxind = find([zeros(IRows,1), roiD] == 1);
npx = length(pxind);
% fill entries with foward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1, :) = [linInd(pxind), linInd(pxind+IRows), linInd(pxind+2*IRows)]; % edge pixel locations in unrolled (masked) image
if strcmp(bc,'neumann')
    valDM(cur:cur+npx-1,:) = repmat([-2 2 0],npx,1);
else % forward
    valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);
end
cur = cur+npx;

% find right sided edges
pxind = find([roiD, zeros(IRows,1)] == -1);
npx = length(pxind);
% fill entries with backward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-2*IRows),linInd(pxind-IRows), linInd(pxind)];    % edge pixel locations in unrolled (masked) image
if strcmp(bc,'neumann')
    valDM(cur:cur+npx-1,:) = repmat([0 2 -2],npx,1);
else % forward
    valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);
end
cur = cur+npx;

% fill interior with central difference
pxind = find(~imdilate(~roi, [1, 1, 1]));
npx = length(pxind);
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind), linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-IRows),linInd(pxind), linInd(pxind+IRows)];
valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);

Dxx = sparse(rowDM,colDM,valDM);

% Do the same for vertical gradient, with substitution (and similar):
% [zeros(IRows,1), maskD] -> [zeros(1,ICols); maskD]
% linInd(pxind+IRows) -> linInd(pxind+1) 

% reset pivot (reuse same matrices for sake of memory efficiency)
cur = 1;

% compute vertical gradient to detect edges in mask
roiD = diff(roi,1,1);

% find top sided edges (0 0 0 1 1 1)
pxind = find([zeros(1,ICols); roiD] == 1);
npx = length(pxind);
% fill entries with foward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind),linInd(pxind)];
colDM(cur:cur+npx-1, :) = [linInd(pxind), linInd(pxind+1), linInd(pxind+2)];    % edge pixel locations in unrolled (masked) image
if strcmp(bc,'neumann')
    valDM(cur:cur+npx-1,:) = repmat([-2 2 0],npx,1);
else
    valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);
end
cur = cur+npx;

% find bottom sided edges
pxind = find([roiD; zeros(1,ICols)] == -1);
npx = length(pxind);
% fill entries with backward difference
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind),linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-2), linInd(pxind-1), linInd(pxind)];    % edge pixel locations in unrolled (masked) image
if strcmp(bc,'neumann')
    valDM(cur:cur+npx-1,:) = repmat([0 2 -2],npx,1);
else
    valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);
end
cur = cur+npx;

% fill interior with central difference
pxind = find(~imdilate(~roi, [1; 1; 1]));
npx = length(pxind);
rowDM(cur:cur+npx-1,:) = [linInd(pxind), linInd(pxind),linInd(pxind)];
colDM(cur:cur+npx-1,:) = [linInd(pxind-1), linInd(pxind), linInd(pxind+1)];
valDM(cur:cur+npx-1,:) = repmat([1 -2 1],npx,1);

Dyy = sparse(rowDM,colDM,valDM);

end

