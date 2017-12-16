function [I, u, v] = interpimwarp( I, u, v, x, y, interp_method )
%INTERPIMWARP Warp an image using Matlab's interp2
%
% SYNOPSIS: [I, u, v] = interpimwarp( I0, u, v, x, y, interp_method )
%
% INPUT I0: image to warp
%       u,v: displacement field used to warp image
%       x,y: column, row, coordinates of displacement vectors 
%       interp_method: interpolation method used for warping 
%                      'linear', 'nearest', 'cubic', or 'spline'
% 
% OUTPUT I: warped image I(x,y) = I0(x-u,y-v)
%        u,v: displacement field used for warping 
%             (original (u,v) interpolated to resolution of image) 
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

% make sure image is doubles between 0 and 1 in calculations
if isa(I, 'uint8') || (isa(I, 'gpuArray') && strcmp(classUnderlying(I),'uint8'))
    I = double(I)/255.;
    orig_uint8 = true;
else
    orig_uint8 = false;
end

% interpolate displacement field to match size of image
[xq, yq] = meshgrid(1:size(I,2), 1:size(I,1));
u = interp2(x,y,u,xq,yq,'linear',NaN);
v = interp2(x,y,v,xq,yq,'linear',NaN);

% warp coordinates
Xw = xq - u;
Yw = yq - v;

% don't warp invalid points
unan = isnan(u);
vnan = isnan(v);
Xw(unan) = xq(unan);
Yw(vnan) = yq(vnan);

% warp
I = interp2(xq, yq, I, Xw, Yw, interp_method, 0);

% restore original class if necessary
if orig_uint8
   I = uint8(I*255); 
end

end

