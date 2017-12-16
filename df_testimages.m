function [ Iref, Idef, h, u, v ] = df_testimages()
%TESTIMAGES Generate pair of test images for DIC/OF analysis
%
% SYNOPSIS: [ Iref, Idef, h, u, v ] = df_testimages()
%
% OUTPUT Iref: Reference dot-pattern
%        Idef: Iref distorted by wave-like disturbance
%        h: Ground-truth height field so that [u,v] = grad(-h);
%        u,v: Ground-truth displacement field
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

brd = 4;
sz = 512;

[x,y] = meshgrid(-brd:sz-1+brd,-brd:sz-1+brd);

xp = x-150;
yp = y-200;

kwave = 2*pi/(32);
r = sqrt(xp.^2 + yp.^2);
a = 20;
h = a*besselj(0,kwave*r);
% [u,v] = gradient(-h);
u = a * kwave * xp .* besselj(1,kwave*r) ./ r;
v = a * kwave * yp .* besselj(1,kwave*r) ./ r;

u(isnan(u)) = 0;
v(isnan(v)) = 0;

Iref = particleimage(sz+2*brd,sz+2*brd,3,1,1);

D = zeros([size(u),2]);
D(:,:,1) = -u;
D(:,:,2) = -v;
Idef = imwarp(Iref,D,'cubic');

h = h(1+brd:end-brd, 1+brd:end-brd);
u = u(1+brd:end-brd, 1+brd:end-brd);
v = v(1+brd:end-brd, 1+brd:end-brd);
Idef = Idef(1+brd:end-brd, 1+brd:end-brd);
Iref = Iref(1+brd:end-brd, 1+brd:end-brd);

end

