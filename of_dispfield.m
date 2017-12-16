function [u,v] = of_dispfield( I1, I2, alpha, roi, Dx, Dy, Lap )
% OF_DISPFIELD Use Horn-Schunck method to estimate the optical flow between
% two images.
%
% SYNOPSIS: [u,v] = of_dispfield( I1, I2, alpha, roi, Dx, Dy, Lap )
%
% INPUT I1,I2: calculate optical flow from image I1 to image I2
%       alpha: regulating parameter HS method
%       (optional) Dx,Dy: precomputed matrix such that Dn*I(:) gives 
%                         n-derivative of image (in unrolled format)
%       (optional) Lap: precomputed matrix such that Lap*I gives laplacian 
%                       of image (in unrolled format)
%
% OUTPUT u,v: displacements in x (column) and y (row) directions
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file

I1 = double(I1);
I2 = double(I2);

if isa(I1, 'gpuArray')
    ongpu = true;
else
    ongpu = false;
end

if nargin < 4
    roi = ~isnan(I1) & ~isnan(I2);
else
    roi = logical(roi);
end

if nargin < 5
    [Dx,Dy] = designgrad(roi);
    [Dxx,Dyy] = designlap(roi,'bc', 'neumann');
    Lap = Dxx + Dyy;
end

N = sum(roi(:));

% image derivatives
Ix = sparse(1:N,1:N, Dx*I1(roi), N,N);
Iy = sparse(1:N,1:N, Dy*I1(roi), N,N);
It = I2(roi)-I1(roi);

% calculate matrix elements
A11 = Ix*Ix-alpha^2*Lap;
A12_21 = Ix*Iy;
A22 = Iy*Iy-alpha^2*Lap;

% gather stuff from the gpu, because matlab does not support horzcat on
% sparse gpuArrays yet...
if ongpu
    A11 = gather(A11);
    A12_21 = gather(A12_21);
    A22 = gather(A22);
end

% Horn-Schunck system in A*x = b form
A = [A11, A12_21;
     A12_21, A22];

b = full([Ix*It;Iy*It]);

if ongpu
   A = gpuArray(A);
   b = gpuArray(b);
end

% Solve using conjugate gradient method
D = cgs(A,b,1e-2,1000);

% Return displacement vectors
D = reshape(D, [N,2]);

u = NaN(size(I1),'like',I1);
v = NaN(size(I1),'like',I1);

u(roi) = -D(:,1);
v(roi) = -D(:,2);

end