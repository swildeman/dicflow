function f = fftinvgrad( fx, fy, impbc )
% FFTINVGRAD Integrate a gradient field (up to a constant) using FFT2
% 
% SYNOPSIS: f = fftinvgrad( fx, fy, impbc )
%
% INPUT fx: gradient accross columns (left to right)
%       fy: gradient accross rows (up to down)
%       impbc: (default: true) add a special impulse function on the 
%              boundaries of the domain to correct for periodicity assumed 
%              in FFT 
%
% OUTPUT f: integrated gradient field such that [fx,fy] = grad(f) in a
%           least squares sense
%
% REMARKS In case large values are present near the boundaries of the 
% integrated field or when some regions in the field are masked, the slower 
% but more flexible/robust <a href="matlab:help invgrad2">invgrad2</a> function is preferred.
%
% Based on <a href="https://doi.org/10.1007/s00348-016-2236-3">Huhn, et al. Exp Fluids (2016), 57, 151</a>
%
% See also:
% INVGRAD2
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file


if nargin < 3
   impbc = true; 
end

[Ny, Nx] = size(fx);

if impbc
    % add impulse to boundaries to compensate for non-periodicity
    impIx = -.5*sum(fx(:,2:end-1),2);
    impIy = -.5*sum(fy(2:end-1,:),1);
    fx(:,[1,end]) = [impIx, impIx];
    fy([1,end],:) = [impIy; impIy];
end

% the fourier method will implicitly subtract mean from gradient to satisfy 
% the periodicity assumption, we will tag it back on later
mx = mean(fx(:));
my = mean(fy(:));

% generate k vectors according to matlab's fft layout
[kx, ky] = meshgrid(kvec(Nx),kvec(Ny));
% pre-compute k^2;
k2 = kx.^2 + ky.^2;

if mod(Nx,2)==0 % if Nx even
    kx(:,Nx/2+1) = 0; % remove degeneracy at kx=Nx/2 leading to imaginary part
end

if mod(Ny,2)==0 % if Ny even
    ky(Ny/2+1,:) = 0; % remove degeneracy at ky=Nx/2 leading to imaginary part
end

% compute fft of gradients
fx_hat = fft2(fx);
fy_hat = fft2(fy);

% integrate in fourier domain
k2(1,1) = 1; % shortcut to prevent division by zero (this effectively subtracts a linear plane)
f_hat = (-1i*kx.*fx_hat + -1i*ky.*fy_hat)./k2;

% transform back to spatial domain
f = real(ifft2(f_hat));

% add mean slope back on
[x,y] = meshgrid(0:Nx-1,0:Ny-1);
f = f + mx*x + my*y;

end

