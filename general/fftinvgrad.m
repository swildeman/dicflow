function f = fftinvgrad( fx, fy, varargin )
% FFTINVGRAD Integrate a gradient field (up to a constant) using FFT2
% 
% SYNOPSIS: f = fftinvgrad( fx, fy, (opt) 'bcFix', ..., (opt) 'gradType', ... )
%
% INPUT fx: gradient accross columns (left to right)
%       fy: gradient accross rows (up to down)
%       'bcFix', {'impulse', 'mirror', 'none'}: (default: impulse) 
%              add a special impulse function on the 
%              boundaries of the domain to correct for periodicity assumed 
%              in FFT
%       'gradType', {'spectral', 'difference'}: (default: 'spectral')
%              when 'difference' is specified, [fx,fy] is assumed to come 
%              from a central difference gradient, instead of being analytical
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

% -- Parse input options -- %
defaultBcFix = 'impulse';
expectBcFix = {'impulse', 'mirror', 'none'};

defaultGradType = 'spectral';
expectGradType = {'spectral', 'difference'};

p = inputParser;
p.addParameter('bcfix', defaultBcFix,...
    @(x) any(validatestring(x,expectBcFix)));
p.addParameter('gradtype', defaultGradType, ...
    @(x) any(validatestring(x,expectGradType)));
p.parse(varargin{:});

impbc = false;
mirbc = false;
switch validatestring(p.Results.bcfix, expectBcFix)
    case 'impulse'
        impbc = true;
    case 'mirror'
        mirbc = true;
end

finitediff = false;
switch validatestring(p.Results.gradtype, expectGradType)
    case 'difference'
        finitediff = true;
end

if mirbc
   fx = [fx, -fliplr(fx);  flipud(fx), -rot90(fx,2)];
   fy = [fy,  fliplr(fy); -flipud(fy), -rot90(fy,2)];
end

% -- Integration algorithm -- %
[Ny, Nx] = size(fx);

if impbc
    % add impulse to boundaries to compensate for non-periodicity
    impIx = -.5*sum(fx(:,2:end-1),2);
    impIy = -.5*sum(fy(2:end-1,:),1);
    fx_edge = fx(:,[1,end]); % save for later recovery
    fy_edge = fy([1,end],:);
    fx(:,[1,end]) = [impIx, impIx];
    fy([1,end],:) = [impIy; impIy];
end

% the fourier method will implicitly subtract mean from gradient to satisfy 
% the periodicity assumption, we will tag it back on later
mx = mean(fx(:));
my = mean(fy(:));

% generate k vectors according to matlab's fft layout
if finitediff
    [kx, ky] = meshgrid( sin(kvec(Nx)), sin(kvec(Ny)) );
else
    [kx, ky] = meshgrid( kvec(Nx), kvec(Ny));
end
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
k2(k2 < eps) = 1; % shortcut to prevent division by zero (this effectively subtracts a linear plane)
f_hat = (-1i*kx.*fx_hat + -1i*ky.*fy_hat)./k2;

% transform back to spatial domain
f = real(ifft2(f_hat));

% add mean slope back on
[x,y] = meshgrid(0:Nx-1,0:Ny-1);
f = f + mx*x + my*y;

if impbc % fix edges
    % use inverted second order difference at edges
    f(:,1) = (4*f(:,2) - f(:,3) - 2*fx_edge(:,1)) / 3;
    f(:,end) = (4*f(:,end-1) - f(:,end-2) + 2*fx_edge(:,2)) / 3;
    f(1,:) = (4*f(2,:) - f(3,:) - 2*fy_edge(1,:)) / 3;
    f(end,:) = (4*f(end-1,:) - f(end-2,:) + 2*fy_edge(2,:)) / 3;
elseif mirbc
    f = f(1:Ny/2,1:Nx/2);
end

end

