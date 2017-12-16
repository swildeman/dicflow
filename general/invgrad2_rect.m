function [ h ] = invgrad2_rect( hx, hy )
%INVGRAD2_RECT Integrate gradient field (up to constant) in rectangular 
%domain by inverting finite difference matrix in least squares sense
% 
% SYNOPSIS: [ h ] = invgrad2_rect( hx, hy )
%
% INPUT hx: gradient accross columns (left to right)
%       hy: gradient accross rows (up to down)
%
% OUTPUT f: integrated gradient field such that [fx,fy] = grad(f) in a
%           least squares sense
%
% Based on <a href="https://doi.org/10.1109/CVPR.2008.4587414">Harker and O'Leary, IEEE CVPR 2008, pp. 1-7</a>
%
% See also:
% FFTINVGRAD
% DESIGNGRAD1D
%
% Copyright (c) 2017 Sander Wildeman
% Distributed under the MIT License, see LICENSE file 

[m,n] = size(hx);

Dx = designgrad1D(n);
Dy = designgrad1D(m);

% Householder reflection matrix
vn = ones(n,1);
vn(1) = 1+sqrt(n);
vm = ones(m,1);
vm(1) = 1+sqrt(m);
Px = housh(vn);
Py = housh(vm);

Dhx = Dx*Px;
Dhy = Dy*Py;

A = Dhy.'*Dhy;
A = A(2:end,2:end);
B = Dhx.'*Dhx;
B = B(2:end,2:end);

C = Dhy.'*hy*Px + Py.'*hx*Dhx;
c01 = C(1,2:end);
c10 = C(2:end,1);
C = C(2:end,2:end);

w01 = c01/B;
w10 = A\c10;
W11 = sylvester(A,B,C);

W = [0 w01;
    w10 W11];

h = Py*W*Px.';

function [ H ] = housh( v )
    %HOUSH Householder reflection matrix
    v = v(:);
    H = eye(length(v)) - 2*(v*v')/(v'*v);
end

end

