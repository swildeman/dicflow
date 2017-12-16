rng(21312)
[ Iref, Idef, hgt, ugt, vgt ] = df_testimages();

% make uint8 to fake real image data
% Iref = uint8(Iref*255);
% Idef = uint8(Idef*255);

blockspacing = 4;   % resolution
borderoverlap = 4;  % blocksize = blockspacing + 2*borderoverlap
maxdisp = 4;        % sets size of search window = blocksize + 2*maxdisp
minquality = 0.01;  % minimum quality of the correlation to be a valid vector, NaN otherwise
normcorr = true;    % use normalized cross correlation? relatively slow but robust

tic
[u, v, br, bc, q] = dic_dispfield(Iref, Idef, blockspacing, borderoverlap, maxdisp, [], minquality, normcorr);

u = medfilt2(u,'symmetric');
v = medfilt2(v,'symmetric');

% warp cycle(s) if required
Iref_w = interpimwarp(Iref, u, v, bc, br, 'cubic');
smalldisp = 2; % allow small displacement around current solution
[du, dv, ~, ~, q] = dic_dispfield(Iref_w, Idef, blockspacing, borderoverlap, maxdisp, smalldisp, minquality, normcorr);
u = u + du;
v = v + dv;

u = medfilt2(u,'symmetric');
v = medfilt2(v,'symmetric');

% final interpolation for display
[Xq, Yq] = meshgrid(1:size(Iref,2), 1:size(Iref,1));
% u = interp2(bc,br,u,Xq,Yq,'linear',0);
% v = interp2(bc,br,v,Xq,Yq,'linear',0);
[Iref_w, u, v] = interpimwarp(Iref, u, v, bc, br, 'cubic');
roi = ~isnan(u);
[du, dv] = of_dispfield(Iref_w,Idef, .1, roi);
u = u + du;
v = v + dv;

toc

figure(1)
imshow(Idef.*.3);
hold on
quiver(Xq, Yq, u, v, 0, 'y')
axis image
axis ij
hold off

figure(2)
imshowpair(ugt, u,'montage','scaling','joint')

figure(3)
imshowpair(vgt, v,'montage','scaling','joint')

figure(4)
imagesc(q)
axis image

figure(5)
plot(vgt(:,150))
hold on
plot(v(:,150));
hold off