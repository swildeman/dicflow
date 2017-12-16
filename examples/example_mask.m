rng(234235);
Iref = particleimage(513,513,3,1,1);
Idef = Iref(1:end-1,1:end-1); % shift reference image one pixel right and one pixel down
Iref = Iref(2:end,2:end);

blockspacing = 4;   % resolution
borderoverlap = 4;  % blocksize = blockspacing + 2*borderoverlap
maxdisp = 4;        % sets size of search window = blocksize + 2*maxdisp
minquality = 0.02;  % minimum quality of the correlation to be a valid vector, NaN otherwise
normcorr = true;    % use normalized cross correlation? relatively slow but robust

wgt = true(size(Iref));
[X,Y] = meshgrid(1:512);
wgt((X-150).^2+(Y-200).^2 < 100^2) = false;

figure(1)
imshow(wgt)

Iref = Iref.*wgt;
Idef = Idef.*wgt;

[u, v, br, bc, q] = dic_dispfield(Iref, Idef, blockspacing, borderoverlap, maxdisp,[], minquality, normcorr, wgt);

[Iref_w, u, v] = interpimwarp(Iref, u, v, bc, br, 'cubic');
roi = ~isnan(u);
[du, dv] = of_dispfield(Iref_w,Idef, .1, roi);
u = u + du;
v = v + dv;

figure(2)
imagesc(v)
colormap(parula(256))
axis image