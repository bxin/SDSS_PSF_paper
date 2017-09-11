M = importdata('../SDSSdata/masterTXT/run3427.txt', ' ', 1);

%#field 	 camCol  filter FWHMvK psf_width airmass mjd		 psf_nstar  neff_psf	 sky_frames
idxu=(M.data(:,3)==0);

subplot(2,1,1);
xu = M.data(idxu,6);

yu = M.data(idxu,4);
nbin = 8;
xblim = linspace(min(xu)-0.01,max(xu)+0.01,nbin+1); %bin limits
xb = (xblim(1:nbin)+xblim(2:end))/2;
yb = zeros(1,nbin);
eb = zeros(1,nbin);
for ib=1:nbin
   idx = (xu>=xblim(ib) & xu<=xblim(ib+1) );
   yb(ib) = mean(yu(idx));
   eb(ib) = std(yu(idx));
end
plot(xu, yu,'ro');
hold on;
errorbar(xb, yb, eb, 'b.', 'markersize',50,'linewidth',2);
