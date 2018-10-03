
calibration_PSF_nominal = fitsread('normalized_PSF_nominal.fits');
calibration_PSF_nominal_wp = fitsread('normalized_PSF_nominal_wp.fits');
calibration_PSF_opt0 = fitsread('normalized_PSF_opt0.fits');
calibration_PSF_opt0_wp = fitsread('normalized_PSF_opt0_wp.fits');

figure;
subplot(2,2,1);
imagesc(log10(calibration_PSF_nominal));
axis equal;
axis tight;
caxis([-10 -3]);
title('Dan''s optimization, no planet');
colorbar;
subplot(2,2,2);
imagesc(log10(calibration_PSF_nominal_wp));
axis equal;
axis tight;
caxis([-10 -3]);
title('Dan''s optimization, with planet');
colorbar;
subplot(2,2,3);
imagesc(log10(calibration_PSF_opt0));
axis equal;
axis tight;
caxis([-10 -3]);
title('Steve''s optimization, no planet');
colorbar;
subplot(2,2,4);
imagesc(log10(calibration_PSF_opt0_wp));
axis equal;
axis tight;
caxis([-10 -3]);
title('Steve''s optimization, with planet');
colorbar