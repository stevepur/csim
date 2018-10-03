hexNum = fitsread('segmented/8Rings/luvior_base_mask_data_hexNum.fits');
nSubPix = fitsread('segmented/8Rings/luvior_base_mask_data_nSubPix.fits');
lambda = 1.235e-06;
lambdaRef = 1.3e-6;
NBsubPix = 64;
N = 2048;
beamrad = 0.0188474;
pscale = 0.00003916814414;
fpmFRatio = 85.9;
fpscale = 2*(beamrad/(pscale*N*16))*lambdaRef*fpmFRatio;
sags = fitsread('segmented/8Rings/luvior_olivier_sags.fits'); % Olivier's original sags

sagVals = zeros(size(hexNum));
% set the 3D sag array based on the hex number
for r = 1:size(hexNum, 1)
    for c = 1:size(hexNum, 2)
        for s = 1:size(hexNum, 3)
            if (hexNum(r, c, s) > -1)
                sagVals(r, c, s) = sags(hexNum(r, c, s) + 1); % hexNum is indexed from 0
            else
                sagVals(r, c, s) = 0;
            end
        end
    end
end
mask.M = sum(nSubPix.*(exp(-2*pi*1i*2*sagVals/lambda)), 3)./(NBsubPix*NBsubPix);
%%
figure('Color', 'white');
subplot(2,2,1);
imagesc(abs(mask.M));
title('abs(M)');
axis equal;
subplot(2,2,2);
imagesc(angle(mask.M));
title('angle(M)');
axis equal;
%%
olivierMaskAmp = fitsread('segmented/eval0/fpm_CAampmap2D.s2_l1300_sr10_nbr008_mr200_ssr20_ssm0_wb10.minsag-00500_maxsag000500_fpmregc1000000000_fpmrega001000_Mirror.fits');
olivierMaskPh = fitsread('segmented/eval0/fpm_CAphamap2D.s2_l1300_sr10_nbr008_mr200_ssr20_ssm0_wb10.minsag-00500_maxsag000500_fpmregc1000000000_fpmrega001000_Mirror.fits');
subplot(2,2,3);
imagesc(olivierMaskAmp(:,:,1));
title('abs(M)');
axis equal;
subplot(2,2,4);
imagesc(olivierMaskPh(:,:,1));
title('angle(M)');
axis equal;

