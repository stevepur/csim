%%
respRe = fitsread('resp_m20_20_1_3/testResponse/responseData_s0_p0_re.fits');
respIm = fitsread('resp_m20_20_1_3/testResponse/responseData_s0_p0_im.fits');
resp = respRe + 1i*respIm;

fpmSags = fitsread('wfirst_sags.fits');
fpmSags(end+1) = 0;

maskAmp = fitsread('resp_m20_20_1_3/useOnlyThisHexMask_-1_am.fits');
maskPh = fitsread('resp_m20_20_1_3/useOnlyThisHexMask_-1_ph.fits');

lambda = 5.62175e-07; 

maskAmpReconstructed = zeros(size(maskAmp));
mask = complex(maskAmpReconstructed, maskAmpReconstructed);
for i=1:length(fpmSags)
    fprintf(1,'%d\n', i);
    fn = ['resp_m20_20_1_3/useOnlyThisHexMask_' num2str(i-1) '_am.fits'];
    maskComp = fitsread(fn);
    maskAmpReconstructed = maskAmpReconstructed + maskComp;
    
    mask = mask + maskComp*exp(-4*1i*pi*fpmSags(i)/lambda);
end

figure();
subplot(1,3,1);
imagesc(maskAmp);
title('original mask amplitude');
colorbar;
subplot(1,3,2);
imagesc(abs(mask));
colorbar;
title('reconstructed mask amplitude');
subplot(1,3,3);
imagesc(maskAmp - abs(mask));
colorbar;
title('maskAmp - reconstructed mask amplitude');

figure();
subplot(1,3,1);
imagesc(maskPh);
title('original mask phase');
colorbar;
subplot(1,3,2);
imagesc(angle(mask));
colorbar;
title('reconstructed mask phase');
subplot(1,3,3);
imagesc(maskPh - angle(mask));
colorbar;
title('maskPh - reconstructed mask phase');

%%

% dataDir = 'resp_m20_20_1_3';
% sagFile = 'wfirst_sags.fits';

dataDir = 'resp_m180_180_1_3';
sagFile = 'wfirst_sags.fits';

% dataDir = 'resp_full_0sag';
% sagFile = 'zeroSag.fits';

% dataDir = 'resp_hex0_sag1_all';
% sagFile = 'oneSag.fits';

respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);
resp = respRe + 1i*respIm;

eRe = fitsread([dataDir '/trueEVecRe.fits'])';
eIm = fitsread([dataDir '/trueEVecIm.fits'])';
e = eRe + 1i*eIm;

fpmSags = fitsread(sagFile)';
fpmSags(end+1) = 0;

lambda = 5.62175e-07; 

sv = exp(-4*1i*pi*fpmSags/lambda);

pixelE = resp*sv;

xi = 1:length(pixelE);

respE = sum(resp,2);

figure;
subplot(1,2,1);
plot(xi, real(pixelE), 'd-', xi, real(e), '+-', xi, real(respE), 'o-');
title('real(E)');
legend('reconstructed E', 'original E', 'zero-sag E');
subplot(1,2,2);
plot(xi, imag(pixelE), 'd-', xi, imag(e), '+-', xi, imag(respE), 'o-');
title('imag(E)');

figure;
subplot(1,2,1);
plot(xi, abs(pixelE), 'd-', xi, abs(e), '+-', xi, abs(respE), 'o-');
title('abs(E)');
legend('reconstructed E', 'original E', 'zero-sag E');
subplot(1,2,2);
plot(xi, angle(pixelE), 'd-', xi, angle(e), '+-', xi, angle(respE), 'o-');
title('angle(E)');

figure;
subplot(1,2,1);
plot(xi, abs(pixelE) - abs(e), '+-');
title('abs(reconstructed E) - abs(original E)');
subplot(1,2,2);
plot(xi, angle(pixelE) - angle(e), '+-');
title('angle(reconstructed E) - angle(original E)');

figure;
subplot(1,2,1);
plot(xi, real(pixelE) - real(e), '+-');
title('real(reconstructed E) - real(original E)');
subplot(1,2,2);
plot(xi, imag(pixelE) - imag(e), '+-');
title('imag(reconstructed E) - imag(original E)');


calibMaxIntensity = 1.71144e+07;
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;

regionContrast = mean(pixelN);

disp(['constrast for this region: ' num2str(regionContrast)]);

%% test 5 random sag values

sagMat = [];
trueE = [];
for i=1:5
    dataDir = ['randSag' num2str(i)];
    sagFile = [dataDir '/randSag' num2str(i) '.fits'];
    
    % add an extra column with zero sag for the component outside the mask
    sagMat = [sagMat; [fitsread(sagFile), 0] ];
    
    eRe = fitsread([dataDir '/trueEVecRe.fits'])';
    eIm = fitsread([dataDir '/trueEVecIm.fits'])';
    e = eRe + 1i*eIm;
    trueE = [trueE, e];
end

respRe = fitsread('resp_m20_20_1_3/testResponse/responseData_s0_p0_re.fits');
respIm = fitsread('resp_m20_20_1_3/testResponse/responseData_s0_p0_im.fits');
resp = respRe + 1i*respIm;

lambda = 5.62175e-07; 

sv = exp(4*1i*pi*sagMat/lambda)';

pixelE = resp*sv;

xi = 1:size(pixelE, 1);

figure;
subplot(1,2,1);
plot(xi, real(pixelE), '+-', xi, real(trueE), 'o-');
title('real(E)');
legend('reconstructed E', 'original E');
subplot(1,2,2);
plot(xi, imag(pixelE), '+-', xi, imag(trueE), 'o-');
title('angle(E)');

figure;
subplot(1,2,1);
plot(xi, abs(pixelE), '+-', xi, abs(trueE), 'o-');
title('abs(E)');
legend('reconstructed E', 'original E');
subplot(1,2,2);
plot(xi, angle(pixelE), '+-', xi, angle(trueE), 'o-');
title('angle(E)');

figure;
subplot(1,2,1);
plot(xi, abs(pixelE) - abs(trueE), '+-');
title('abs(reconstructed E) - abs(original E)');
subplot(1,2,2);
plot(xi, angle(pixelE) - angle(trueE), '+-');
title('angle(reconstructed E) - angle(original E)');

figure;
subplot(1,2,1);
plot(xi, real(pixelE) - real(trueE), '+-');
title('real(reconstructed E) - real(original E)');
subplot(1,2,2);
plot(xi, imag(pixelE) - imag(trueE), '+-');
title('imag(reconstructed E) - imag(original E)');


%% play with varying a couple sags to get a sense of the landscape

% dataDir = 'resp_m20_20_1_3';
% sagFile = 'wfirst_sags.fits';

dataDir = 'resp_m180_180_1_3';
sagFile = 'sag_m180_180_1_3.fits';

% dataDir = 'resp_m180_180_1_3';
% sagFile = 'wfirst_sags.fits';

respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);
resp = respRe + 1i*respIm;

fpmSags = fitsread(sagFile)';
fpmSags(end+1) = 0;

lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

% compute baseline contrast
sv = exp(-4*1i*pi*fpmSags/lambda);
pixelE = resp*sv;
respE = sum(resp,2);
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
baseContrast = mean(pixelN);

sagIdx = [100, 508, 526, 786];
nSagSteps = 1000;
sagStep = (max(fpmSags) - min(fpmSags))/nSagSteps;
sagRange = min(fpmSags):sagStep:max(fpmSags);

regionContrast = zeros(length(sagIdx), length(sagRange));
for si=1:length(sagIdx)
    sagMat = repmat(fpmSags, 1, length(sagRange));
    sagMat(sagIdx(si), :) = sagRange';

    sv = exp(-4*1i*pi*sagMat/lambda);
    pixelE = resp*sv;

    pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
    regionContrast(si, :) = mean(pixelN);
end
figure;
semilogy(regionContrast');
hold on;
line([0 length(sagRange)], [baseContrast baseContrast]);
hold off;

%% play with varying all sags to get a sense of the landscape

% dataDir = 'resp_m20_20_1_3';
% sagFile = 'wfirst_sags.fits';

% dataDir = 'resp_m20_20_1_3';
% sagFile = 'sag_m20_20_1_3.fits';

dataDir = 'resp_m180_180_1_3';
sagFile = 'sag_m180_180_1_3.fits';
% 
% dataDir = 'resp_m180_180_1_3';
% sagFile = 'wfirst_sags.fits';

respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);
resp = respRe + 1i*respIm;

fpmSags = fitsread(sagFile)';
fpmSags(end+1) = 0;

lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

% compute baseline contrast
sv = exp(-4*1i*pi*fpmSags/lambda);
pixelE = resp*sv;
respE = sum(resp,2);
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
baseContrast = mean(pixelN);

nSagSteps = 1000;
sagStep = (max(fpmSags) - min(fpmSags))/nSagSteps;
sagRange = min(fpmSags):sagStep:max(fpmSags);

% regionContrast = zeros(length(fpmSags)-1, length(sagRange));
regionContrast = [];


% partition fpmSags
nParts = 12;
partSize = floor((length(fpmSags)-1)/nParts);
for p = 1:nParts
    sagStart = (p-1)*partSize + 1;
    sagEnd = sagStart + partSize;
    disp(['sags from ' num2str(sagStart) ' to ' num2str(sagEnd)]);
    sagIdx = sagStart:sagEnd;
    for si=1:length(sagIdx)
        sagMat = repmat(fpmSags, 1, length(sagRange));
%         sagMat(sagIdx(si), :) = fpmSags(sagIdx(si)) + sagRange';
        sagMat(sagIdx(si), :) = sagRange';

        sv = exp(-4*1i*pi*sagMat/lambda);
        pixelE = resp*sv;

        pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
%         regionContrast(sagIdx(si), :) = mean(pixelN);
        regionContrast = [regionContrast; mean(pixelN)];
    end
end

figure;
semilogy(sagRange, regionContrast');
hold on;
line([sagRange(1), sagRange(end)], [baseContrast baseContrast], 'Color', 'k', 'LineWidth', 1);
hold off;
title(dataDir, 'interpreter', 'none');
axis tight;
%% play with varying two sags to get a sense of the landscape

% dataDir = 'resp_m20_20_1_3';
% sagFile = 'wfirst_sags.fits';

dataDir = 'resp_m180_180_1_3';
sagFile = 'sag_m180_180_1_3.fits';

% dataDir = 'resp_m180_180_1_3';
% sagFile = 'wfirst_sags.fits';

respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);
resp = respRe + 1i*respIm;

fpmSags = fitsread(sagFile)';
fpmSags(end+1) = 0;

lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

% compute baseline contrast
sv = exp(-4*1i*pi*fpmSags/lambda);
pixelE = resp*sv;
respE = sum(resp,2);
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
baseContrast = mean(pixelN);

nSagSteps = 100;
sagStep = (max(fpmSags) - min(fpmSags))/nSagSteps;
sagRange = min(fpmSags):sagStep:max(fpmSags);

% regionContrast = zeros(length(fpmSags)-1, length(sagRange));
regionContrast = [];

for t=1:5
    sagIdx = fix(rand(1,2)*length(fpmSags)-2) + 1;
    for s1 = 1:length(sagRange)
        for s2 = 1:length(sagRange)
            sagMat = fpmSags;
            sagMat(sagIdx(1)) = sagRange(s1);
            sagMat(sagIdx(2)) = sagRange(s2);

            sv = exp(-4*1i*pi*sagMat/lambda);
            pixelE = resp*sv;

            pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
            regionContrast(s1, s2) = mean(pixelN);
        end
    end
    figure;
    imagesc(sagRange, sagRange, log10(regionContrast));
    hold on;
    plot(fpmSags(sagIdx(2)), fpmSags(sagIdx(1)), 'r*');
    hold off;
    title([dataDir ' sag indices ' num2str(sagIdx)], 'interpreter', 'none');
    xlabel('sag values');
    ylabel('sag values');
    colorbar;
end
%% play with varying three sags to get a sense of the landscape

dataDir = 'resp_m20_20_1_3';
sagFile = 'wfirst_sags.fits';

respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);
resp = respRe + 1i*respIm;

fpmSags = fitsread(sagFile)';
fpmSags(end+1) = 0;

lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

% compute baseline contrast
sv = exp(-4*1i*pi*fpmSags/lambda);
pixelE = resp*sv;
respE = sum(resp,2);
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
baseContrast = mean(pixelN);

nSagSteps = 100;
sagStep = (max(fpmSags) - min(fpmSags))/nSagSteps;
sagRange = min(fpmSags):sagStep:max(fpmSags);

regionContrast = zeros([length(sagRange), length(sagRange), length(sagRange)]);
sagIdx = fix(rand(1,3)*length(fpmSags)-2) + 1;
for s1 = 1:length(sagRange)
    disp(s1);
    for s2 = 1:length(sagRange)
        for s3 = 1:length(sagRange)
                sagMat = fpmSags;
                sagMat(sagIdx(1)) = sagRange(s1);
                sagMat(sagIdx(2)) = sagRange(s2);
                sagMat(sagIdx(3)) = sagRange(s3);

                sv = exp(-4*1i*pi*sagMat/lambda);
                pixelE = resp*sv;

                pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;
                regionContrast(s1, s2, s3) = mean(pixelN);
        end
    end
end
%%
figure;
imagesc(log10(squeeze(regionContrast(:,:,45))));
title(['sag indices ' num2str(sagIdx)]);
colorbar;
