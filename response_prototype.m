respRe = fitsread('resp_nom/testResponse/responseData_s0_p0_re.fits');
respIm = fitsread('resp_nom/testResponse/responseData_s0_p0_im.fits');
resp = respRe + 1i*respIm;

eRe = fitsread('resp_nom/trueEVecRe.fits');
eIm = fitsread('resp_nom/trueEVecIm.fits');
e = eRe + 1i*eIm;

fpmSags = fitsread('wfirst_sags.fits');

lambda = 5.62175e-07; 

sv = exp(-4*1i*pi*fpmSags/lambda)';

pixelE = resp*sv;

xi = 1:length(pixelE);
figure;
subplot(1,2,1);
plot(xi, real(pixelE), 'd-', xi, real(e), '+-');
title('nominal real(E)');
legend('reconstructed E', 'original E');
subplot(1,2,2);
plot(xi, imag(pixelE), 'd-', xi, imag(e), '+-');
title('nominal imag(E)');

figure;
subplot(1,2,1);
plot(xi, abs(pixelE), 'd-', xi, abs(e), '+-');
title('nominal abs(E)');
legend('reconstructed E', 'original E');
subplot(1,2,2);
plot(xi, angle(pixelE), 'd-', xi, angle(e), '+-');
title('nominal angle(E)');



calibMaxIntensity = 1.71144e+07;
pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;

constrast = sum(pixelN);

%%
respRe = fitsread('resp_hex0_sagnom/testResponse/responseData_s0_p0_re.fits');
respIm = fitsread('resp_hex0_sagnom/testResponse/responseData_s0_p0_im.fits');

resp = respRe + 1i*respIm;

eRe = fitsread('resp_hex0_sagnom/trueEVecRe.fits');
eIm = fitsread('resp_hex0_sagnom/trueEVecIm.fits');
e = eRe + 1i*eIm;

fpmSags = fitsread('wfirst_sags.fits');

lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

sv = exp(4*1i*pi*fpmSags/lambda)';
sv(2:end) = 0;

pixelE = resp*sv;

xi = 1:length(pixelE);
figure;
plot(xi, real(pixelE), xi, real(e), xi, real(resp(:,1)));
legend('reconstructed E', 'original E', 'zero-sag E');



pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;

constrast = sum(pixelN);

%% test controlled sag values
sagVal = 1;
dataDir = ['resp_hex0_sag' num2str(sagVal)];
respRe = fitsread([dataDir '/testResponse/responseData_s0_p0_re.fits']);
respIm = fitsread([dataDir '/testResponse/responseData_s0_p0_im.fits']);

resp = respRe + 1i*respIm;

eRe = fitsread([dataDir '/trueEVecRe.fits']);
eIm = fitsread([dataDir '/trueEVecIm.fits']);
e = eRe + 1i*eIm;

% fpmSags = fitsread('wfirst_sags.fits');
% fpmSags = zeros(size(fpmSags));
% fpmSags(1) = sagVal*1e-7;
fpmSags = fitsread('oneSag.fits');
lambda = 5.62175e-07; 
calibMaxIntensity = 1.71144e+07;

sv = exp(-4*1i*pi*fpmSags/lambda)';

pixelE = resp*sv;

xi = 1:length(pixelE);
figure;
subplot(1,2,1);
plot(xi, real(pixelE), 'd-', xi, real(e), '+-', xi, real(resp(:,1)), 'o-');
title(['sagVal = ' num2str(sagVal) ' real(E)']);
legend('reconstructed E', 'original E', 'zero-sag E');
subplot(1,2,2);
plot(xi, imag(pixelE), 'd-', xi, imag(e), '+-', xi, imag(resp(:,1)), 'o-');
title(['sagVal = ' num2str(sagVal) ' imag(E)']);

figure;
subplot(1,2,1);
plot(xi, abs(pixelE), 'd-', xi, abs(e), '+-', xi, abs(resp(:,1)), 'o-');
title(['sagVal = ' num2str(sagVal) ' abs(E)']);
legend('reconstructed E', 'original E', 'zero-sag E');
subplot(1,2,2);
plot(xi, angle(pixelE), 'd-', xi, angle(e), '+-', xi, angle(resp(:,1)), 'o-');
title(['sagVal = ' num2str(sagVal) ' angle(E)']);


pixelN = pixelE.*conj(pixelE)/calibMaxIntensity;

constrast = sum(pixelN);
