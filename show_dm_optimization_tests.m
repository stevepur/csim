nelderMeadeBest = textread('dmOptimization/optBestVal_n10_nelderMead.txt');
subPlexBest = textread('dmOptimization/optBestVal_n10_subPlex.txt');
cobylaBest = textread('dmOptimization/optBestVal_n10_cobyla.txt');
bobyqaBest = textread('dmOptimization/optBestVal_n10_bobyqa.txt');
newuoaBest = textread('dmOptimization/optBestVal_n10_newuoa.txt');
figure;
semilogy(nelderMeadeBest(:,1)/60, nelderMeadeBest(:,2), ...
    subPlexBest(:,1)/60, subPlexBest(:,2), ...
    cobylaBest(:,1)/60, cobylaBest(:,2), ...
    bobyqaBest(:,1)/60, bobyqaBest(:,2), ...
    newuoaBest(:,1)/60, newuoaBest(:,2));
title('NLOpt Best Value, N=10');
xlabel('time (minutes)');
ylabel('contrast');
legend('Nelder-Mead', 'subPlex', 'cobyla', 'bobyqa', 'newuoa');

nelderMeadeHistory = textread('dmOptimization/optValHistory_n10_nelderMead.txt');
subPlexHistory = textread('dmOptimization/optValHistory_n10_subPlex.txt');
cobylaHistory = textread('dmOptimization/optValHistory_n10_cobyla.txt');
bobyqaHistory = textread('dmOptimization/optValHistory_n10_bobyqa.txt');
newuoaHistory = textread('dmOptimization/optValHistory_n10_newuoa.txt');
figure;
semilogy(nelderMeadeHistory(:,1)/60, nelderMeadeHistory(:,2), nelderMeadeBest(:,1)/60, nelderMeadeBest(:,2), ...
    subPlexHistory(:,1)/60, subPlexHistory(:,2), subPlexBest(:,1)/60, subPlexBest(:,2), ...
    cobylaHistory(:,1)/60, cobylaHistory(:,2), cobylaBest(:,1)/60, cobylaBest(:,2), ...
    bobyqaHistory(:,1)/60, bobyqaHistory(:,2), bobyqaBest(:,1)/60, bobyqaBest(:,2), ...
    newuoaHistory(:,1)/60, newuoaHistory(:,2), newuoaBest(:,1)/60, newuoaBest(:,2));
title('NLOpt History, N=10');
xlabel('time (minutes)');
ylabel('contrast');
legend('Nelder-Mead history', 'Nelder-Mead best', 'subPlex history', 'subPlex best', ...
    'cobyla history', 'cobyla best', 'bobyqa history', 'bobyqa best', 'newuoa history', 'newuoa best');
grid on;

semilogy(1:size(nelderMeadeHistory, 1), nelderMeadeHistory(:,2), ...
    1:size(subPlexHistory, 1), subPlexHistory(:,2), ...
    1:size(cobylaHistory, 1), cobylaHistory(:,2), ...
    1:size(bobyqaHistory, 1), bobyqaHistory(:,2), ...
    1:size(newuoaHistory, 1), newuoaHistory(:,2));
title('NLOpt Best Value, N=10');
xlabel('evaluation');
ylabel('contrast');
legend('Nelder-Mead', 'subPlex', 'cobyla', 'bobyqa', 'newuoa');
grid on;

%%
bobyqaBestn100 = textread('dmOptimization/optBestVal_n100_bobyqa.txt');
bobyqaHistoryn100 = textread('dmOptimization/optValHistory_n100_bobyqa.txt');
figure;
subplot(1,2,1);
semilogy(bobyqaBestn100(:,1)/60/60, bobyqaBestn100(:,2));
title(['NLOpt BOBYQA Best Value, N=100: ' num2str(min(bobyqaBestn100(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
grid on;
subplot(1,2,2);
semilogy(1:size(bobyqaHistoryn100, 1), bobyqaHistoryn100(:,2));
title(['NLOpt BOBYQA Hisitory,  N=100: ' num2str(min(bobyqaHistoryn100(:,2)))]);
xlabel('evaluation');
ylabel('contrast');
grid on;

%%
currentBest = textread('dmOptimization/optBestVal.txt');
currentHistory = textread('dmOptimization/optValHistory.txt');
figure;
subplot(1,3,1);
semilogy(currentBest(:,1)/60/60, currentBest(:,2));
title(['NLOpt Best Value, current run: ' num2str(min(currentBest(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
grid on;
subplot(1,3,2);
semilogy(currentHistory(:,1)/60/60, currentHistory(:,2), currentBest(:,1)/60/60, currentBest(:,2));
title(['NLOpt History, current run: ' num2str(min(currentBest(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
grid on;
subplot(1,3,3);
loglog(1:size(currentHistory, 1), currentHistory(:,2));
title(['NLOpt History, current run: ' num2str(min(currentHistory(:,2)))]);
xlabel('evaluation');
ylabel('contrast');
grid on;

%%
% optDir = 'dmOptimization/mono_650nm_mma_full_R17_fromRandom/';
% optDir = 'dmOptimization/poly_3w_bobyqa_opt_6/';
optDir = 'dmOptimization/linearOptLargeStar/';
% optDir = 'dmOptTest/';
currentBest = textread([optDir 'optBestVal.txt']);
currentHistory = textread([optDir 'optValHistory.txt']);

figure;
subplot(1,3,1);
semilogy(currentHistory(:,1)/60/60, currentHistory(:,2), '+-', currentBest(:,1)/60/60, currentBest(:,2), 'o-');
title(['NLOpt Broadband History, current run: ' num2str(min(currentBest(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
legend('all evaluations', 'best value');
grid on;
subplot(1,3,2);
loglog(1:size(currentHistory, 1), currentHistory(:,2), '+-');
title(['NLOpt Broadband History, current run: ' num2str(min(currentHistory(:,2)))]);
xlabel('evaluation');
ylabel('contrast');
grid on;
subplot(1,3,3);
plot(1:size(currentHistory, 1), currentHistory(:,2), '+-');
title(['NLOpt Broadband History, current run: ' num2str(min(currentHistory(:,2)))]);
xlabel('evaluation');
ylabel('contrast');
grid on;

%%
figure;
dmAct = fitsread([optDir2 'optimalActuators.fits']);
imagesc(dmAct);
title(['NLOpt Broadband Actuators, current run: ' num2str(min(currentBest2(:,2)))]);
axis equal;
axis tight;
colorbar;

%% linear plot

currentBest = textread('dmOptimization/poly_3w_bobyqa_opt_6/optBestVal.txt');
currentHistory = textread('dmOptimization/poly_3w_bobyqa_opt_6/optValHistory.txt');
figure;
t0 = 7.9;
tcut = 33.9;
tcut2 = 47;
t = currentHistory(:,1)/60/60;
tt1 = find(t>t0 & t < tcut);
ch1 = currentHistory(tt1,2);
tt2 = find(t>tcut & t<tcut2);
ch2 = currentHistory(tt2,2);
tt3 = find(t>tcut2);
ch3 = currentHistory(tt3,2);
b = currentBest(:,1)/60/60;
bb1 = find(b>t0 & b<tcut);
cb1 = currentBest(bb1,2);
bb2 = find(b>tcut & b<tcut2);
cb2 = currentBest(bb2,2);
bb3 = find(b>tcut2);
cb3 = currentBest(bb3,2);
plot(t(tt1), ch1, b(bb1), cb1, t(tt2), ch2, b(bb2), cb2, t(tt3), ch3, b(bb3), cb3);
title(['NLOpt Broadband History, current run: ' num2str(min(currentBest(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
legend('all evaluations', 'best value', 'all evaluations', 'best value', 'all evaluations', 'best value');
grid on;

%%
optDir1 = 'dmOptimization/linearOpt/';
optDir2 = 'dmOptimization/linearOptTest1/';
% optDir = 'dmOptTest/';
currentBest1 = textread([optDir1 'optBestVal.txt']);
currentBest2 = textread([optDir2 'optBestVal.txt']);

figure;
subplot(1,2,1);
semilogy(currentBest1(:,1)/60/60, currentBest1(:,2), '+-', currentBest2(:,1)/60/60, currentBest2(:,2), 'o-');
title(['NLOpt Broadband History, current run: ' num2str(min(currentBest2(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
legend('\mu = 5e7', '\mu = 3e7');
grid on;
subplot(1,2,2);
plot(currentBest1(:,1)/60/60, currentBest1(:,2), '+-', currentBest2(:,1)/60/60, currentBest2(:,2), 'o-');
title(['NLOpt Broadband History, current run: ' num2str(min(currentBest2(:,2)))]);
xlabel('time (hours)');
ylabel('contrast');
grid on;




