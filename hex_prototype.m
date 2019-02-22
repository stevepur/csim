clear;

fid = fopen('../../../piaa_optimization/PIAACMCdesigns/propoly_test/piaacmcconf_i000/fpmsectors2_032.txt', 'r');
C = textscan(fid, '%d %d %f %f');
fclose(fid);

hexNum = C{1};
hexRing = C{2};
hexX = C{3};
hexY = C{4};
hexInnerRad = (hexX(5) - hexX(1))/2;

% hex coordinates from PIAACMCsimul.c line 2338:
% ii1max = (long) (piaacmc[0].NBrings/3+2) = 32/3 + 2 = 12.6667
% jj1max = (long) (piaacmc[0].NBrings/sqrt(3.0)+2) = 20.4752
% hexX ranges from -ii1max*3 to ii1max*3 or -38 to 38 
% hexy ranges from -jj1max to jj1max -35.4641 to 35.4641
% with each step adding a second point offset by 1.5 in x and sqrt(3)/2 =
% 0.8660 in y

figure('Color', 'white');
plot(hexX, hexY, '+');
axis equal
axis xy
%% recreate the hex positions from scratch
NBrings = 32; % 32 rings
ii1max = fix(NBrings/3+2);
jj1max = fix(NBrings/sqrt(3.0)+2);
hexstep = 1;

% first compute the hex coordinates relative to the ring number
hindex = 1;
badHexNum = 1e8;
for ii1 = -ii1max:ii1max-1
    for jj1 = -jj1max:jj1max-1
        hx = hexstep*ii1*3;
        hy = hexstep*sqrt(3)*jj1;
        ring = fix(sqrt(hx*hx+hy*hy));
        hexXTest(hindex) = hx;
        hexYTest(hindex) = hy;
        hexRingTest(hindex) = ring;
        if ring < NBrings
            hexNumTest(hindex) = 0;
        else
            hexNumTest(hindex) = badHexNum + hindex;
        end
        hindex = hindex + 1;
        
        hx = hx + hexstep*1.5;
        hy = hy + hexstep*sqrt(3)/2;
        ring = fix(sqrt(hx*hx+hy*hy));
        hexXTest(hindex) = hx;
        hexYTest(hindex) = hy;
        hexRingTest(hindex) = ring;
        if ring < NBrings
            hexNumTest(hindex) = 0;
        else
            hexNumTest(hindex) = badHexNum + hindex;
        end
        hindex = hindex + 1;
    end
end
% now organize the list by ring
hCount = 1;
% hexNumTest(hCount) = 0;
% hCount = hCount + 1;
for r=0:NBrings
    rr = find(hexRingTest == r & hexNumTest == 0); % find the hexes on this ring
    for rri = 1:length(rr)
        hexNumTest(rr(rri)) = hCount - 1;
        hCount = hCount + 1;
    end
end

[hexNumTest, hsIdx] = sort(hexNumTest);
hexXTest = hexXTest(hsIdx);
hexYTest = hexYTest(hsIdx); % hexYTest differs from hexY because the hexY text file only has 6 decimal places
hexRingTest = hexRingTest(hsIdx);

hexNumTest = hexNumTest(:);
hexXTest = hexXTest(:);
hexYTest = hexYTest(:);
hexRingTest = hexRingTest(:);

hexXTest = hexXTest(hexNumTest < badHexNum);
hexYTest = hexYTest(hexNumTest < badHexNum);
hexNumTest = hexNumTest(hexNumTest < badHexNum);

% test
if any(hexNum ~= hexNumTest(hexNumTest < badHexNum))
    disp('hexNum does not match hexNumTest!!');
end
if any(hexRing ~= hexRingTest(hexNumTest < badHexNum))
    disp('hexRing does not match hexRingTest!!');
end
dhx = hexX - hexXTest(hexNumTest < badHexNum);
if any(abs(dhx) > 0)
    disp('hexX does not match hexXTest!!');
end
% hexYTest differs from hexY because the hexY text file only has 6 decimal places
dhy = hexY - hexYTest(hexNumTest < badHexNum);
if any(abs(dhy) > 1e-6)
    disp('hexY does not match hexYTest!!');
end

hexNumTest(hexNumTest >= badHexNum) = -1;

% hexNumTest = hexNum;
% hexXTest = hexX;
% hexYTest = hexY;
% hexRingTest = hexRing;
%%
figure('Color', 'white');
text(hexX, hexY, num2str(hexNum));
axis([min(hexX) max(hexX) min(hexY) max(hexY)]);
axis equal
axis xy

figure('Color', 'white');
text(hexXTest, hexYTest, num2str(hexNumTest));
axis([min(hexXTest) max(hexYTest) min(hexXTest) max(hexYTest)]);
axis equal
axis xy

figure('Color', 'white');
text(hexXTest(hexNumTest>-1), hexYTest(hexNumTest>-1), num2str(hexNumTest(hexNumTest>-1)));
axis([min(hexXTest) max(hexYTest) min(hexXTest) max(hexYTest)]);
axis equal
axis xy

%%

fpmPhase = fitsread('tmp_fpmCA_pha.fits');
fpmAmp = fitsread('tmp_fpmCA_ampl.fits');
figure('Color', 'white');
imagesc(fpmPhase(:,:,5));
axis equal
axis xy
%%

sagArray = fitsread('fpm_slice5.fits');
figure('Color', 'white');
imagesc(sagArray);
axis equal
axis xy

%%
fpmSags = fitsread('../../matlab_examples/matlab-wfirstPIAAValidation/proppoly_0c/optNoAstigNoAstig/eval1/piaacmcconf_i000/fpm_zonez_s2_l0565_sr10_nbr032_mr300_minsag-02000_maxsag002000_fpmreg001000_ssr50_ssm0_Mirror_wb10.best.fits');
figure('Color', 'white');
plot(fpmSags);

%%
% cooridnates in the hexagon creation array from PIAACMCsimul.c line 2488:
% FPMSCALEFACTOR = 0.9
% piaacmc[0].fpmarraysize = 2048
% hexsteppix = 0.5*piaacmc[0].fpmarraysize/piaacmc[0].NBrings * FPMSCALEFACTOR 
%               = 0.5*2048/32 * 0.9 = 28.8
% hexagons are made at 0.5*piaacmc[0].fpmarraysize + hexX * hexsteppix = , 
%       0.5*piaacmc[0].fpmarraysize + hexY * hexsteppix
% 
FPMSCALEFACTOR = 0.9; % the size of the mask relative to the array
fpmarraysize = 2048;
hexsteppix = 0.5*fpmarraysize/32 * FPMSCALEFACTOR;

figure('Color', 'white');
plot(0.5*fpmarraysize + hexX * hexsteppix, 0.5*fpmarraysize + hexY * hexsteppix, '+');
axis equal
axis xy
%% Duplicate Olivier's hex reconstruction
% first render hexes onto the fpm 2048 x 2048 construction array, which
% contains hex number.  This rendering is order dependent.  
vPlus = [sqrt(3)/2, 0.5]; % [cos(pi/6), sin(pi/6)]
vMinus = [sqrt(3)/2, -0.5]; % [cos(-pi/6), sin(-pi/6)]
hexgap = -0.0001;

hexXFpmIndexCoords = 0.5*fpmarraysize + hexXTest * hexsteppix;
hexYFpmIndexCoords = 0.5*fpmarraysize + hexYTest * hexsteppix;
hexRadius = hexsteppix*(1.0-hexgap)*(sqrt(3.0)/2.0); % inner radius when hexes at same y are separated in x by 3 units
hexOuterRadius = 2*hexRadius/sqrt(3);
closeHexSq = hexOuterRadius^2;

fpmDesignArray = -1*ones(fpmarraysize, fpmarraysize);
fpmDesignSags = zeros(fpmarraysize, fpmarraysize);
for h = 1:length(hexNumTest)
    disp(['fpmDesignArray: h = ' num2str(h) ' of ' num2str(length(hexNumTest))]);
    for ii=1:fpmarraysize
        for jj=1:fpmarraysize
            dx = 1.0*jj-hexXFpmIndexCoords(h);
            dy = 1.0*ii-hexYFpmIndexCoords(h);
            r = dx.*dx + dy.*dy;
            if r > closeHexSq
                continue;
            end
            
            dv = [dx, dy];
            if abs(dy) < hexRadius ...
                && abs(dv*vPlus') < hexRadius ...
                && abs(dv*vMinus') < hexRadius;
                fpmDesignArray(ii,jj) = hexNumTest(h);
                fpmDesignSags(ii,jj) = fpmSags(h);
            end            
        end
    end
end

%%
% Map the hex creation array to the mask array in fpm from PIAACMCsimul.c starting line 2730:
% size = 1024;
% NBsubPix = 64;
% fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[0]*piaacmc[0].Fratio;
%           = 1.05386e-06
% piaacmc[0].fpmRad = 0.0001356
%
% we're going to subsample the fpm array by NBsubPix
% each grid node (ii, jj) + (iii, jjj) in fpm is given the physical coordinates
% x = (1.0*ii - size/2 + 1.0*(0.5+iii)/NBsubPix-0.5)*fpscale;
% y = (1.0*jj - size/2 + 1.0*(0.5+jjj)/NBsubPix-0.5)*fpscale;
% this is mapped to an index into the hexagon creation array as 
% ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);
% jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);

fpmSize = 1024;
% fpscale = 1.05386e-06;
% fpscale = 1.098e-06
lambda = 5.62175e-07;
pscale = 0.00011;
beamrad = 0.022;
fpscale = 2*(beamrad/(pscale*fpmSize*16))*lambda*80;
NBsubPix = 64;
fpmRad = 0.0001356;
ac1 = -1 - fpmSize/2 - 0.5;
xc1 = 0.5/fpmRad*FPMSCALEFACTOR;

testSagArray = zeros(fpmSize);
nSubPix = zeros(fpmSize, fpmSize, 3);
sagVals = zeros(fpmSize, fpmSize, 3);
subHexNum = zeros(fpmSize, fpmSize, 3);
for i=1:size(testSagArray, 1)
    disp(['i = ' num2str(i)]);
    iii = fix(NBsubPix/2);
    x = (i-1 - fpmSize/2 + (0.5+iii)/NBsubPix-0.5)*fpscale;
    disp(['x = ' num2str(x)]);
    ii1 = fix((0.5 + 0.5*x/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5);
    if ii1 < 1
        disp('off C');
    end
    for j=1:size(testSagArray, 2)
        jjj = fix(NBsubPix/2);
        y = (j-1 - fpmSize/2 + (0.5+jjj)/NBsubPix-0.5)*fpscale;
        jj1 = fix((0.5 + 0.5*y/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5);
        if ii1 < 1 || ii1 > fpmarraysize || jj1 < 1 || jj1 > fpmarraysize
            continue;
        end
        
        hexList = [];
        for iii = 0:NBsubPix-1
            for jjj = 0:NBsubPix - 1
                % physical coordinates of each subpixel
                y = (j + ac1 + (0.5+jjj)/NBsubPix)*fpscale;
                x = (i + ac1 + (0.5+iii)/NBsubPix)*fpscale;
                % coordinates of each subpixel on the hexagon creation array
                jj1 = fix((0.5 + xc1*y)*fpmarraysize + 0.5);
                ii1 = fix((0.5 + xc1*x)*fpmarraysize + 0.5);
                if ii1 < 1 || ii1 > fpmarraysize || jj1 < 1 || jj1 > fpmarraysize
                    nSubPix(i, j, 1) = nSubPix(i, j, 1) + 1;
                    continue;
                end

                inHexIdx = fpmDesignArray(ii1, jj1) + 1; % hex number + 1, indexes into sag array
                if inHexIdx < 1 % if this is not on a hex
                    nSubPix(i, j, 1) = nSubPix(i, j, 1) + 1;
                    continue;
                end
                if ~ismember(inHexIdx, hexList)
                    hexList = [hexList, inHexIdx]; % add this hex index to the list of the three
                    if length(hexList) > 3
                        error('too many hexes');
                    end
                end
                hexListIdx = find(hexList == inHexIdx); % find which of the three is this hex index
                if ~isempty(hexListIdx)
                    nSubPix(i, j, hexListIdx) = nSubPix(i, j, hexListIdx) + 1; % add to the subPixel count
                    sagVals(i, j, hexListIdx) = fpmSags(inHexIdx); % record the sag for this value
                    subHexNum(i, j, hexListIdx) = inHexIdx - 1; % record the hex number
                end
            end
        end
    end
end
return
%%
testSagArray = sum(nSubPix.*sagVals, 3)./(NBsubPix*NBsubPix);

figure('Color', 'white');
imagesc(testSagArray);
axis equal
axis xy

save hexdata_olivier_1.mat nSubPix sagVals subHexNum testSagArray fpscale lambda NBsubPix hexXFpmIndexCoords hexYFpmIndexCoords fpmDesignArray fpmDesignSags

%% intpolate amplitude and phase
lambda = 5.50875e-07;
ssp = sum(nSubPix,3);
for i=1:size(ssp, 1)
    for j=1:size(ssp, 2);
        if ssp(i,j) == 0
            nSubPix(i,j,1) = NBsubPix*NBsubPix;
        end
    end
end
testFPMInterp = sum(nSubPix.*(exp(2*pi*1i*2*sagVals/lambda)), 3)./(NBsubPix*NBsubPix);

testFPMInterpPh = angle(testFPMInterp);
testFPMInterpAmp = abs(testFPMInterp);

fitswrite(testFPMInterpPh, ['testFPMInterpPh_' num2str(lambda) '.fits']);
fitswrite(testFPMInterpAmp, ['testFPMInterpAmp_' num2str(lambda) '.fits']);

return
%% make hex from exact physical space
% Map the hex creation array to the mask array in fpm from PIAACMCsimul.c starting line 2730:
% size = 1024;
% NBsubPix = 64;
% fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[0]*piaacmc[0].Fratio;
%           = 1.05386e-06
% piaacmc[0].fpmRad = 0.0001356
%
% we're going to subsample the fpm array by NBsubPix
% each grid node (ii, jj) + (iii, jjj) in fpm is given the physical coordinates
% x = (1.0*ii - size/2 + 1.0*(0.5+iii)/NBsubPix-0.5)*fpscale;
% y = (1.0*jj - size/2 + 1.0*(0.5+jjj)/NBsubPix-0.5)*fpscale;
% this is mapped to an index into the hexagon creation array as 
% ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);
% jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);


fpmSize = 1024;
% fpscale = 1.05386e-06;
% fpscale = 1.098e-06
lambda = 5.62175e-07;
pscale = 0.00011;
beamrad = 0.022;
fpscale = 2*(beamrad/(pscale*fpmSize*16))*lambda*80;
NBsubPix = 64;
fpmRad = 0.0001356;
vPlus = [sqrt(3)/2, 0.5]; % [cos(pi/6), sin(pi/6)]
vMinus = [sqrt(3)/2, -0.5]; % [cos(-pi/6), sin(-pi/6)]
hexgap = -0.0001;
hexRadius = hexsteppix*(1.0-hexgap)*(sqrt(3.0)/2.0); % inner radius when hexes at same y are separated in x by 3 units
hexOuterRadius = 2*hexRadius/sqrt(3);
hexXFpmIndexCoords = 0.5*fpmarraysize + hexXTest * hexsteppix;
hexYFpmIndexCoords = 0.5*fpmarraysize + hexYTest * hexsteppix;
ac1 = -1 - fpmSize/2 - 0.5;
xc1 = 0.5/fpmRad*FPMSCALEFACTOR;
closeHexSq = hexOuterRadius^2;

testSagArray = zeros(fpmSize);
nSubPix = zeros(fpmSize, fpmSize, 3);
sagVals = zeros(fpmSize, fpmSize, 3);
for i=1:size(testSagArray, 1)
    disp(['i = ' num2str(i)]);
    for j=1:size(testSagArray, 2)
        iii = fix(NBsubPix/2);
        jjj = fix(NBsubPix/2);
        x = (j-1 - fpmSize/2 + (0.5+jjj)/NBsubPix-0.5)*fpscale;
        y = (i-1 - fpmSize/2 + (0.5+iii)/NBsubPix-0.5)*fpscale;
        jj1 = (0.5 + 0.5*x/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
        ii1 = (0.5 + 0.5*y/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
        if ii1 < 1 || ii1 > fpmarraysize || jj1 < 1 || jj1 > fpmarraysize
            continue;
        end
        
        % sub pixel row and column coordinates
        iii = 0:NBsubPix-1;
        jjj = 0:NBsubPix - 1;
        % physical coordinates of each subpixel
        x = (j + ac1 + (0.5+jjj)/NBsubPix)*fpscale;
        y = (i + ac1 + (0.5+iii)/NBsubPix)*fpscale;
        % coordinates of each subpixel on the hexagon creation array
        jj1 = (0.5 + xc1*x)*fpmarraysize + 0.5;
        ii1 = (0.5 + xc1*y)*fpmarraysize + 0.5;
        
        % grids of coodinates on the hexagon creation array
        [jj1Grid, ii1Grid] = meshgrid(jj1, ii1);
%         % make a linearlized list for later multiplying by vectors to
%         % see if it's in the interior of a hex
%         subGridCoords = [jj1Grid(:), ii1Grid(:)]; 
        
        % collect how many sub-pixels are inside a hexagon along with the 
        % sag of that hexagon for this i,j
        % There can be at most three hexagons intersecting this i,j
        % loop through the hex centers
        hexCount = 1;
        for h=1:length(hexXFpmIndexCoords)
            % find the subgrid points close enough to each hex
            dx = jj1Grid - hexXFpmIndexCoords(h);
            dy = ii1Grid - hexYFpmIndexCoords(h);
            dv = [dx(:), dy(:)];
            r = dx.*dx + dy.*dy;
            
            if all(r > closeHexSq)
                continue;
            end
            % we're in a hex that intersects some sub-pixels at this (i,j)
            % Now find the sub-pixels that are actually in this hex
            inHexIdx = find(abs(dy(:)) < hexRadius ...
                & abs(dv*vPlus') < hexRadius ...
                & abs(dv*vMinus') < hexRadius);
            
            % record the number of sub-pixels are inside this hexagon
            nSubPix(i, j, hexCount) = length(inHexIdx);
            sagVals(i, j, hexCount) = fpmSags(h);
            hexCount = hexCount + 1;
            if hexCount > 3
                break;
            end
        end
    end
end
testSagArray = sum(nSubPix.*sagVals, 3)./(NBsubPix*NBsubPix);

figure('Color', 'white');
imagesc(testSagArray);
axis equal
axis xy


return

%%
% Map the hex creation array to the mask array in fpm from PIAACMCsimul.c starting line 2730:
% size = 1024;
% NBsubPix = 64;
% fpscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size/piaacmc[0].fpzfactor*optsyst[0].lambdaarray[0]*piaacmc[0].Fratio;
%           = 1.05386e-06
% piaacmc[0].fpmRad = 0.0001356
%
% we're going to subsample the fpm array by NBsubPix
% each grid node (ii, jj) + (iii, jjj) in fpm is given the physical coordinates
% x = (1.0*ii - size/2 + 1.0*(0.5+iii)/NBsubPix-0.5)*fpscale;
% y = (1.0*jj - size/2 + 1.0*(0.5+jjj)/NBsubPix-0.5)*fpscale;
% this is mapped to an index into the hexagon creation array as 
% ii1 = (long) ( (0.5 + 0.5*x/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);
% jj1 = (long) ( (0.5 + 0.5*y/piaacmc[0].fpmRad*FPMSCALEFACTOR)*piaacmc[0].fpmarraysize + 0.5);


fpmSize = 1024;
% fpscale = 1.05386e-06;
% fpscale = 1.098e-06
lambda = 5.62175e-07;
pscale = 0.00011;
beamrad = 0.022;
fpscale = 2*(beamrad/(pscale*fpmSize*16))*lambda*80;
NBsubPix = 64;
fpmRad = 0.0001356;
vPlus = [sqrt(3)/2, 0.5]; % [cos(pi/6), sin(pi/6)]
vMinus = [sqrt(3)/2, -0.5]; % [cos(-pi/6), sin(-pi/6)]
hexgap = -0.0001;
hexRadius = hexsteppix*(1.0-hexgap)*(sqrt(3.0)/2.0); % inner radius when hexes at same y are separated in x by 3 units
hexXFpmIndexCoords = 0.5*fpmarraysize + hexXTest * hexsteppix;
hexYFpmIndexCoords = 0.5*fpmarraysize + hexYTest * hexsteppix;
ac1 = -1 - fpmSize/2 - 0.5;
xc1 = 0.5/fpmRad*FPMSCALEFACTOR;

testSagArray = zeros(fpmSize);
testSagArray2 = zeros(fpmSize);
for i=1:size(testSagArray, 1)
    disp(['i = ' num2str(i)]);
    for j=1:size(testSagArray, 2)
        iii = fix(NBsubPix/2);
        jjj = fix(NBsubPix/2);
        x = (j-1 - fpmSize/2 + (0.5+jjj)/NBsubPix-0.5)*fpscale;
        y = (i-1 - fpmSize/2 + (0.5+iii)/NBsubPix-0.5)*fpscale;
        jj1 = (0.5 + 0.5*x/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
        ii1 = (0.5 + 0.5*y/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
        if ii1 < 1 || ii1 > fpmarraysize || jj1 < 1 || jj1 > fpmarraysize
            continue;
        end
        for iii = 0:NBsubPix-1
            for jjj = 0:NBsubPix - 1
%                 x = (j-1 - fpmSize/2 + (0.5+jjj)/NBsubPix-0.5)*fpscale;
%                 y = (i-1 - fpmSize/2 + (0.5+iii)/NBsubPix-0.5)*fpscale;
                x = (j + ac1 + (0.5+jjj)/NBsubPix)*fpscale;
                y = (i + ac1 + (0.5+iii)/NBsubPix)*fpscale;
%                 jj1 = (0.5 + 0.5*x/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
%                 ii1 = (0.5 + 0.5*y/fpmRad*FPMSCALEFACTOR)*fpmarraysize + 0.5;
                jj1 = (0.5 + xc1*x)*fpmarraysize + 0.5;
                ii1 = (0.5 + xc1*y)*fpmarraysize + 0.5;
                if ii1 < 1 || ii1 > fpmarraysize || jj1 < 1 || jj1 > fpmarraysize
                    continue;
                end
                d = sqrt((jj1 - hexXFpmIndexCoords).^2 + (ii1 - hexYFpmIndexCoords).^2);
                [dum, zoneNumber] = min(d);
                if hexNumTest(zoneNumber) == -1
                    continue;
                end
                testSagArray2(i,j) = testSagArray2(i,j) + fpmSags(zoneNumber);
                dv = [jj1 - hexXFpmIndexCoords(zoneNumber), ii1 - hexYFpmIndexCoords(zoneNumber)];
                if abs(dv(2)) < hexRadius && abs(dv*vPlus') < hexRadius && abs(dv*vMinus') < hexRadius
                    testSagArray(i,j) = testSagArray(i,j) + fpmSags(zoneNumber);
                end
            end
        end
%         if i == 512 && j == 512
%             keyboard;
%         end
    end
end
testSagArray = testSagArray./(NBsubPix*NBsubPix);
testSagArray2 = testSagArray2./(NBsubPix*NBsubPix);

figure('Color', 'white');
imagesc(testSagArray);
axis equal
axis xy



