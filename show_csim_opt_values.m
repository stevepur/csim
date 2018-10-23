figure;

tgtGeomHexNum = fitsread('../../../donutPiaa/mask_data_hexNum.fits');
tgtGeomNSubPix = fitsread('../../../donutPiaa/mask_data_nSubPix.fits');
% tgtGeomHexNum = fitsread('../../wfirst_base_mask_data_1024_hexNum.fits');
% tgtGeomNSubPix = fitsread('../../wfirst_base_mask_data_1024_nSubPix.fits');
NBsubPix = 64;
tgtGeomSagVals = zeros(size(tgtGeomHexNum));
while (1) 

    fid = fopen('optValHistory.txt', 'r');
    a = fscanf(fid, '%d %f');
    fclose(fid);
    
    subplot(1,3,1);
    c = a(3:2:end);
    v = a(4:2:end);
    semilogy(c,v,'ro-');
    grid on;
    ylim([prctile(v, 0), prctile(v, 99.5)]);
    ax = axis;

    fid = fopen('optBestVal.txt', 'r');
    a = fscanf(fid, '%d %f');
    fclose(fid);
    
    subplot(1,3,2);
    c = a(3:2:end);
    v = a(4:2:end);
    ii = find(c < 100);
    if isempty(ii)
        ii = 0;
    end
    semilogy(c(ii+1:end),v(ii+1:end),'ro-');
    ylim([min(v), prctile(v, 50)]);
    title([num2str(v(end))]);
    grid on;

    tgtSags = fitsread('bestOptVector.fits');

    for r = 1:size(tgtGeomHexNum, 1)
        for c = 1:size(tgtGeomHexNum, 2)
            for s = 1:size(tgtGeomHexNum, 3)
                if (tgtGeomHexNum(r, c, s) > -1)
                    tgtGeomSagVals(r, c, s) = tgtSags(tgtGeomHexNum(r, c, s) + 1); % hexNum is indexed from 0
                else
                    tgtGeomSagVals(r, c, s) = 0;
                end
            end
        end
    end
    tgtMask = sum(tgtGeomNSubPix.*tgtGeomSagVals, 3)./(NBsubPix*NBsubPix);

    subplot(1,3,3);
    imagesc(tgtMask);
    caxis([-5e-7, 5e-7]);
    title('target mask');
    axis equal;
    axis tight;
    
    pause(10);
end