%% 2016_12_11
% segment holes in cytosolic filled neurons
% see GU_Expansion.m
rt =  'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\';
cd(rt)
im = readtiff([rt 'c0.tif']);
% im2 = readtiff([rt 'c1.tif']);
% im3 = readtiff([rt 'c2.tif']);
px = .097;
pz = .18;
ef = 3.95;
zAniso = pz/px;
edgepxX = [5, 5]; % to remove
edgepxY = [5, 5]; % to remove
edgepxZ = [30, 30]; % to remove
MinHoleVolume = 50;
%%
% segment cellbody
% opts = {'GaussSigma', 3, 'PixelSize', px, 'PixelSize', pz,'LargestObject', false,'DistanceTransform', false,'Otsu', false,'Bk', 1700,'sdMultiplier', 0};
% [filledDM, unfilledSeg, T] = GU_segmentCytosolicVol3D(im, opts{:});
% filledDM = single(filledDM);

% % %
filledDM = readtiff([rt 'Mask_Filled_T1700.tif']);
unfilledSeg = readtiff([rt 'Mask_T1700.tif']);
%

holes = logical(filledDM.*(~logical(unfilledSeg)));
writetiff(uint8(holes) , [rt 'Mask_holes.tif']);
% writetiff(mask , [rt 'MaskedData.tif']);
% % mask = readtiff([rt 'MaskedData.tif']);
% holes = filterGauss3D(mask,2) < T;
% holes = logical(filledDM>3).*holes;


holes(:,:,[1:25,end-25:end])=0;
writetiff(uint8(holes) , [rt 'Mask_holes_z26-404.tif']);

[lholes, ~] = bwlabeln(holes);
clear holes
% 
% discard small holes
CC = bwconncomp(lholes, 26);
csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinHoleVolume;
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
lholes = labelmatrix(CC)~=0;
% 
[lholes, nholes] = bwlabeln(logical(lholes));
% if nholes < 2^16-1
%     lholes = uint16(lholes);
% else
%     lholes = single(lholes);
% end
% 
% clear mask
% lholes = readtiff([rt 'holes.tif']);
cellidx = unique(lholes(:));
    ucl = unique(lholes(:));
ucl = ucl(ucl>0);
rmidx =unique([unique(lholes(:,:,[1:edgepxZ(1), end-edgepxZ(2):end])); unique(lholes(:,[1:edgepxX(1), end-edgepxX(2):end],:)); unique(lholes([1:edgepxY(1), end-edgepxY(2):end],:,:))]);

% remove edge labels
idx = ~ismember(cellidx,rmidx);
cellidx = cellidx(idx);
cellidx = sort(cellidx);
lholes(ismember(cellidx,rmidx)) = 0;
writetiff(uint16(lholes) , [rt 'holes_T1700.tif']);


%%
lholes = readtiff([rt 'holes_T1700_withOldMaskToRemoveCRAP.tif']);
manualRMIDX = [2515,2511,549,2377,1250,2476,1432,2332,1559,1939,1915,1938,2086,2018,2133,2406,2260];
lholes(ismember(lholes(:),manualRMIDX)) = 0;
idx = lholes>0;
lx = lholes(idx);
nlx = zeros(size(lx), 'single');
%
u = unique(lx(lx>0));
parfor_progress(numel(u));
for i = 1:numel(u)
    npx = lx==u(i);
   nlx(npx) = single(log10(sum(npx)*(px^2*pz)/ef^3*10^5));
   parfor_progress;
end
parfor_progress(0)

%
volHoles = zeros(size(lholes), 'single');
volHoles(idx) = nlx;


% oldHoles = logical(readtiff([rt 'holes.tif']));
% volHoles(~oldHoles) = 0;
writetiff(volHoles, [rt 'holes_T1700_volumecoded_x1000log10_manualCurrated_corrected.tif']);
%%
edgepxZ = [2,2];
edgepxX = [2,2];
edgepxY = [2,2];

lholes = readtiff([rt 'holes_T1700_withOldMaskToRemoveCRAP.tif']);
 rmidx = unique([unique(lholes(:,:,[1:edgepxZ(1), end-edgepxZ(2):end])); unique(lholes(:,[1:edgepxX(1), end-edgepxX(2):end],:)); unique(lholes([1:edgepxY(1), end-edgepxY(2):end],:,:))]);
% get all cell labels
cellidx = unique(lholes(:));
% remove edge labels
idx = ~ismember(cellidx,rmidx);
cellidx = cellidx(idx);
cellidx = sort(cellidx);
lholes(~ismember(lholes,cellidx)) = 0; % set the edge labels to 0

CC = bwconncomp(lholes, 26);
csize = cellfun(@numel, CC.PixelIdxList);
holeVol = csize*0.097^2*0.18;%%%%%%% calc volume of mito holes
rp3 = regionprops3(logical(lholes));

 save([rt date 'parameters_rp3_allholes'], 'holeVol', 'rp3')
%% initialize vectors

lholes(~oldHoles) = 0;
writetiff(uint16(lholes) , [rt 'holes_T1700_withOldMaskToRemoveCRAP.tif']);
ucl = unique(lholes);
nc = numel(ucl);%numel(cellidx);
holevol = zeros(nc,1);
holearea = zeros(nc,1);
minDistFromMem = zeros(nc,1);
meanDistFromMem = zeros(nc,1);
maxDistFromMem = zeros(nc,1);
cX = zeros(nc,1);
cY = zeros(nc,1);
cZ = zeros(nc,1);
% loop through each of the cells in the volume
tic
parfor_progress(nc);
for j = 1:nc % par
    mask = lholes==ucl(j);
    rp3(j)  = regionprops3(mask);
%     dt = filledDM(mask);
%     ch2pxInt(j) = sum(im2(mask>0));
%      ch3pxInt(j) = sum(im3(mask>0));
%     minDistFromMem(j) = min(dt(dt>0));
%     maxDistFromMem(j) = max(dt(dt>0));
%     meanDistFromMem(j) = mean(dt(dt>0));
    idx = find(mask);
    [ny,nx,nz] = size(mask);
    [yi,xi,zi] = ind2sub([ny,nx,nz], idx);
    b = 5; % boundary pixels for cropping
    xa = max(min(xi)-b,1):min(max(xi)+b,nx);
    ya = max(min(yi)-b,1):min(max(yi)+b,ny);
    za = max(min(zi)-b,1):min(max(zi)+b,nz);
    mask = mask(ya,xa,za);
    % interpolate
    [ny,nx,nz] = size(mask);
    [y,x,z] = ndgrid(1:ny,1:nx,1:nz);
    [Y,X,Z] = ndgrid(1:ny,1:nx,1:1/zAniso:nz);
    imask = interp3(x,y,z,double(mask),X,Y,Z,'nearest');
    perim = bwperim(imask);
    holevol(j) = sum(imask(:)) * px^3;
    holearea(j) = sum(perim(:)) * px^2;
    s  = regionprops(imask,'centroid');
    centroids = cat(1, s.Centroid);
    cX(j) = centroids(:,1);
    cY(j) = centroids(:,2);
    cZ(j) = centroids(:,3);
    parfor_progress;
end
toc
parfor_progress(0);

%write parmas
% params.minDistFromMem = minDistFromMem;
% params.meanDistFromMem = meanDistFromMem;
% params.maxDistFromMem = maxDistFromMem;
params.intCh2Int = ch2pxInt;
params.intCh3Int = ch3pxInt;
params.holeVol = holevol;
params.holeArea = holearea;
params.holeLabel = ucl;
params.cX = cX;
params.cY = cY;
params.cZ = cZ;
params.rp3 = rp3;

save([rt filesep date 'holeParameters_rp3.mat'], 'params','rp3');
% params.ThresholdOtsu = T;
% save([rt filesep 'holeParameters.mat'], 'params');
