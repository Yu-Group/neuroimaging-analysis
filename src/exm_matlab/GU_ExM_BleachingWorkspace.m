%% GU_ExM_BleachingWorkspace

pr.OTSUGaussKernel = 0.5;
pr.OTSUMaxPer = 99;
MinVoxelVolume = 246;
MaxVoxelVolume = 1500000;
% spinning disk, Airy Scan
rt = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Raw/488/unmixed/';
% rt = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/488/unmixed/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};
imf1 = readtiff([rt fnlist{1}]);
im = zeros([size(imf1), numel(fnlist)], 'uint16');

parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    im(:,:,i) = readtiff([rt fnlist{i}]);
    parfor_progress;
end
parfor_progress(0);

% filter and threshold
tic
imG = filterGauss3D(im(1:3000, 1:3000,:),pr.OTSUGaussKernel);
toc
PP = prctile(im(:),pr.OTSUMaxPer);
T = thresholdOtsu(im(im > 0 & im < PP));
imG = imG-T;
imG(imG<0) = 0;

imG = uint16(imG);


tic
% writetiff(imG, '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/unmixedCleanedSD_sig5.tif')
writetiff(imG, '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/unmixedCleanedSD_sig5.tif');
toc

CC = bwconncomp(logical(imG), 26);
csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinVoxelVolume & csize<=MaxVoxelVolume;
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
imGL = labelmatrix(CC)~=0;
imG = uint16(imG) .* uint16(imGL);



%% read spinning disc data
% imG = readtiff('/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/unmixedCleanedSD.tif');
imG = readtiff('/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/unmixedCleanedSD_sig5.tif');
%%
MaxI = zeros(1,size(imG,3));
AvgI = zeros(1,size(imG,3));
nVox = zeros(1,size(imG,3));
SumI = zeros(1,size(imG,3));
parfor_progress(size(imG,3));
parfor i = 1:size(imG,3)
    slice = imG(:,:,i);
    slice(slice<0) = 0;
    %     slice = slice(slice>0);
    if ~isempty(slice)
        %         MaxI(i) = max(slice(slice>0));
        %         AvgI(i) = mean(slice(slice>0));
        %         SumI(i) = sum(slice(slice>0));
        %         nVox(i) = sum(slice(:)>0);
        
        MaxI(i) = max(slice(:));
        AvgI(i) = mean(slice(:));
        SumI(i) = sum(slice(:));
        nVox(i) = sum(slice(:)>0);
    end
    parfor_progress;
end
parfor_progress(0);
FN = 1:size(imG,3);
MI = MaxI;%SumI./nVox;
AI = AvgI;%SumI./nVox;
sf = 10;
figure, plot(FN, smooth(AI,sf)/max(smooth(AI,sf))), hold on
plot(FN, smooth(MI,sf)/max(smooth(MI,sf))), hold on
plot(smooth(SumI./nVox,sf)/max(smooth(SumI./nVox,sf)))

%% AiryScan

rt = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/Individualtiles/488/unmixed/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name}';
MaxVoxelVolume = 8800000;
for i = 1:numel(fnlist)
    fn = [rt fnlist{i}];
    submit = ['bsub -n 12 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "bleach' num2str(i) '" -o ~/logs/bleach' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_BleachingOverStack ' fn ' ' num2str(MaxVoxelVolume) ''''];
    [~, ~] = system(submit,'-echo');
end


%% Airy bleaching per imaged volume analysis
rt = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/Individualtiles/488/unmixed/';
fnlist = dir([rt '*.mat']);
fnlist = {fnlist.name}';

% get each tile size
parfor_progress(numel(fnlist))
for i = 1:numel(fnlist)
    load([rt fnlist{i}]);
    f = fieldnames(params);
    for n = 1:numel(f)
        allParams(i).(f{n}) = params.(f{n});
    end
    parfor_progress;
end
parfor_progress(0)
%%
load('/groups/betzig/betziglab/4Stephan/171105_unmixing/Airy_bleachingPerTile.mat')

%
tileSize = [1948, 1948, 1525];
meanSumInt = arrayfun(@(x) mean(x.SumI), allParams, 'unif',1);
meanAvgInt = arrayfun(@(x) mean(x.AvgI), allParams, 'unif',1);
meanMaxInt = arrayfun(@(x) mean(x.MaxI), allParams, 'unif',1);
meanVox = arrayfun(@(x) mean(x.nVox), allParams, 'unif',1);
allVox = arrayfun(@(x) sum(x.nVox), allParams, 'unif',1);
meanAvgInt_2 = arrayfun(@(x) mean(x.AvgI_ofValsMorethanZero), allParams, 'unif',1);

sf = 1;
figure,
plot(smooth(meanSumInt./meanVox,sf)/max(smooth(meanSumInt./meanVox,sf))); hold on
% plot(smooth(meanSumInt,sf)/max(smooth(meanSumInt,sf))); hold on
plot(smooth(meanAvgInt_2,sf)/max(smooth(meanAvgInt_2,sf))); hold on
% plot(smooth(meanMaxInt,sf)/max(smooth(meanMaxInt,sf))); hold on
plot(smooth(allVox,sf)/max(smooth(allVox,sf))); hold on


%% llsm

rt = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/matlab_decon/unmixed/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name}';
MaxVoxelVolume = 3800000;
for i = 1:numel(fnlist)
    fn = [rt fnlist{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "llsmBleach' num2str(i) '" -o ~/logs/llsmBleach' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_BleachingOverStack ' fn ' ' num2str(MaxVoxelVolume) ''''];
    [~, ~] = system(submit,'-echo');
end
%% LLSM bleaching per imaged volume analysis
rt = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/matlab_decon/unmixed/';
fnlist = dir([rt '*.mat']);
fnlist = {fnlist.name}';

% get each tile size


parfor_progress(numel(fnlist))
for i = 1:numel(fnlist)
    load([rt fnlist{i}]);
    f = fieldnames(params);
    
    for n = 1:numel(f)
        allParams(i).(f{n}) = params.(f{n});
    end
    parfor_progress;
end
parfor_progress(0)

%% LLSM - stage sweep mode
% data is deskewed, deconvolved and rotated
rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/allRenamedTiles/export.n5';
opath = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/allRenamedTiles/n5ExtractedSubVols/';

if ~exist(opath)
    mkdir(opath)
end

mx = '11393'; %1
my = '10663'; %7
mz = '5301'; % 41
csx = '5697';
csy = '1524';
csz = '130';
OL = '0';
ch = '0';
% subPanelIdx = 1;

% GU_ExM_JaneliaCluster_Extract_n5_subVols(rt, opath, mx,my,mz,csx,csy,csz, OL,ch, subPanelIdx)
% crop data
for i = 1:574
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "llsmCrop' num2str(i) '" -o ~/logs/llsmCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Extract_n5_subVols ' rt ' ' opath ' ' mx ' ' my ' ' mz ' ' csx ' ' csy ' ' csz ' ' OL ' ' ch ' ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end

%% run llSM stage sweep bleach calculations
fnlist = dir([opath '*.tif']);
fnlist = {fnlist.name}';
MaxVoxelVolume = 3040000;
for i = 1:numel(fnlist)
    fn = [opath fnlist{i}];
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "llsmBleach' num2str(i) '" -o ~/logs/llsmBleach' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_BleachingOverStack ' fn ' ' num2str(MaxVoxelVolume) ''''];
    [~, ~] = system(submit,'-echo');
end
% % test
% i = 72
% fn = [opath fnlist{i}];
% GU_BleachingOverStack(fn, num2str(MaxVoxelVolume))

%% merge params
fnlist = dir([opath '*.tif']);
fnlist = {fnlist.name}';
paramRT = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/allRenamedTiles/n5ExtractedSubVols/BleachParams/';
parfor_progress(numel(fnlist));
for i = 1:numel(fnlist)
    fn = [opath fnlist{i}];
    [p, f, ~] = fileparts(fn);
    nfn = [paramRT f '.mat'];
    if exist(nfn, 'file')
        load(nfn)
        allParams(i) = params;
    end
    parfor_progress;
end
parfor_progress(0);
%%
load('/groups/betzig/betziglab/4Stephan/171105_unmixing/LLSM_bleachingPerTile.mat')
%%
%
tileSize = [360, 704, 501];
meanSumInt = arrayfun(@(x) mean(x.SumI), allParams, 'unif',1);
meanAvgInt = arrayfun(@(x) mean(x.AvgI), allParams, 'unif',1);
meanMaxInt = arrayfun(@(x) mean(x.MaxI), allParams, 'unif',1);
meanVox = arrayfun(@(x) mean(x.nVox), allParams, 'unif',1);
allVox = arrayfun(@(x) sum(x.nVox), allParams, 'unif',1);
meanAvgInt_2 = arrayfun(@(x) mean(x.AvgI_ofValsMorethanZero), allParams, 'unif',1);

sf = 10;
figure,
plot(smooth(meanSumInt./meanVox,sf)/max(smooth(meanSumInt./meanVox,sf))); hold on
% plot(smooth(meanSumInt,sf)/max(smooth(meanSumInt,sf))); hold on
plot(smooth(meanAvgInt_2,sf)/max(smooth(meanAvgInt_2,sf))); hold on
% plot(smooth(meanMaxInt,sf)/max(smooth(meanMaxInt,sf))); hold on
plot(smooth(allVox,sf)/max(smooth(allVox,sf))); hold on

%% spinning disk

rt = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/BleachingAnalysis_SD/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name}';
MaxVoxelVolume = 1000000;
fn = [rt fnlist{1}];
submit = ['bsub -n 32 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "sdbleach' num2str(i) '" -o ~/logs/sdbleach' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_BleachingOverStack ' fn ' ' num2str(MaxVoxelVolume) ''''];
[~, ~] = system(submit,'-echo');


%% crop regions for Figure 1 - need MIPs along the xy. yz. xz

%LLSM(O):
rtLLSMO = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch0/unmixed/';
rt = rtLLSMO;
cd(rt);
fname = dir([rt '*.tif']);
fname = {fname.name}';
% generate cropped tifs
xrmin = 2599;
xrmax = 3531;
yrmin = 1;
yrmax = 5461;
k = 1;
for i = 1:numel(fname)
    fn = [rt fname{i}];
    %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
    [stat, res] = system(submit,'-echo');
end

% LLSM(SS):
rtLLSMss = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/allRenamedTiles/slice-tiff/ch0/';
rt = rtLLSMss;
cd(rt);
fname = dir([rt '*.tif']);
fname = {fname.name}';
% generate cropped tifs
xrmin = 5328;
xrmax = 6054;
yrmin = 1;
yrmax = 10663;
k = 1;
for i = 1:numel(fname)
    fn = [rt fname{i}];
    %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
    [stat, res] = system(submit,'-echo');
end

% SD:
rtSD = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Deconed/488/';
rt = rtSD;
cd(rt);
fname = dir([rt '*.tif']);
fname = {fname.name}';
% generate cropped tifs
xrmin = 1;
xrmax = 3584;
yrmin = 1510;
yrmax = 2074;
k = 1;
for i = 1:numel(fname)
    fn = [rt fname{i}];
    %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
    [stat, res] = system(submit,'-echo');
end


% A(F):
rtAS = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/488/unmixed/';
rt = rtAS;
cd(rt);
fname = dir([rt '*.tif']);
fname = {fname.name}';
% generate cropped tifs
xrmin = 1;
xrmax = 5788;
yrmin = 875;
yrmax = 2625;
k = 1;
for i = 1:numel(fname)
    fn = [rt fname{i}];
    %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
    [stat, res] = system(submit,'-echo');
end
%
% LLSM(O):
% \\dm11\betziglab\Ved\Data\20161211_ExM\YFP2_comparison2\Stitch_Igor\slice-tiff\ch0\unmixed
%  
% 
% LLSM(SS): 
% \\dm11\betziglab\4Stephan\171129_YFPsample4\Processing\ch1\DSR_dx0.1229_dz0.16515\allRenamedTiles\slice-tiff\ch0
% 
% 
% SD: 
% \\dm11\betziglab\4Stephan\171105_unmixing\Spinningdisk\Deconed\488\unmixed
% 
% 
% A(F): 
% \\dm11\betziglab\4Stephan\171105_unmixing\Airyscan\488\unmixed
