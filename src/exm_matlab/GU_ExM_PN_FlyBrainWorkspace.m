
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PN data set: nc82 + Syd1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.Alpha = 0.05;
opts.CellMask = [];
opts.RemoveRedundant = true;
opts.WindowSize = [];
opts.Mode = 'xyzAsrc';
sigma = [1.3, 1.5];
% chromatic offset correction
for k = 1:numel(data)
    for j = 1:3
        %         [sigmaXY{k}{j}, sigmaZ{k}{j}, XYZ{k}{j}, pstruct{k}{j}, mask{k}{j}] = GU_estimateSigma3D(data(k).framePaths{1}{j},[]);
        frame = double(readtiff(data(k).framePaths{1}{j})); %#ok<PFBNS>
        frame(frame==0) = NaN; % to avoid border effects
        frame = frame-700;
        [pstruct, mask] = pointSourceDetection3D(frame, sigma, 'Mode', opts.Mode, 'Alpha', opts.Alpha,...
            'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
            'RefineMaskLoG', false, 'WindowSize', opts.WindowSize); %#ok<PFBNS>
        idx = find(pstruct.A==max([pstruct.A]));
        sigmaXY{k}{j} = pstruct.s(1,idx);
        sigmaZ{k}{j} = pstruct.s(2,idx);
        XYZ{k}{j}.x = pstruct.x(idx);
        XYZ{k}{j}.y = pstruct.y(idx);
        XYZ{k}{j}.z = pstruct.z(idx);
    end
    %
end
% z offset
for k = 1:numel(data)
    Off_rfp_alexa(k,:) =  [XYZ{k}{1,2}.x - XYZ{k}{1,3}.x, XYZ{k}{1,2}.y - XYZ{k}{1,3}.y, XYZ{k}{1,2}.z - XYZ{k}{1,3}.z];
    Off_gfp_alexa(k,:) =  [XYZ{k}{1,1}.x - XYZ{k}{1,3}.x, XYZ{k}{1,1}.y - XYZ{k}{1,3}.y, XYZ{k}{1,1}.z - XYZ{k}{1,3}.z];
    Off_gfp_rfp(k,:) =  [XYZ{k}{1,1}.x - XYZ{k}{1,2}.x, XYZ{k}{1,1}.y - XYZ{k}{1,2}.y, XYZ{k}{1,1}.z - XYZ{k}{1,2}.z];
    
end

%%




sigVox(1) = GU_calcGaussianIntegral(1, [2.5,2.5,2.5]);
% fnv1 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch1\channel 1 [8775, 4298, 4954].tif']; % nc82
% fnv2 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch2\channel 2 [8775, 4298, 4954].tif']; % syd1
% fnv3 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch0\channel 0 [8775, 4298, 4954].tif']; % PN

fnv1 = [data.framePaths{2}{1}]; % nc82
fnv2 = [data.framePaths{3}{1}]; % syd1
fnv3 = [data.framePaths{1}{1}]; % PN

%     GU_ExM_zOffset3Dframe(fnv1, 1,1);
%     GU_ExM_zOffset3Dframe(fnv2, 1,6);


% depuncate and calculate the local maxima
for ni = 1:numel(sigVox)
    % nc82 + PN
    GU_ExM_calcPunctaDensityVol(fnv1, fnv3,...
        'FindLocalMaxima', true, 'MinThreshold', [576,756],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99,'minAZdia', 0.1);
    % Syd1 + PN
    GU_ExM_calcPunctaDensityVol(fnv2, fnv3,...
        'FindLocalMaxima', true, 'MinThreshold', [119,756],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99,'minAZdia', 0.1);
end


%% chromatic offset correction
rt = 'D:\Gokul\Wholebrainsynapse\x17_PN+nc82+Syd1_testforchrome_z0p18\Copy_of_data_for_FAKEpsfs_z0p18\';
ch0 = 'ch0[4307, 5619, 5902].tif';
ch1 = 'ch1[4307, 5619, 5902].tif';
ch2 = 'ch2[4307, 5619, 5902].tif';

% GU_ExM_zOffset3Dframe(ch0, 0,6);
% GU_ExM_zOffset3Dframe(ch1, 1,1);
% GU_ExM_zOffset3Dframe(ch1, 0,5);
% GU_ExM_zOffset3Dframe(ch2, 1,6);

im0 = readtiff([rt ch0]);
im1 = readtiff([rt ch1]);
im2 = readtiff([rt ch2]);


%% determine PN + nc82 + Syd1 data's thresholds
MinHoleVolume = 1008;% for membrane
MinHoleVolume = 516; % for puncta

for i = 1%1:numel(data)
    for j = 1:numel(data(i).framePaths)
        fprintf('reading image...')
        tic
        im = readtiff(data(i).framePaths{1}{j});
        toc
        
        fprintf('filtering image...')
        tic
        imG = filterGauss3D(im,1);
        T = thresholdOtsu(im(im > 0 & im < prctile(im(:),99)));
        imG = imG-T;
        imG(imG<0) = 0;
        toc
        
        % discard small holes
        fprintf('discarding small holes...')
        tic
        CC = bwconncomp(logical(imG), 26);
        csize = cellfun(@numel, CC.PixelIdxList);
        idx = csize>=MinHoleVolume;
        CC.NumObjects = sum(idx);
        CC.PixelIdxList = CC.PixelIdxList(idx);
        imGL = labelmatrix(CC)~=0;
        toc
        imG = uint16(imG) .* uint16(imGL);
        [p, ~, ~] = fileparts([data(i).framePathsDS{j}{1}]);
        mkdir(p);
        
        
        fprintf('writing image...')
        tic
        writetiff(imG, [data(i).framePathsDS{j}{1}]);
        
        toc
        
        clear im imG imGL;
    end
end
% check signal intensity of the clean mask in 1px eroded region
% imG = im;
se = strel('sphere',3);
cleanDAN_E = imerode(logical(imG), se);
figure, histogram(imG(imG>0 & ~cleanDAN_E ),1000);
muCut = prctile(imG(imG>0 & ~cleanDAN_E),99.9);

% ch0 - PN threshold = 756-->otsu+99.9 prctle of the edge voxels (from 1)
% ch0 - PN threshold = 1960-->otsu+99.9th of the 2 edge voxels; 4749 for 3px
% ch1 - nc82 - 576
% ch2 - Syd1 - 119
%% submit PN for sub-vol extraction,
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/';
mkdir([rt 'ch0'])
mkdir([rt 'ch1'])
mkdir([rt 'ch2'])
mx = 9664;
my = 6721;
mz = 5186;
OL = 50;
s = repmat(750, [1,3]);

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

% generate x y z cropping coordinates
clear xmin xmax ymin zmin ymax zmax

xmin = zeros(nx*ny*nz,1);
xmax = zeros(nx*ny*nz,1);
ymin = zeros(nx*ny*nz,1);
ymax = zeros(nx*ny*nz,1);
zmin = zeros(nx*ny*nz,1);
zmax = zeros(nx*ny*nz,1);

nn = 1;
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            if i <nx
                xmax(nn) = (s(1)*i-(OL*(i-1)));
            else
                xmax(nn) = mx;
            end
            
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            if j <ny
                ymax(nn) = (s(2)*j-(OL*(j-1)));
            else
                ymax(nn) = my;
            end
            
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            if k < nz
                zmax(nn) = (s(3)*k-(OL*(k-1)));
            else
                zmax(nn) = mz;
            end
            
            nn = nn + 1;
        end
    end
end

for i = 1:nn-1
    submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "PN' num2str(i) '" -o ~/logs/PN' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_PN_ExtractSubVol_ChromaticCorr_Puncta ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end

rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch1/Analysis/246/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};
for ex = 1:nn-1
    fn = [filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_246_clean.tif'];
    CompletedJobsIdx(ex) = ismember(fn(2:end),fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);
for i = job2resubmit
    submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "PN' num2str(i) '" -o ~/logs/PN' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_PN_ExtractSubVol_ChromaticCorr_Puncta ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end

%%
%% combine all the nc82 data - from the recalculated set
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch1/Analysis/246/DataStructures';
rt2 = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch2/Analysis/246/DataStructures';
allnc82filteredData = [];
allnc82Vol2PunctaIdx=[];
allSyd1filteredData = [];
allSyd1Vol2PunctaIdx = [];

parfor_progress(nn-1);
for i = 1:nn-1
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_densityData.mat'];
    % nc82
    variableInfo = who('-file', [rt filesep fn]);
    PF = ismember('filteredData', variableInfo) & ismember('Vol2PunctaIdx', variableInfo);
    if PF
        t = load([rt filesep fn], 'Vol2PunctaIdx', 'filteredData');
        idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
            t.filteredData(:,1)<=725 &  t.filteredData(:,2)<=725 & t.filteredData(:,3)<=725;
        t.filteredData = t.filteredData(idx,:);
        t.filteredData(:,1) = t.filteredData(:,1) + xmin(i)-1;
        t.filteredData(:,2) = t.filteredData(:,2) + ymin(i)-1;
        t.filteredData(:,3) = t.filteredData(:,3) + zmin(i)-1;
        t.Vol2PunctaIdx = t.Vol2PunctaIdx(idx);
        allnc82filteredData = [allnc82filteredData; t.filteredData];
        allnc82Vol2PunctaIdx = logical([allnc82Vol2PunctaIdx; t.Vol2PunctaIdx]);
    end
    clear t variableInfo PF
    
    % Syd1
    variableInfo = who('-file', [rt2 filesep fn]);
    PF = ismember('filteredData', variableInfo) & ismember('Vol2PunctaIdx', variableInfo);
    if PF
        t = load([rt2 filesep fn], 'Vol2PunctaIdx', 'filteredData');
        idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
            t.filteredData(:,1)<=725 &  t.filteredData(:,2)<=725 & t.filteredData(:,3)<=725;
        t.filteredData = t.filteredData(idx,:);
        t.filteredData(:,1) = t.filteredData(:,1) + xmin(i)-1;
        t.filteredData(:,2) = t.filteredData(:,2) + ymin(i)-1;
        t.filteredData(:,3) = t.filteredData(:,3) + zmin(i)-1;
        t.Vol2PunctaIdx = t.Vol2PunctaIdx(idx);
        allSyd1filteredData = [allSyd1filteredData; t.filteredData];
        allSyd1Vol2PunctaIdx = logical([allSyd1Vol2PunctaIdx; t.Vol2PunctaIdx]);
    end
    
    parfor_progress;
end
parfor_progress(0);

%% get mask coordinates
mz = 5186;
% rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/Mask';
rt = 'X:\4Stephan\171120_PNnc82Syd1\Mask';
p = recursiveDir(rt);
allCoord = cell(1,numel(p)-1);

% rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/Mask/';
rt = 'X:\4Stephan\171120_PNnc82Syd1\Mask\';
fn = dir([rt '*.mat']);
fn = {fn.name}';
clear coord
for k =12
    coord(1:mz) = struct('x',[], 'y', [], 'z', [], 'label', []);
    %     FE = zeros(1, mz);
    load([rt fn{k-8}], 'FE');
    parfor_progress(mz);
    parfor j = 1:mz
        if exist([p{k} num2str(j) '.nrrd.tif'],'file') && FE(j) == 1
            [mask] = readtiff([p{k} num2str(j) '.nrrd.tif']);
            u = unique(mask(mask>0));
            if ~isempty(u)
                for i = 1:numel(u)
                    [y,x,z] = ind2sub(size(mask), find(mask == u(i)));
                    coord(j).x{i} = x;
                    coord(j).y{i} = y;
                    coord(j).z{i} = repmat(j, numel(x),1);
                    coord(j).label{i} = repmat(u(i),numel(x),1);
                end
            end
            FE(j) = 1;
        end
        parfor_progress;
    end
    parfor_progress(0);
    save([p{k} 'data.mat'], 'coord', 'FE');
    allCoord{k} = coord;
end

%%
rt = 'X:\4Stephan\171120_PNnc82Syd1\Mask\';
fn = dir([rt '*.mat']);
fn = {fn.name}';
mz = 5186;

load(['X:\4Stephan\171120_PNnc82Syd1\allnc82_PNIdx.mat']);
load(['X:\4Stephan\171120_PNnc82Syd1\allSyd1_PNIdx.mat']);

PNnc82XYZ = allnc82filteredData(allnc82Vol2PunctaIdx,:);
PNSyd1XYZ = allSyd1filteredData(allSyd1Vol2PunctaIdx,:);

for i = 1:numel(fn)
    load([rt fn{i}]);
    l = arrayfun(@(x) x.label, coord, 'unif',0);
    l = horzcat(l{:});
    ul{i} = unique([vertcat(l{:})]); % unique labels
    
    allX = arrayfun(@(x) x.x, coord, 'unif',0);
    allX = horzcat(allX{:})';
    allX = vertcat(allX{:});
    
    allY = arrayfun(@(x) x.y, coord, 'unif',0);
    allY = horzcat(allY{:})';
    allY = vertcat(allY{:});
    
    allZ = arrayfun(@(x) x.z, coord, 'unif',0);
    allZ = horzcat(allZ{:})';
    allZ = vertcat(allZ{:});
    
    allXYZ = horzcat(allX, allY, allZ);
    alllabel = arrayfun(@(x) x.label, coord, 'unif',0);
    alllabel = horzcat(alllabel{:})';
    alllabel = vertcat(alllabel{:});
    uniqlabel = ul{1};
    
    [PNnc82_AssIdx, idx_nc82] = ismember(PNnc82XYZ, allXYZ, 'rows');
    [PNsyd1_AssIdx, idx_syd1] = ismember(PNSyd1XYZ, allXYZ, 'rows');
    PNnc82_AssLabelIdx = alllabel(idx_nc82(idx_nc82>0));
    PNsyd1_AssLabelIdx = alllabel(idx_syd1(idx_syd1>0));
    [~, fn2, ext] = fileparts(fn{i});
    save([rt fn2 '_mergedIdx.mat'], 'allXYZ', 'alllabel', 'uniqlabel','PNnc82_AssIdx','PNsyd1_AssIdx','PNnc82_AssLabelIdx','PNsyd1_AssLabelIdx');
    i
end
%%
mz = 5186;
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/Mask';
% rt = 'X:\4Stephan\171120_PNnc82Syd1\Mask';
p = recursiveDir(rt);
allCoord = cell(1,numel(p)-1);

for i = [11,13]
    submit = ['bsub -n16 -R"affinity[core(1)]"  -J' ' "PN' num2str(i) '" -o ~/logs/PN' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_getPNMaskCoord ' p{i} ''''];
    [~, ~] = system(submit,'-echo');
end

%%
rt = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\';
fn = dir([rt '*.mat']);
fn = {fn.name}';
mz = 5186;

load(['X:\4Stephan\171120_PNnc82Syd1\allnc82_PNIdx.mat']);
load(['X:\4Stephan\171120_PNnc82Syd1\allSyd1_PNIdx.mat']);

PNnc82XYZ = allnc82filteredData(allnc82Vol2PunctaIdx,:);
PNSyd1XYZ = allSyd1filteredData(allSyd1Vol2PunctaIdx,:);


for i = 1:numel(fn)
    load([rt fn{i}]);
    for n = 1:numel(uniqlabel)
        uniqlabel(n)
        sum(PNnc82_AssLabelIdx== uniqlabel(n))
        sum(PNsyd1_AssLabelIdx== uniqlabel(n))
        pause
    end
    
end

%% mask region PN surface area

% get all xyz coord from all masks
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/MaskDerivedCalculations/';
% rt = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\';
fn = dir([rt '*.mat']);
fn = {fn.name}';
mergedAllXYZ = [];
parfor_progress(numel(fn));
for i = 1:numel(fn)
    load([rt fn{i}], 'allXYZ');
    if i == 1
        mergedAllXYZ = [allXYZ];
    end
    mergedAllXYZ = [mergedAllXYZ; allXYZ(~ismember(allXYZ, mergedAllXYZ, 'rows'),:)];
    parfor_progress;
end
parfor_progress(0);

% get mix max for bounding box
xrange = [min(mergedAllXYZ(:,1)), max(mergedAllXYZ(:,1))];
yrange = [min(mergedAllXYZ(:,2)), max(mergedAllXYZ(:,2))];
zrange = [min(mergedAllXYZ(:,3)), max(mergedAllXYZ(:,3))];

%% load cropped image stack
load('/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/dataBoundingBox.mat')

% crop PN image stack if necessary
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/improved-slice-tiff-chroma-cleaned/ch0';
dataStack = zeros(numel(yrange(1):yrange(2)), numel(xrange(1):xrange(2)),numel(zrange(1):zrange(2)), 'uint16');

z = zrange(1):zrange(2);
parfor_progress(numel(z));
parfor i = 1:numel(z)
    im = readtiff([rt filesep num2str(z(i)) '.tif']);
    dataStack(:,:,i) = im(yrange(1):yrange(2), xrange(1):xrange(2));
    parfor_progress;
end
parfor_progress(0);

rtsave = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/improved-slice-tiff-chroma-cleaned.tif';
% save 3d cropped stack
writetiff(dataStack, rtsave);

% load file if processing is already done above
tic
dataStack = readtiff(rtsave);
toc
% crop data at mask coordinates
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/MaskDerivedCalculations/';
rtcrop = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/MaskDerivedCalculations/croppedMasks/';
mkdir(rtcrop);
% rt = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\';
fn = dir([rt '*.mat']);
fn = {fn.name}';

% input params
px = 0.097; % pixel size in um
zAniso = 0.18/px;
T = 756; % PN mem singal OTSU
[ny, nx, nz] = size(dataStack);
HSIZE = 15;
vol = cell(1, numel(fn));
SA = cell(1, numel(fn));

parfor_progress(numel(fn));
for i = 7%1:numel(fn)
    mkdir([rtcrop num2str(i) filesep]);
    load([rt fn{i}], 'allXYZ', 'alllabel');
    uniqlabel = unique(alllabel);
    % reset coordinates based on the min and max
    allXYZ(:,1) = allXYZ(:,1) - xrange(1)+1;
    allXYZ(:,2) = allXYZ(:,2) - yrange(1)+1;
    allXYZ(:,3) = allXYZ(:,3) - zrange(1)+1;
    
    Volume = zeros(1,numel(uniqlabel));
    SurfaceArea = zeros(1,numel(uniqlabel));
    for n = 1:numel(uniqlabel)
        xi = allXYZ(alllabel==uniqlabel(n),1);
        yi = allXYZ(alllabel==uniqlabel(n),2);
        zi = allXYZ(alllabel==uniqlabel(n),3);
        
        b = 0; % boundary pixels for cropping
        xa = max(min(xi)-b,1):min(max(xi)+b,nx);
        ya = max(min(yi)-b,1):min(max(yi)+b,ny);
        za = max(min(zi)-b,1):min(max(zi)+b,nz);
        cropData = dataStack(ya,xa,za);
        cleanMask = zeros(size(cropData), 'logical');
        idx = sub2ind(size(cropData), yi-min(yi)+1,xi-min(xi)+1,zi-min(zi)+1);
        cleanMask(idx) = 1;
        mask =   cropData > T;
        mask(~cleanMask) = 0; % remove data in the bounding box outside the mask
        
        % fill gaps and holes
        mask = imclose(mask, binarySphere(HSIZE));
        %         for nclose = 1:2
        %             for z = 1:size(mask,3)
        %                 mask(:,:,z) = imfill(mask(:,:,z), 8, 'holes');
        %             end
        %             for x = 1:size(mask,2)
        %                 mask(:,x,:) = imfill(squeeze(mask(:,x,:)), 8, 'holes');
        %             end
        %             for y = 1:size(mask,1)
        %                 mask(y,:,:) = imfill(squeeze(mask(y,:,:)), 8, 'holes');
        %             end
        %             for z = 1:size(mask,3)
        %                 mask(:,:,z) = imfill(mask(:,:,z), 8, 'holes');
        %             end
        %             for x = 1:size(mask,2)
        %                 mask(:,x,:) = imfill(squeeze(mask(:,x,:)), 8, 'holes');
        %             end
        %             for y = 1:size(mask,1)
        %                 mask(y,:,:) = imfill(squeeze(mask(y,:,:)), 8, 'holes');
        %             end
        %         end
        % interpolate
        [mny,mnx,mnz] = size(mask);
        [y,x,z] = ndgrid(1:mny,1:mnx,1:mnz);
        [Y,X,Z] = ndgrid(1:mny,1:mnx,1:1/zAniso:mnz);
        imask = logical(interp3(x,y,z,double(mask),X,Y,Z,'nearest'));
        clear y x z Y X Z
        % calculate surface area and vol
        perim = bwperim(imask);
        Volume(n) = sum(imask(:)) * px^3;
        SurfaceArea(n) =sum(perim(:)) * px^2;
        
        % save cropped regions
        writetiff(uint16(cropData), [rtcrop num2str(i) filesep fn{i}(1:end-4) '_cropData_label_' num2str(uniqlabel(n)) '.tif'])
        writetiff(uint8(imask), [rtcrop num2str(i) filesep fn{i}(1:end-4) '_intrpMask_label_' num2str(uniqlabel(n)) '.tif'])
        
    end
    
    vol{i} = Volume;
    SA{i} = SurfaceArea;
    save([rtcrop num2str(i) filesep fn{i}(1:end-4) '_Vol_SA_' '.mat'], 'Volume', 'SurfaceArea', 'uniqlabel');
    parfor_progress;
end
parfor_progress(0);
%% get surface area and vol based on the mask only
% crop data at mask coordinates
rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/MaskDerivedCalculations/';
rtcrop = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/MaskDerivedCalculations/croppedMasks/';
mkdir(rtcrop);
% rt = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\';
fn = dir([rt '*.mat']);
fn = {fn.name}';

% input params
px = 0.097; % pixel size in um
zAniso = 0.18/px;
T = 756; % PN mem singal OTSU
[ny, nx, nz] = size(dataStack);
HSIZE = 15;
vol = cell(1, numel(fn));
SA = cell(1, numel(fn));

parfor_progress(3);
for i = 4:6%1:numel(fn)
    mkdir([rtcrop num2str(i) filesep]);
    load([rt fn{i}], 'allXYZ', 'alllabel');
    uniqlabel = unique(alllabel);
    % reset coordinates based on the min and max
    allXYZ(:,1) = allXYZ(:,1) - xrange(1)+1;
    allXYZ(:,2) = allXYZ(:,2) - yrange(1)+1;
    allXYZ(:,3) = allXYZ(:,3) - zrange(1)+1;
    nSec = 32;
    Volume = zeros(1,nSec);
    SurfaceArea = zeros(1,nSec);
    n = 1
    xi = allXYZ(alllabel==uniqlabel(n),1);
    yi = allXYZ(alllabel==uniqlabel(n),2);
    zi = allXYZ(alllabel==uniqlabel(n),3);
    
    mask = zeros(size(dataStack), 'logical');
    idx = sub2ind(size(dataStack), yi,xi,zi);
    mask(idx) = 1;
    
    for sv = 1:nSec
        if sv == 1
            submask = mask(:,1:251, :);
        elseif sv < nSec
            submask = mask(:,(250*(sv-1)):(250*(sv)+1), :);
        else
            submask = mask(:,(250*(sv-1)+1):end, :);
        end
        
        idx = find(submask);
        [ny,nx,nz] = size(submask);
        [yi,xi,zi] = ind2sub([ny,nx,nz], idx);
        b = 5; % boundary pixels for cropping
        xa = max(min(xi)-b,1):min(max(xi)+b,nx);
        ya = max(min(yi)-b,1):min(max(yi)+b,ny);
        za = max(min(zi)-b,1):min(max(zi)+b,nz);
        submask = submask(ya,xa,za);
        if ~isempty(submask)
            % interpolate
            [mny,mnx,mnz] = size(submask);
            [y,x,z] = ndgrid(1:mny,1:mnx,1:mnz);
            [Y,X,Z] = ndgrid(1:mny,1:mnx,1:1/zAniso:mnz);
            imask = logical(interp3(x,y,z,double(submask),X,Y,Z,'nearest'));
            clear y x z Y X Z
            % calculate surface area and vol
            perim = bwperim(imask);
            cutperim =  perim(:, 2:end-1,:);
            Volume(sv) = sum(imask(:)) * px^3;
            SurfaceArea(sv) =sum(cutperim(:)) * px^2;
        end
        sv
    end
    % save cropped regions
    %         writetiff(uint16(cropData), [rtcrop num2str(i) filesep fn{i}(1:end-4) '_cropData_label_' num2str(uniqlabel(n)) '.tif'])
    %         writetiff(uint8(imask), [rtcrop num2str(i) filesep fn{i}(1:end-4) '_intrpMask_label_' num2str(uniqlabel(n)) '.tif'])
    VolParts = Volume;
    SAParts = SurfaceArea;
    Volume = sum(Volume(:));
    SurfaceArea = sum(SurfaceArea(:));
    
    
    vol{i} = Volume;
    SA{i} = SurfaceArea;
    save([rtcrop num2str(i) filesep fn{i}(1:end-4) '_Vol_SA_' '.mat'], 'Volume', 'SurfaceArea', 'VolParts', 'SAParts', 'uniqlabel');
    parfor_progress;
end
parfor_progress(0);

%%
rt = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\';
rt_vol_sa = 'X:\4Stephan\171120_PNnc82Syd1\MaskDerivedCalculations\croppedMasks\';
fn = dir([rt '*.mat']);
fn = {fn.name}';
EF = 3.99;
for i = 1:numel(fn)
    
    load([rt fn{i}]);
    load([rt_vol_sa num2str(i) filesep fn{i}(1:end-4) '_Vol_SA_.mat']);
    
    for n = 1:numel(uniqlabel)
        PNnc82_n{i}(n) = sum(PNnc82_AssLabelIdx == uniqlabel(n));
        PNnc82_n_SA{i}(n) = PNnc82_n{i}(n)/SurfaceArea(n)/EF^2;
        PNnc82_n_Vol{i}(n) = PNnc82_n{i}(n)/Volume(n)/EF^3;
        
        PNSyd1_n{i}(n) = sum(PNsyd1_AssLabelIdx == uniqlabel(n));
        PNSyd1_n_SA{i}(n) = PNSyd1_n{i}(n)/SurfaceArea(n)/EF^2;
        PNSyd1_n_Vol{i}(n) = PNSyd1_n{i}(n)/Volume(n)/EF^3;
    end
    PNnc82_n{i} = PNnc82_n{i}';
    PNnc82_n_SA{i} = PNnc82_n_SA{i}';
    PNnc82_n_Vol{i} = PNnc82_n_Vol{i}';
    
    PNSyd1_n{i} = PNSyd1_n{i}';
    PNSyd1_n_SA{i} = PNSyd1_n_SA{i}';
    PNSyd1_n_Vol{i} = PNSyd1_n_Vol{i}';
end

%% dendrite data test
sigVox = GU_calcGaussianIntegral(1, [2.5,2.5,2.5]);

fnv1 = data.framePaths{2}{1}; % nc82
fnv2 = data.framePaths{3}{1}; % syd1
fnv3 = data.framePaths{1}{1}; % PN


% GU_ExM_zOffset3Dframe(fnv1, 1,1);
% GU_ExM_zOffset3Dframe(fnv2, 1,6);
memT = 1960;


GU_ExM_calcPunctaDensityVol(fnv1, fnv3,...
    'FindLocalMaxima', true, 'MinThreshold', [576,memT],'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [round(sigVox), 1008],'ExpansionFactor', 3.99,'MinVolumeOccupancy', 0.0);
% Syd1 + PN
GU_ExM_calcPunctaDensityVol(fnv2, fnv3,...
    'FindLocalMaxima', true, 'MinThreshold', [119,memT],'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [round(sigVox), 1008],'ExpansionFactor', 3.99,'MinVolumeOccupancy', 0.0);

%%
nc82 = load('D:\Gokul\Wholebrainsynapse\Ex22_PN+nc82+Syd1_testneardendrites_z0p18\ch1\Analysis_1px756\246\DataStructures\channel 0 [1270, 4008, 4792]_246_densityData.mat', 'filteredData', 'Vol2PunctaIdx');
syd1 = load('D:\Gokul\Wholebrainsynapse\Ex22_PN+nc82+Syd1_testneardendrites_z0p18\ch2\Analysis_1px756\246\DataStructures\channel 0 [1270, 4008, 4792]_246_densityData.mat', 'filteredData', 'Vol2PunctaIdx');
zAniso = .18/0.097;
syd1.filteredData(:,3) = syd1.filteredData(:,3)*zAniso;
nc82.filteredData(:,3) = nc82.filteredData(:,3)*zAniso;

Mdl = KDTreeSearcher(syd1.filteredData(syd1.Vol2PunctaIdx)); % all nonDAN synapse positions
[~,D] = knnsearch(Mdl, nc82.filteredData(nc82.Vol2PunctaIdx), 'k',1,'NSMethod', 'kdtree');
figure, ecdf(D)



Mdl = KDTreeSearcher(nc82.filteredData(nc82.Vol2PunctaIdx)); % all nonDAN synapse positions
[~,D] = knnsearch(Mdl, syd1.filteredData(syd1.Vol2PunctaIdx), 'k',1);
figure, ecdf(D)

%%
for i = 1:7
    nc82{i} = [PNnc82_n{i},PNnc82_n_SA{i}, PNnc82_n_Vol{i}];
    syd1{i} = [PNSyd1_n{i},PNSyd1_n_SA{i}, PNSyd1_n_Vol{i}];
end

%% old PN calc thresh

MinHoleVolume = 1008;% for membrane
% MinHoleVolume = 516; % for puncta

for i = 1%1:numel(data)
    for j = 2%1:numel(data(i).framePaths)
        fprintf('reading image...')
        tic
        im = readtiff(data(i).framePaths{j}{1});
        toc
        
        fprintf('filtering image...')
        tic
        imG = filterGauss3D(im,1);
        T = thresholdOtsu(imG(imG > 0 & imG < prctile(im(:),99)));
        imG = imG-T;
        imG(imG<0) = 0;
        toc
        
        % discard small holes
        fprintf('discarding small holes...')
        tic
        CC = bwconncomp(logical(imG), 26);
        csize = cellfun(@numel, CC.PixelIdxList);
        idx = csize>=MinHoleVolume;
        CC.NumObjects = sum(idx);
        CC.PixelIdxList = CC.PixelIdxList(idx);
        imGL = labelmatrix(CC)~=0;
        toc
        imG = uint16(imG) .* uint16(imGL);
        [p, ~, ~] = fileparts([data(i).framePathsDS{j}{1}]);
        mkdir(p);
        
        
        fprintf('writing image...')
        tic
        writetiff(imG, [data(i).framePathsDS{j}{1}]);
        
        toc
        
        clear im imG imGL;
    end
end
% check signal intensity of the clean mask in 1px eroded region
% imG = im;
se = strel('sphere',1);
cleanDAN_E = imerode(logical(imG), se);
figure, histogram(im(imG>0 & ~cleanDAN_E ),1000);
muCut = prctile(im(imG>0 & ~cleanDAN_E),99.9);

% PN mem T = 441 (346+1px)
% PN nc82 T = 435 (307+1px)
%% old PN

sigVox(1) = GU_calcGaussianIntegral(1, [2.5,2.5,2.5]);
% fnv1 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch1\channel 1 [8775, 4298, 4954].tif']; % nc82
% fnv2 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch2\channel 2 [8775, 4298, 4954].tif']; % syd1
% fnv3 = ['D:\Gokul\Wholebrainsynapse\Ex16_PN+nc82+Syd1_test_z0p18\ch0\channel 0 [8775, 4298, 4954].tif']; % PN

for i = 1:numel(data)
    % fnv1 = [data(i).framePaths{1}{1}]; % nc82
    % fnv2 = [data(i).framePaths{2}{1}]; %  PN
    
    %     GU_ExM_zOffset3Dframe(fnv1, 1,1);
    %     GU_ExM_zOffset3Dframe(fnv2, 1,6);
    
    
    % depuncate and calculate the local maxima
    for ni = 1:numel(sigVox)
        % nc82 + PN
        GU_ExM_calcPunctaDensityVol(fnv1, fnv2,...
            'FindLocalMaxima', true, 'MinThreshold', [435,441],'Verbose', true,'OTSUGaussKernel', 0.5,...
            'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99,'minAZdia', 0.1);
    end
end

%% old PN mask read out
k = 1;
d = load([data(k).source 'Analysis' filesep '246' filesep 'DataStructures' filesep '2membrane_246_densityData.mat'], 'filteredData');
filteredData = d.filteredData;
mask{1} = '/groups/betzig/betziglab/4Stephan/171128_oldPNnc82/Ex1_AL1_z0p18/Mask/3mask.tif';

%% PN stereotypy analysis


% get OTSU threshold values
for i = 1:numel(data)
    for j = 1:numel(data(i).framePaths{1})
        fprintf('reading image...')
        tic
        im = readtiff(data(i).framePaths{1}{j});
        toc
        
        fprintf('filtering image...')
        tic
        T{i}(j) = thresholdOtsu(im(im > 0 & im < prctile(im(:),99)));
        toc
    end
end






% load the data from:
% rt = '/groups/betzig/betziglab/4Stephan/171211_PN_stereo/';
for i = 1:numel(data)
    submit = ['bsub -n 32 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "PNL' num2str(i) '" -o ~/logs/PNL' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_PunctaDensity_LARGEVOL ' data(i).framePaths{1}{1} ' ' data(i).framePaths{2}{1} ' ' num2str(T{i}(2)) ' ' num2str(T{i}(1)) ''''];
    [stat, res] = system(submit,'-echo');
end


%% stereotypy
rt = '/groups/betzig/betziglab/4Stephan/171211_PN_stereo/';
cd(rt)
load('loadCond_IS.mat')
ef = 3.95;
se = strel('disk',3);
zAniso = 0.18/0.097;
px = 0.097/ef;
for k = 1:numel(data)
    
    % read membrane channel
%     imrt = dir([data(k).source 'Analysis' filesep '1008' filesep '*.tif']);
%     tic
%     im = readtiff([data(k).source 'Analysis' filesep '1008' filesep imrt.name]);
%     toc
    % read mask
    m = dir([data(k).source 'Analysis' filesep '1008' filesep 'Boutons' filesep 'Masks' filesep '*.tif']);
    % loop over each cell
    
    for j = 1:numel(m)
        savert  = [m(j).folder filesep m(j).name(1:end-4) filesep];
        mkdir(savert);
        cd (savert)
        if ~exist([savert 'BoutonSurfaceArea_Vol.mat'],'file')
            
%             % read mask
%             mask =  logical(readtiff([m(j).folder filesep m(j).name]));
%             
%             % break mask into 2 parts and upsample to match image
%             [mny,mnx,mnz] = size(mask);
%             mask1 = mask(:,:,1:mnz/2);
%             mask2 = mask(:,:,mnz/2+1:end);
%             [y,x,z] = ndgrid(1:mny,1:mnx,1:mnz/2);
%             [Y,X,Z] = ndgrid(1:1/2.001:mny,1:1/2.001:mnx,1:1/2.003:mnz/2);
%             imask1 = logical(interp3(x,y,z,mask1,X,Y,Z,'nearest'));
%             imask2 = logical(interp3(x,y,z,mask2,X,Y,Z,'nearest'));
%             clear y x z Y X Z
%             imask = cat(3, imask1, imask2);
%             clear  imask1 imask2 mask1  mask2 mask
%             imask = imdilate(imask,se);
%             maskedIm = zeros(size(im), 'uint16');
%             maskedIm(imask) = im(imask);
%             [l,n] = bwlabeln(logical(maskedIm));
%             writetiff(uint8(l), [savert 'labeledBoutons.tif']);
%             [y,x,z] = ind2sub(size(maskedIm),find(maskedIm));
%             Coord_xyz_id = [x,y,z, l(maskedIm>0)];
%             save([savert 'BoutonCoord_id.mat'], 'Coord_xyz_id');
%             clear x y z Coord_xyz_id
                n = dir([savert 'CroppedBouton_*.tif']);
                n = numel(n);
            % crop each bouton, calculate surface area and volume
            for nl = 1:n
%                 if ~exist([savert 'CroppedBouton_' num2str(nl) '.tif'], 'file')
%                 idx = find(l == nl);
%                 [ny,nx,nz] = size(maskedIm);
%                 [yi,xi,zi] = ind2sub([ny,nx,nz], idx);
%                 b = 5; % boundary pixels for cropping
%                 xa = max(min(xi)-b,1):min(max(xi)+b,nx);
%                 ya = max(min(yi)-b,1):min(max(yi)+b,ny);
%                 za = max(min(zi)-b,1):min(max(zi)+b,nz);
%                 cm = logical(maskedIm(ya,xa,za));
%                 % close gaps
%                 cm = imclose(cm, binarySphere(3));
%                 for z = 1:size(cm,3)
%                     cm(:,:,z) = imfill(cm(:,:,z), 8, 'holes');
%                 end
%                 for x = 1:size(cm,2)
%                     cm(:,x,:) = imfill(squeeze(cm(:,x,:)), 8, 'holes');
%                 end
%                 for y = 1:size(cm,1)
%                     cm(y,:,:) = imfill(squeeze(cm(y,:,:)), 8, 'holes');
%                 end
%                 writetiff(uint8(cm), [savert 'CroppedBouton_' num2str(nl) '.tif']);
%                 else
                    cm = readtiff([savert 'CroppedBouton_' num2str(nl) '.tif']);
%                 end
                % interpolate
                [mny,mnx,mnz] = size(cm);
                [y,x,z] = ndgrid(1:mny,1:mnx,1:mnz);
                [Y,X,Z] = ndgrid(1:mny,1:mnx,1:1/zAniso:mnz);
                icm = logical(interp3(x,y,z,cm,X,Y,Z,'nearest'));
                perim = bwperim(icm);
                vol(nl) = sum(icm(:)) * px^3;
                area(nl) = sum(perim(:)) * px^2;
                clear icm cm y x z Y X Z
                fprintf(['Completed bouton %d of %d \n'],nl, n);
            end
            save([savert 'BoutonSurfaceArea_Vol.mat'], 'vol', 'area');
            clear vol area
            fprintf('Completed cell %d of %d \n',j, numel(m));
        end
    end
    fprintf('Completed PN %d of %d \n',k, numel(data));
    
end

%%  figures for the stereotypy

rt = '/groups/betzig/betziglab/4Stephan/171211_PN_stereo/';
cd(rt)
load('loadCond_IS.mat')
figure
hold on
for k = 1:numel(data)
    m = dir([data(k).source 'Analysis' filesep '1008' filesep 'Boutons' filesep 'Masks' filesep '*.tif']);
    
   x = k;
    for j = 1:numel(m)
        statsrt  = [m(j).folder filesep m(j).name(1:end-4) filesep];
        PN(k).cell(j) = load([statsrt 'BoutonSurfaceArea_Vol.mat']);
        plot(x, PN(k).cell(j).area, '.', 'MarkerSize',30);
        hold on
        x = x+0.2;
    end
    
end

%% stereotypy - number of nc82 in every region
rt = '/groups/betzig/betziglab/4Stephan/171211_PN_stereo/';
cd(rt)
load('loadCond_IS.mat')
for k = 1:numel(data)
    m = dir([data(k).source 'Analysis' filesep '1008' filesep 'Boutons' filesep 'Masks' filesep '*.tif']);
    nc82rt = [data(k).channels{2} 'Analysis' filesep '246' filesep 'DataStructures' filesep];
    load([nc82rt 'ch0_246_densityData.mat'], 'Vol2PunctaIdx', 'filteredData');
   
    for j = 1:numel(m)
        statsrt  = [m(j).folder filesep m(j).name(1:end-4) filesep];
        load([statsrt 'BoutonCoord_id.mat']);
        coord = Coord_xyz_id(:,1:3);
        label = Coord_xyz_id(:,4);
        clear Coord_xyz_id
        nc82coord = filteredData(Vol2PunctaIdx,:);
        idx = ismember(filteredData, coord, 'rows');
        PN(k).nc82(j) = load([statsrt 'BoutonCoord_id.mat']);
        plot(x, PN(k).cell(j).area, '.', 'MarkerSize',30);
        hold on
        x = x+0.2;
    end
    
end
%% determine the distance of nc82 and syd1 in the masked regions
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/PN_nc82Syd1/';
cd(rt)
load('allnc82_PNIdx.mat')
load('allSyd1_PNIdx.mat')

ef = 3.95;
zAniso = 0.18/0.097;
px = 0.097/ef;

nc82 = allnc82filteredData(allnc82Vol2PunctaIdx,:);
Syd1 = allSyd1filteredData(allSyd1Vol2PunctaIdx,:);

nc82(:,3) = nc82(:,3) * zAniso;
Syd1(:,3) = Syd1(:,3) * zAniso;

fl = dir('*_PN-*.mat');
fl = {fl.name};

allMasked_nc82 = cell(1, numel(fl));
allMasked_syd1 = cell(1, numel(fl));
allDist_zcorr = cell(1, numel(fl));

ha = setupFigure(4,2, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

for i = 1:numel(fl)
    load(fl{i}, 'PNnc82_AssIdx', 'PNsyd1_AssIdx')
    allMasked_nc82{i} = nc82(PNnc82_AssIdx,:);
    allMasked_syd1{i} = Syd1(PNsyd1_AssIdx,:);
    
    tic
    Mdl = KDTreeSearcher(allMasked_syd1{i}); % 
    toc
    tic
    % IdxKDT = rangesearch(Mdl,allMasked_syd1{i},bandwidth2);% Search radius = bandwidth2;
    [~,D] = knnsearch(Mdl, allMasked_nc82{i}, 'k',1); %distance of the closest 
    toc
    allDist_zcorr{i} = D*px*1000; % in nm
    
    axes(ha(i))
    [f,x] = ecdf(allDist_zcorr{i});
    plot(x,f)
    
    % add text labels
    for j = [100,150,200]
        idx = find(abs(x-j)<10 & abs(x-j)==min(abs(x-j)));
        x0 = x(idx);
        tl = f(idx);
        text(x(idx),tl, ['x=' num2str(x(idx)) ',' num2str(tl)]);
        line([x(idx), x(idx)], [0 tl]);
    end
    title(fl{i});
    xlabel('nc82 distance to closest Syd1 (nm)');
    ylabel('Cumulative Frequency');
    xlim([0 500]);
    ylim([0 1]);
end

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' IndividualRegions_closestnc82_syd1.eps']);
saveas(f0,[date ' IndividualRegions_closestnc82_syd1.fig'], 'fig')



ha = setupFigure(3,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

r = {1:3, 4:6,7};
l = {'Boutons', 'Axons', 'Dendrites'};
for i = 1:numel(r)
    Masked_nc82 = vertcat(allMasked_nc82{r{i}});
    Masked_syd1 = vertcat(allMasked_syd1{r{i}});
    Dist_zcorr = vertcat(allDist_zcorr{r{i}}); % in nm
    
    axes(ha(i))
    [f,x] = ecdf(Dist_zcorr);
    plot(x,f)
    
    % add text labels
    for j = [100,150,200]
        idx = find(abs(x-j)<10 & abs(x-j)==min(abs(x-j)));
        x0 = x(idx);
        tl = f(idx);
        text(x(idx),tl, ['x=' num2str(x(idx)) ',' num2str(tl)]);
        line([x(idx), x(idx)], [0 tl]);
    end
    title(l{i});
    xlabel('nc82 distance to closest Syd1 (nm)');
    ylabel('Cumulative Frequency');
    xlim([0 500]);
    ylim([0 1]);
end

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'Grouped_closestnc82_syd1.eps']);
saveas(f0,[date 'Grouped_closestnc82_syd1.fig'], 'fig')


[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(3,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

r = {1:3, 4:6,7};
l = {'Boutons', 'Axons', 'Dendrites'};
for i = 1:numel(r)
    Masked_nc82 = vertcat(allMasked_nc82{r{i}});
    Masked_syd1 = vertcat(allMasked_syd1{r{i}});
    Dist_zcorr = vertcat(allDist_zcorr{r{i}}); % in nm
    
    axes(ha(i))
    histogram(Dist_zcorr, 150, 'EdgeColor', cfK, 'FaceColor', cfK,'BinWidth',5)
  legend(['n = ' num2str(numel(Dist_zcorr))])
    legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica')
    title(l{i});
    xlabel('nc82 distance to closest Syd1 (nm)');
    ylabel('Cumulative Frequency');
    xlim([0 500]);
%     ylim([0 1]);
end

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'Grouped_Histogram_closestnc82_syd1.eps']);
saveas(f0,[date 'Grouped_Histogram_closestnc82_syd1.fig'], 'fig')
