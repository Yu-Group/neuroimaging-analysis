rt = 'I:\Gokul\GDrive\JaneliaSharedwithGokul\ExpansionRelated\Nc82_synapseData';
load([rt filesep 'dataLoad.mat']);
sigma = [1.26, 1.4; 1.4, 1.5];
overwrite = true;
apath =[];
Mode = 'xyzAc';
MinHoleVolume = 50;
FitMixtures = false;
MaxMixtures = 1;
if isempty(apath) % store results locally in cell/Analysis directory
    apath = arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0);
else
    % expand results path for all data sets
    apath = arrayfun(@(i) [apath getShortPath(i,3) filesep 'Analysis' filesep], data, 'unif', 0);
end
[~,~] = cellfun(@mkdir, apath, 'unif', 0);

AZdia = .2;% diameter of activezone
clusterSizeMin = 0;
sigmas = [1.3, 1.3];
ef = 4;
px = 0.097; %in um
bandwidth = AZdia/2/.097*ef;
V = 1;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % expanded radius in um
bandwidth2 = re/px;
strctmask{1} = 'I:\Gokul\GDrive\JaneliaSharedwithGokul\ExpansionRelated\Nc82_synapseData\Ex1_Mushroombodyalpha3_z0p18\Mushroombodyalpha3_a3-mask.tif';
strctmask{2} = 'I:\Gokul\GDrive\JaneliaSharedwithGokul\ExpansionRelated\Nc82_synapseData\Ex2_Calyxpart1_z0p18\nc82-a3-BrpV5 samples_calyx-mask.tif';
strctmask{3} = 'I:\Gokul\GDrive\JaneliaSharedwithGokul\ExpansionRelated\Nc82_synapseData\Ex3_Mushrioombodyalpha3_C1_z0p18\a3-mask.tif';

%% clean up and generate local max map

for i = 3%1:numel(data)
    for j = 1:numel(data(i).framePaths)
        fprintf('reading image...')
        tic
        im = readtiff(data(i).framePaths{j}{1});
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
        
        fprintf('calculating Local Maxima image...')
        tic;
        allMax = locmax3d(imG, 2*ceil(sigma([1 1 2]))+1, 'ClearBorder', false);
        toc
        
        fprintf('writing image...')
        tic
        writetiff(imG, [data(i).framePathsDS{j}{1}]);
        writetiff(uint8(allMax), [data(i).source 'Analysis' filesep 'LocalMax.tif']);
        toc
        
        clear im imG imGL;
    end
end
%% test detection section
fprintf('running Detection3D on image...')
tic
GU_runDetection3D(data, 'Sigma', sigma, 'Overwrite', overwrite,...
    'WindowSize', [], 'ResultsPath', apath,'Mode', Mode,...
    'FitMixtures',FitMixtures, 'MaxMixtures',MaxMixtures);
toc

% opts = {'GaussSigma', 3, 'PixelSize', 0.097, 'PixelSizeZ', 0.18,'LargestObject', false,'DistanceTransform', false,'MinBk', 300,...
%         'MaxPercentileCutoff',100,'Fill2D', false};
%     [filledDM, unfilledSeg, T] = GU_segmentCytosolicVol3D(im, opts{:});
%

[pstruct, mask] = pointSourceDetection3D(im, sigma(1,:), 'Alpha', 0.05,...
    'Mask', [], 'RemoveRedundant', true,...
    'RefineMaskLoG', false, 'WindowSize', [], 'Mode', Mode); %#ok<PFBNS>


%% cluster objects
se = strel('sphere', 2);
for i = 1:numel(data)
    mkdir([data(i).source 'Analysis'])
    clear clus clusterInfo ptData s
    fprintf('reading image...')
    tic
    im = logical(readtiff([data(i).source 'Analysis' filesep 'LocalMax.tif']));
    outlineMask = logical(readtiff(strctmask{i}));
    toc
    
    im(outlineMask == 0) = 0;
    fprintf('calculating centroids...')
    tic
    %     imMask = logical(im);
    
    s  = regionprops(im,'centroid');
    imcentroids = cat(1, s.Centroid);
    imcX = imcentroids(:,1);
    imcY = imcentroids(:,2);
    imcZ = imcentroids(:,3);
    toc
    
    % first round of clustering to remove duplicate detections
    fprintf('clustering to remove duplicates (diameter pre expanded ~ 200)*expansion factor...')
    tic
    ptData(:,1) = imcX;
    ptData(:,2) = imcY;
    ptData(:,3) = imcZ;
    [clusterInfo, ~, ~] = MeanShiftClustering( ptData, bandwidth);
    toc
    
    % second round of clustering to group
    fprintf('clustering to calculate spot density (diameter -- 1um3)*expansion factor^3...')
    tic
    ClusCenters = cat(1, clusterInfo.ptClusterCenter);
    [clusterInfo2, ~, ~] = MeanShiftClustering( ClusCenters, bandwidth2);
    toc
    
    ClusCenters2 = cat(1, clusterInfo2.ptClusterCenter);
    cm = jet(numel([clusterInfo2.numPoints]));
    cm = flip(cm,2);
    [~,idx] = sort([clusterInfo2.numPoints], 'descend');
    
    % make figures
    fprintf('generating figures...')
    tic
    clear clus x2 y2 z2
    ha = setupFigure(1,1, 'AxesWidth', 10, 'AxesHeight', 10,'SameAxes', false,...
        'XSpace', [1.5 2 1], 'YSpace', [1.25 1 1]);
    set(ha, 'FontSize', 8);
    axes(ha(1)) % cell vol
    %     scatter3(ClusCenters(:,1)*px,ClusCenters(:,2)*px,ClusCenters(:,3)*px,10, 'filled'),hold on;
    scatter3(ClusCenters2(idx,1)*px,ClusCenters2(idx,2)*px,ClusCenters2(idx,3)*px,5*([clusterInfo2.numPoints]+1), cm, 'filled','MarkerFaceAlpha',0.3);
    hold on
    axis equal
    view(3)
    xlabel('x-axis (um^3)')
    ylabel('y-axis (um^3)')
    zlabel('z-axis (um^3)')
    
    for n = 1:numel(idx)
        x2= ClusCenters(clusterInfo2(idx(n)).ptIdData,1);
        y2= ClusCenters(clusterInfo2(idx(n)).ptIdData,2);
        z2= ClusCenters(clusterInfo2(idx(n)).ptIdData,3);
        scatter3(x2*px, y2*px, z2*px, 1, cm(n,:), 'filled'); hold on
    end
    f0 = gcf();
    print(f0, '-dsvg', [data(i).source 'Analysis' filesep 'nc82spotDensity_balls_surveypreexpandedVol' num2str(V) 'um3.svg']);
    
    ha = setupFigure(1,1, 'AxesWidth', 10, 'AxesHeight', 10,'SameAxes', false,...
        'XSpace', [1.5 2 1], 'YSpace', [1.25 1 1]);
    set(ha, 'FontSize', 8);
    %     axes(ha(1))
    for n = 1:numel(idx)
        x2= ClusCenters(clusterInfo2(idx(n)).ptIdData,1);
        y2= ClusCenters(clusterInfo2(idx(n)).ptIdData,2);
        z2= ClusCenters(clusterInfo2(idx(n)).ptIdData,3);
        scatter3(x2*px, y2*px, z2*px, 10, cm(n,:), 'filled','MarkerFaceAlpha',0.3); hold on
    end
    axis equal
    view(3)
    xlabel('x-axis (um^3)')
    ylabel('y-axis (um^3)')
    zlabel('z-axis (um^3)')
    toc
    f0 = gcf();
    print(f0, '-dsvg', [data(i).source 'Analysis' filesep 'nc82spotDensity_surveypreexpandedVol' num2str(V) 'um3.svg']);
    
    fprintf('generating intenisty coded raw mask data...')
    tic
    % write image and them convolve with psf
    [lholes, NUM] = bwlabeln(imdilate(uint8(im), se));
    lholes = (lholes);
    %     CC = bwconncomp(lholes, 26);
    %     csize = cellfun(@numel, CC.PixelIdxList);
    
    mask = zeros(size(im), 'single');
    parfor_progress(numel(idx));
    for n = 1:numel(idx)
        x2= round(ClusCenters(clusterInfo2(idx(n)).ptIdData,1));
        y2= round(ClusCenters(clusterInfo2(idx(n)).ptIdData,2));
        z2= round(ClusCenters(clusterInfo2(idx(n)).ptIdData,3));
        val = clusterInfo2(n).numPoints;%*((2*pi)^(3/2)* prod(sigmas));
        lidx = zeros(numel(x2),'single');
        for nn = 1:numel(x2)
            lidx(nn) =  lholes(y2(nn),x2(nn),z2(nn));
        end
        mask(ismember(lholes, lidx)) = val;
        parfor_progress;
    end
    parfor_progress(0);
    toc
    
    maskVol(i) = sum(logical(outlineMask(:)))*0.097^2*0.18; % in um^3
    nSynapses(i) = numel(clusterInfo);
    fprintf('writing intensity coded image and saving cluster info...')
    tic
    writetiff(mask, [data(i).source 'Analysis' filesep 'DensityMap.tif']);
    save([data(i).source 'Analysis' filesep 'ClusterInfo.mat'], 'clusterInfo', 'ClusCenters', 'clusterInfo2', 'maskVol', 'nSynapses');
    toc
end
synapseDensity = nSynapses./(maskVol / ef^3);

%%
GU_ExM_calcSynapseDensityVol(data(i).framePaths{1}{2}, data(i).framePaths{1}{1},'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);


%%
MinHoleVolume = 246;
load('D:\Gokul\Wholebrainsynapse\Ex13_alpha3_left_z0p18\Analysis\DataStructures\channel 1 [3496, 12063, 11498]-1-a3-L_densityData.mat')
DAN = readtiff(data.framePaths{1}{1});
cleanDAN = readtiff('D:\Gokul\Wholebrainsynapse\Ex13_alpha3_left_z0p18\Analysis\channel 0 [3496, 12063, 11498]-1-a3-L_clean.tif');
mask = logical(readtiff('D:\Gokul\Wholebrainsynapse\Ex13_alpha3_left_z0p18\mask\a3-mask.tif'));
[XYZ(:,1),XYZ(:,2),XYZ(:,3)] = GU_getMaskCoordinates(mask);

DANSynapsesXYZ = filteredData(DANSynapseIdx,:);
SynapseinMask_idx = ismember(filteredData,XYZ, 'rows');
DANSynapseinMask_idx = ismember(DANSynapsesXYZ,XYZ, 'rows');
nDANSynapseinMask = sum(DANSynapseinMask_idx);
nDANSynPerofTotal = nDANSynapseinMask/sum(SynapseinMask_idx);


% remove small debris
 CC = bwconncomp(cleanDAN, 26);
        csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinHoleVolume;
        CC.NumObjects = sum(idx);
        CC.PixelIdxList = CC.PixelIdxList(idx);
        cleanDAN = labelmatrix(CC)~=0;
   
for i = 1:numel(filteredData(:,1))
    idx_smalldebris(i) = cleanDAN(filteredData(i,2), filteredData(i,1), filteredData(i,3));
end
        
        
% check signal intensity of the clean mask in 1px eroded region
cleanDAN_G = filterGauss3D(DAN,0.5)-T(2);
se = strel('sphere',1);
cleanDAN_E = imerode(logical(cleanDAN), se);
figure, histogram(cleanDAN_G(cleanDAN>0 & ~cleanDAN_E ),1000);
muCut = prctile(cleanDAN_G(cleanDAN>0 & ~cleanDAN_E),99.9);

nDANSynapseinMask_aboveThresh = sum(SynapseinMask_idx(DANSynapseIdx) & DANchSignal(DANSynapseIdx) >  muCut & idx_smalldebris((DANSynapseIdx))');
nDANSynPerofTotal_aboveThresh = nDANSynapseinMask_aboveThresh/sum(SynapseinMask_idx);

% % erosion test readout of DAN syn idx
% cleanDAN_Emasked = cleanDAN_E .* mask;
% for i = 1:numel(DANSynapsesXYZ(:,1))
%     idx_erroded(i) = cleanDAN_Emasked(DANSynapsesXYZ(i,2), DANSynapsesXYZ(i,1), DANSynapsesXYZ(i,3));
% end


%%
% f = DANchSignal > 22;
for i = 1:numel(DANSynapsesXYZ(:,1))
    idxD(i) = mask(DANSynapsesXYZ(i,2), DANSynapsesXYZ(i,1), DANSynapsesXYZ(i,3));
end
sum(idxD)
for i = 1:numel(filteredData(:,1))
    idxA(i) = mask(filteredData(i,2), filteredData(i,1), filteredData(i,3));
end
sum(idxA)

figure,
imagesc(mask(:,:,500)), hold on
% plot(filteredData(filteredData(idxA,3)==500,1),filteredData(filteredData(idxA,3)==500,2), 'g.')
plot(filteredData([filteredData(:,3)==500 & idxA'],1),filteredData([filteredData(:,3)==500 & idxA'],2), 'g.')
plot(DANSynapsesXYZ([DANSynapsesXYZ(:,3)==500 & idxD'],1),DANSynapsesXYZ([DANSynapsesXYZ(:,3)==500 & idxD'],2), 'ro')

figure, scatter3(filteredData(idxA,1),filteredData(idxA,2),filteredData(idxA,3),3,filteredsynapseDensity(idxA), 'filled'), axis equal; xlabel('xaxis'); ylabel('yaxis'); colormap jet

GU_stairsHistPlot({[],[],filteredsynapseDensity(DANSynapseIdx),filteredsynapseDensity(~DANSynapseIdx)},'ShowECDF', true, 'BinSize',1)


%%
tic
clc
i =1;
GU_ExM_calcSynapseDensityVol(data(i).framePaths{1}{2}, data(i).framePaths{1}{1},'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);
[p2,fn2,~] = fileparts(data(i).framePaths{1}{1});
[p1,fn1,~] = fileparts(data(i).framePaths{1}{2});
GU_ExM_correctDANSynapseDensityCalc([p1 filesep 'Analysis' filesep 'DataStructures' filesep fn1 '_densityData.mat'],[p2 filesep 'Analysis' filesep fn2 '_clean.tif']);
toc

%% load data ex13 as data(1); ex14 as data(2)
s = [0.0237, 0.0237, 0.044];
% left
leftMask = readtiff([data(1).source 'mask' filesep 'a3-mask.tif']);
rightMask = readtiff([data(2).source 'mask' filesep 'channel 1 [3956, 17001, 10986]-1-mask.tif']);
lData = load([[data(1).source 'Analysis' filesep 'Analysis' filesep 'DataStructures_Recalc' filesep 'channel 0 [3496, 12063, 11498]-1-a3-L_clean_densityDataRC.mat']]);
rData = load([[data(2).source 'Analysis' filesep 'Analysis' filesep 'DataStructures_Recalc' filesep 'channel 0 [3956, 17001, 10986]-1_clean_densityDataRC.mat']]);

x = ceil(lData.filteredData(:,1));
y = ceil(lData.filteredData(:,2));
z = ceil(lData.filteredData(:,3));

[sy, sx, sz] = size(leftMask);
[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
lData.maskedIdx = interp3(X,Y,Z,leftMask,x,y,z,'nearest');
toc

x = ceil(rData.filteredData(:,1));
y = ceil(rData.filteredData(:,2));
z = ceil(rData.filteredData(:,3));
[sy, sx, sz] = size(rightMask);
[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
rData.maskedIdx = interp3(X,Y,Z,rightMask,x,y,z,'nearest');
toc
%%
% load('DataStructures_ex13_14.mat')
load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/Wholebrainsynapse/DataStructures_ex13_14.mat')
cd('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/Ex14_alpha3_right')
%% distance - any 2 any
tic
fData(:,1) = lData.filteredData(logical(lData.maskedIdx),1)*0.0237;
fData(:,2) = lData.filteredData(logical(lData.maskedIdx),2)*0.0237;
fData(:,3) = lData.filteredData(logical(lData.maskedIdx),3)*0.044;
Mdl = KDTreeSearcher(fData);
toc
tic
[~,D] = knnsearch(Mdl, fData, 'k',2);
toc
Dist_all2all_left = D(:,2);

clear fData
tic
fData(:,1) = rData.filteredData(logical(rData.maskedIdx),1)*0.0237;
fData(:,2) = rData.filteredData(logical(rData.maskedIdx),2)*0.0237;
fData(:,3) = rData.filteredData(logical(rData.maskedIdx),3)*0.044;
Mdl = KDTreeSearcher(fData);
toc
tic
[~,D] = knnsearch(Mdl, fData, 'k',2);
toc
Dist_all2all_right = D(:,2);

%% distance - DAN 2 DAN
clear fData
tic
fData(:,1) = lData.filteredData(logical(lData.maskedIdx)& lData.DANSynapseIdx,1)*0.0237;
fData(:,2) = lData.filteredData(logical(lData.maskedIdx)& lData.DANSynapseIdx,2)*0.0237;
fData(:,3) = lData.filteredData(logical(lData.maskedIdx)& lData.DANSynapseIdx,3)*0.044;
Mdl = KDTreeSearcher(fData);
toc
tic
[~,D] = knnsearch(Mdl, fData, 'k',2);
toc
Dist_DAN2DAN_left = D(:,2);

clear fData
tic
fData(:,1) = rData.filteredData(logical(rData.maskedIdx)& rData.DANSynapseIdx,1)*0.0237;
fData(:,2) = rData.filteredData(logical(rData.maskedIdx)& rData.DANSynapseIdx,2)*0.0237;
fData(:,3) = rData.filteredData(logical(rData.maskedIdx)& rData.DANSynapseIdx,3)*0.044;
Mdl = KDTreeSearcher(fData);
toc
tic
[~,D] = knnsearch(Mdl, fData, 'k',2);
toc
Dist_DAN2DAN_right = D(:,2);
%% recalculate density

ef = 4.09;
px = .097; %in um
zAniso = 0.18/px;
V = 1;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % radius in expanded samples; in um
bandwidth2 = re/px;


rData.filteredData_zcorr(:,1) = rData.filteredData(:,1);
rData.filteredData_zcorr(:,2) = rData.filteredData(:,2);
rData.filteredData_zcorr(:,3) = rData.filteredData(:,3)*zAniso;

tic
Mdl = KDTreeSearcher(rData.filteredData_zcorr); % all nonDAN synapse positions
toc
tic
IdxKDT = rangesearch(Mdl,rData.filteredData_zcorr,bandwidth2);% Search radius = bandwidth2;
toc
rData.filteredsynapseDensity = cellfun(@numel, IdxKDT)/V;



%%
ha = setupFigure(2,1, 'AxesWidth', 4, 'AxesHeight', 1.5,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


axes(ha(1))
Didx = rData.DANSynapseIdx==1 & logical(rData.maskedIdx);
NDidx = rData.DANSynapseIdx==0 & logical(rData.maskedIdx);
GU_stairsHistPlot({single(rData.filteredsynapseDensity(NDidx)),[],[],single(rData.filteredsynapseDensity(Didx))}, 'ShowECDF', true, 'BinSize',1,'CreateFigure', false,...
    'PercentileEnd', 99.9, 'LineWidth',0.75)

title('Alpha 3 Density of nc82')
legend(['n = ' num2str(numel(single(rData.filteredsynapseDensity(NDidx))))], 'Non-DAN', ['n = ' num2str(numel(single(rData.filteredsynapseDensity(Didx))))], 'DAN')
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 25]);
yyaxis left
ylim([0 0.15]);
xticks(0:10:25)
grid off
ylabel('Relative Frequency');
xlabel('Density (# nc82 punca per µm^3)');

axes(ha(2))
Didx = rData.DANSynapseIdx(logical(rData.maskedIdx))==1;
NDidx = rData.DANSynapseIdx(logical(rData.maskedIdx))==0;

GU_stairsHistPlot({single(Dist_all2all_right(NDidx))*1000,[],[],single(Dist_all2all_right(Didx))*1000,single(Dist_DAN2DAN_right)*1000}, 'ShowECDF', true, 'BinSize',20,'CreateFigure', false,...
    'PercentileEnd', 99.9, 'LineWidth',0.75)
title('Alpha 3 Distance to closest nc82 puncta')
legend(['n = ' num2str(numel(Dist_all2all_right(NDidx)))], 'Non-DAN', ['n = ' num2str(numel(Dist_all2all_right(Didx)))], 'DAN2Any',...
    ['n = ' num2str(numel(Dist_DAN2DAN_right))], 'DAN2DAN')
legend('boxoff')
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000]);
yyaxis left
ylim([0 0.008]);
xticks(0:250:1000)
grid off
ylabel('Relative Frequency');
xlabel('Distance (nm)');

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' alpha3DensityDistancePlots_rightONLY.eps']);

%% write amira spots

tic
GU_amiraWriteSpots('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\alpha3Points',rData.filteredData(logical(rData.maskedIdx),:), ...
    rData.filteredsynapseDensity(logical(rData.maskedIdx)), ...
    rData.DANSynapseIdx(logical(rData.maskedIdx)), [],'scales', [0.0237,0.0237,0.044]);
toc

%% make masked clean data
rightMask = readtiff(['D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex14_alpha3_right_z0p18\mask' filesep 'channel 1 [3956, 17001, 10986]-1-mask.tif']);
imDAN = readtiff('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex14_alpha3_right_z0p18\Analysis\channel 0 [3956, 17001, 10986]-1_clean.tif');
imnc82 = readtiff('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex14_alpha3_right_z0p18\Analysis\channel 1 [3956, 17001, 10986]-1_clean.tif');

rightMask = logical(rightMask);
imDAN_masked = zeros(size(rightMask),'uint16');
imDAN_masked(rightMask) = imDAN(rightMask);
imnc82_masked = zeros(size(rightMask),'uint16');
imnc82_masked(rightMask) = imnc82(rightMask);
 % write masked clean data
 writetiff(imDAN_masked, 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex14_alpha3_right_z0p18\Analysis\channel 0 [3956, 17001, 10986]-1_clean_masked.tif');
 writetiff(imnc82_masked, 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex14_alpha3_right_z0p18\Analysis\channel 1 [3956, 17001, 10986]-1_clean_masked.tif');