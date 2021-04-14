%% check test samples

s(1) = GU_calcGaussianIntegral(1, [2.5,2.5,2.5]);
s(2) = GU_calcGaussianIntegral(1, [2.9,2.9,2.9]);
s(3) = GU_calcGaussianIntegral(1, [3.2 3.2 3.2]);
s(4) = GU_calcGaussianIntegral(1, [3.6 3.6 3.6]);
s(5) = GU_calcGaussianIntegral(1, [4,4,4]);

clc
tic
% parfor_progress(5);
parfor j = 1:4
    for i = 5%1:12
        tic
        GU_ExM_calcPunctaDensityVol(data(i).framePaths{1}{3}, data(i).framePaths{1}{1},...
            'FindLocalMaxima', false, 'MinThreshold', [1450,2200],'Verbose', true,'OTSUGaussKernel', 0.5,...
            'MinVoxelVolume', [round(s(j)), 1008],'ExpansionFactor', 3.62); % sigma 3.9
        toc
        %         parfor_progress;
    end
    %     parfor_progress(0);
end
toc
% T = [1450, 2200]

%%  calculate the number of x,y,z tiles
rt_vol1 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c2';
rt_vol2 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c0';
mx = 10580;
my = 66496;
mz = 4488;
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

%% janelia cluster submit - mouse post-syn detection and yfp cell fill

for i = 5321:10640
    submit = ['bsub -n 4 -R"affinity[core(1)]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mouse_PunctaDensity ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

for i = 1:5320
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mouse_PunctaDensity ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end
%% find failed jobs
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/DataStructures/';
% rt = 'U:\igor\Rui\Wholebrainsynapse\ch1\Analysis\DataStructures\';
fnlist = dir([rt '*.mat']);
fnlist = {fnlist.name};
CompletedJobsIdx = zeros(1, 10640);
for ex = 1:10640
    fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_516_densityData.mat'];
    CompletedJobsIdx(ex) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);

for i = job2resubmit
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mouse_PunctaDensity ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% combine the homer signals

% rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/DataStructures/';
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/DataStructures/';
allYFPsynIDX = [];

allData = [];
allDensity=[];
allDANidx = [];

parfor_progress(10640);
for i = 1:10640
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_densityData.mat'];
    variableInfo = who('-file', [rt filesep fn]);
    PF = ismember('filteredsynapseDensity', variableInfo) & ismember('filteredData', variableInfo) & ismember('Vol2PunctaIdx', variableInfo);
    if PF
        t = load([rt filesep fn], 'filteredsynapseDensity', 'filteredData','Vol2PunctaIdx');
        idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
            t.filteredData(:,1)<=725 &  t.filteredData(:,2)<=725 & t.filteredData(:,3)<=725;
        
        t.filteredData = t.filteredData(idx,:);
        t.filteredsynapseDensity = t.filteredsynapseDensity(idx);
        t.Vol2PunctaIdx = t.Vol2PunctaIdx(idx);
        
        t.filteredData(:,1) = t.filteredData(:,1) + xmin(i)-1;
        t.filteredData(:,2) = t.filteredData(:,2) + ymin(i)-1;
        t.filteredData(:,3) = t.filteredData(:,3) + zmin(i)-1;
        allData = [allData; t.filteredData];
        allDensity = uint8([allDensity; t.filteredsynapseDensity]);
        allDANidx = logical([allDANidx; t.Vol2PunctaIdx]);
    end
    parfor_progress;
end
parfor_progress(0);



%% generate amira point cloud file
load('20171104_CombinedHomorCoordinates_YFPisDAN.mat');
tic
GU_amiraWriteSpots('/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/MouseBrainHomerDensity',allData, allDensity, allDANidx, allDANidx+1,'scales', [0.0237,0.0237,0.044]);
toc

% tic
% GU_amiraWriteSpots('/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/WholeBraineDANSynpaseDensity',allData(allDANidx,:), allDensity(allDANidx), allDANidx(allDANidx),labelidx(allDANidx), 'scales', [0.0237,0.0237,0.044]);
% toc

%% spectral mixing of YFP, YFP-associated Homer, and all Homer

saveFolder = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed';
yfpRT = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008';
yfpSuffix = '_clean';
HyRT = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2Ass';
HySuffix = '_516_clean_V2Ass';
HnyRT = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2nonAss';
HnyRTSuffix = '_516_clean_V2nonAss';

% ex = 5020;
% GU_ExM_CombineMouseYFPHomerCHS(yfpRT, yfpSuffix, HyRT, HySuffix, HnyRT, HnyRTSuffix, saveFolder, ex)
% cluster submit
for i = 1:10640
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "Mix' num2str(i) '" -o ~/logs/Mix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_CombineMouseYFPHomerCHS ' yfpRT ' ' yfpSuffix ' ' HyRT ' ' HySuffix ' ' HnyRT ' ' HnyRTSuffix ' ' saveFolder ' ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end


%% janelia cluster submit - mouse post-syn detection and yfp cell fill
i = 1;
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/Ex7_SingleLargeCell/';
fnv1 = [rt 'ch2.tif'];
fnv2 = [rt 'ch0.tif'];
submit = ['bsub -n 32 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "BigJob' num2str(i) '" -o ~/logs/BigJob' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mosuse_PunctaDensity_LARGEVOL ' fnv1 ' ' fnv2 ''''];
[stat, res] = system(submit,'-echo');

%% dendridic spine analysis
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};
for i = 1:numel(fnlist)
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Spine' num2str(i) '" -o ~/logs/Spine' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_extractSpineParameters ' [rt fnlist{i}] ''''];
    [stat, res] = system(submit,'-echo');
    
end

%%
%% rotated size
% see GU_ExM_MouseDataRotation.m


%% rotate tiles of the spectrally mixed mouse volume
mx = 10580;
my = 66496;
mz = 4488;
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
% rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/';
% for i = 1:nn-1
%     fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
%     submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(i) '" -o ~/logs/Rot' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt fn]  ''''];
%     [stat, res] = system(submit,'-echo');
% end
clear rt
% re-running Nov 28, cubic interpolation
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2Ass/';
a = 0;
for k = 1:2
    for i = 1:nn-1
        a = a +1;
        if k == 1
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_clean_V2Ass.tif'];
        end
        %         if k == 2
        %             fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_clean.tif'];
        %         end
        if k == 2
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean_V2Ass.tif'];
        end
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(a) '" -o ~/logs/Rot' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end
%% find failed jobs
clear rt
% rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/RotatedStacks/';
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/RotatedStacks/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2Ass/RotatedStacks/';
CompletedJobsIdx = cell(1, 3);
for k = 1:2
    
    fnlist = dir([rt{k} '*.tif']);
    fnlist = {fnlist.name};
    
    for i = 1:nn-1
        if k == 1
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_clean_V2Ass.tif'];
        end
        if k==2
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean_V2Ass.tif'];
        end
        CompletedJobsIdx{k}(i) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
    end
    job2resubmit{k} = find(CompletedJobsIdx{k}==0);
end


clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/';
a = 0;
for k = 1:2
    
    for i = job2resubmit{k}
        a = a +1;
        if k == 1
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_clean_V2Ass.tif'];
        end
        if k ==2
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean_V2Ass.tif'];
        end
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(a) '" -o ~/logs/Rot' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end
%%
mx = 10580;
my = 66496;
mz = 4488;
omz = mz;
zxRatio = 0.18/0.097;
theta = 32;
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/RotatedStacks_cubicInterp/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/RotatedStacks_cubicInterp/';
fileSufix = {'_246_clean', '_clean', '_516_clean'};
for k = 1:3
    GU_ExM_WriteRotatedCoordinates(rt{k}, mx, my, mz, zxRatio, theta, fileSufix{k})
end

%% copy only non-empty tiles from the spectrally mixed data
rtmat = 'X:\4Stephan\171104_Mousebrainsynapse\SpectralMixed\RotatedStacks\params\';
rtO = 'X:\4Stephan\171104_Mousebrainsynapse\SpectralMixed\RotatedStacks\';

rt = 'X:\4Stephan\171104_Mousebrainsynapse\ch0\Analysis\1008\RotatedStacks_cubicInterp\allRenamedTiles\';
rt2 = 'X:\4Stephan\171104_Mousebrainsynapse\ch2\Analysis\516\RotatedStacks_cubicInterp\allRenamedTiles\';
rtmv = 'X:\4Stephan\171104_Mousebrainsynapse\ch0\Analysis\1008\RotatedStacks_cubicInterp\nonEmpty\';
rtmv2 = 'X:\4Stephan\171104_Mousebrainsynapse\ch2\Analysis\516\RotatedStacks_cubicInterp\nonEmpty\';
mkdir(rtmv)
mkdir(rtmv2)
parfor_progress(nn-1);
parfor i = 1:nn-1
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
    fnMat = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.mat'];
    if exist([rtO fnO], 'file') && exist([rtmat fnO(1:end-3) 'mat'], 'file')
        Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        if Proceed.nVoxOccupied > 0
            IMinfo = imfinfo([rtO fnO]);
            hrsx = IMinfo(1).Width/2;
            hrsy = IMinfo(1).Height/2;
            hrsz = numel(IMinfo)/2;
            fnN = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_Rotated.tif'];
            copyfile([rt fnN], [rtmv fnN]);
            copyfile([rt2 fnN], [rtmv2 fnN]);
            FE(i) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);

%% copy only non-empty tiles from the spectrally mixed data
rtmat = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/params/';
rtO = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/';

rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/RotatedStacks/';
rt2 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2Ass/RotatedStacks/';
% rt3 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/RotatedStacks_cubicInterp/allRenamedTiles/';
rtmv = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/RotatedStacks/nonEmpty/';
rtmv2 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516_V2Ass/RotatedStacks/nonEmpty/';
% rtmv3 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/RotatedStacks_cubicInterp/nonEmpty/';
mkdir(rtmv)
mkdir(rtmv2)
% mkdir(rtmv3)
FE = zeros(1,nn-1);
parfor_progress(nn-1);
parfor i = 1:nn-1
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
    fnMat = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.mat'];
    fn1 = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_clean_V2Ass.tif'];
    fn2 = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean_V2Ass.tif'];
    if exist([rtO fnO], 'file') && exist([rtmat fnO(1:end-3) 'mat'], 'file')
        Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        if Proceed.nVoxOccupied > 0
            IMinfo = imfinfo([rtO fnO]);
            hrsx = IMinfo(1).Width/2;
            hrsy = IMinfo(1).Height/2;
            hrsz = numel(IMinfo)/2;
            fnN = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_Rotated.tif'];
            copyfile([rt fn1], [rtmv fnN]);
            copyfile([rt2 fn2], [rtmv2 fnN]);
            %             copyfile([rt3 fnN], [rtmv3 fnN]);
            FE(i) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);
%% delete tiles from 246 not in 516 set
rtmat = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/params/';
rtO = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/RotatedStacks_cubicInterp/nonEmpty/';
rtmv2 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/nonEmpty/';
parfor_progress(numel(fnlist2));

fnlist = dir([rtO '*.tif']);
fnlist = {fnlist.name}';

fnlist2 = dir([rtmv2 '*.tif']);
fnlist2 = {fnlist2.name}';
FD = zeros(1, numel(fnlist));

mkdir([rtmv2  'todelete'])

parfor i = 1:numel(fnlist2)
    if exist([rtmv2 fnlist2{i}], 'file') && ~exist([rtO fnlist2{i}], 'file')
        FD(i) = 1;
        movefile([rtmv2  fnlist2{i}], [rtmv2 'todelete' filesep fnlist2{i}])
    end
    parfor_progress;
end
parfor_progress(0);

%% find missing file
rtmat = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/params/';
rtO = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/516/RotatedStacks_cubicInterp/nonEmpty/';
rtmv2 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/nonEmpty/';
parfor_progress(numel(fnlist2));

fnlist = dir([rtO '*.tif']);
fnlist = {fnlist.name}';

fnlist2 = dir([rtmv2 '*.tif']);
fnlist2 = {fnlist2.name}';
FD = zeros(1, numel(fnlist));

mkdir([rtmv2  'todelete'])

parfor i = 1:numel(fnlist)
    if ~exist([rtmv2 fnlist{i}], 'file') && exist([rtO fnlist{i}], 'file')
        FD(i) = 1;
        %         movefile([rtmv2  fnlist2{i}], [rtmv2 'todelete' filesep fnlist2{i}])
    end
    parfor_progress;
end
parfor_progress(0);


%% copy missing file

rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/';
rtmv2 = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/nonEmpty/';
mi = zeros(10640,1);
missingFN = '2191.5368,63001,1636.6116_3355.5368,63750,2621.6116_Rotated.tif';
parfor i = 1:10640
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_246_clean.tif'];
    testfn = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_Rotated.tif'];
    if strcmp(missingFN, testfn)
        mi(i) = 1;
        copyfile([rt fnO], [rtmv2 missingFN]);
    end
    
end

%%
% generate cropped tifs
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/nonEmpty/slice-tiff/ch0/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/RotatedStacks_cubicInterp/nonEmpty/slice-tiff/ch0/';

xrmin = 6000;
xrmax = 8000;
yrmin = 60000;
yrmax = 65000;
z = 2000:3000;
for k = 1:numel(rt)
    for i = z
        fn = [rt{k} num2str(i) '.tif'];
        %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

%% cubic smoothed cleaning of stitched data
minVol = 5000;
for i = 500:3200
    
    %    GU_ExM_JaneliaCluster_MergeMouseSomato_stitchedSlices(ex, minVol)
    submit = ['bsub -n 32 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "mouse' num2str(i) '" -o ~/logs/Rot' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MergeMouseSomato_stitchedSlices ' num2str(i) ' ' num2str(minVol)  ''''];
    [stat, res] = system(submit,'-echo');
end

%%
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed_cubicInterp_Clean/';
minVol = 5000;
for i = 500:3200
    if ~exist([rt num2str(i) '.tif'], 'file')
        submit = ['bsub -n 32 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "mouse' num2str(i) '" -o ~/logs/mouse' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MergeMouseSomato_stitchedSlices ' num2str(i) ' ' num2str(minVol)  ''''];
        [stat, res] = system(submit,'-echo');
        % FE(i) = 1
    end
end

%%
% generate cropped tifs
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed_cubicInterp_Clean/';

xrmin = 8500;
xrmax = 9100;
yrmin = 42200;
yrmax = 43500;
z = 1750:2030;
for k = 1:numel(rt)
    for i = z
        fn = [rt{k} num2str(i) '.tif'];
        %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

%%

for i = 5321:10640
    submit = ['bsub -n 4 -R"affinity[core(1)]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mouse_ReCalcYFPChannel ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

for i = 1:5320
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_Mouse_ReCalcYFPChannel ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% rotate tiles of the spectrally mixed mouse volume
mx = 10580;
my = 66496;
mz = 4488;
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

clear rt
% re-running Dec 6, cubic interpolation
% rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/';
% rt{3} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/5000/';
% rt{4} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/7500/';
% a = 0;
% rtsuffix = {'_1008', '_2500', '_5000', '_7500'};

rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/';
a = 0;
rtsuffix = {'_1008'};

for k = 1:numel(rt)
    for i = 1:nn-1
        a = a +1;
        
        fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) rtsuffix{k} '_clean.tif'];
        
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(a) '" -o ~/logs/Rot' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

%% copy only non-empty tiles from the spectrally mixed data
rtmat = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/params/';
rtO = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/';
clear rt
clear rtmv
% rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/RotatedStacks/';
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/7500/RotatedStacks/';

% rtmv{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/nonEmpty/';
% rtmv{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/RotatedStacks/nonEmpty/';
rtmv{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/nonEmpty/';
% rtmv{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/7500/RotatedStacks/nonEmpty/';

for k = 1:numel(rt)
    mkdir(rtmv{k})
end

rtsuffix = {'_1008', '_2500', '_5000', '_7500'};
% rtsuffix = {'_5000', '_7500'};
FE = zeros(1,nn-1);

parfor_progress(nn-1);
parfor i = 1:nn-1
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
    fnMat = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.mat'];
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i))];
    if exist([rtO fnO], 'file') && exist([rtmat fnO(1:end-3) 'mat'], 'file')
        Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        if Proceed.nVoxOccupied > 0
            IMinfo = imfinfo([rtO fnO]);
            hrsx = IMinfo(1).Width/2;
            hrsy = IMinfo(1).Height/2;
            hrsz = numel(IMinfo)/2;
            fnN = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_Rotated.tif'];
            for k = 1:numel(rt)
                copyfile([rt{k} fn rtsuffix{k} '_clean.tif'], [rtmv{k} fnN]);
            end
            FE(i) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);
%%
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/RotatedStacks/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/5000/RotatedStacks/';
rt{4} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/7500/RotatedStacks/';
CompletedJobsIdx = cell(1, numel(rt));
rtsuffix = {'_1008', '_2500', '_5000', '_7500'};
for k = 1:numel(rt)
    
    fnlist = dir([rt{k} '*.tif']);
    fnlist = {fnlist.name};
    
    for i = 1:nn-1
        fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) rtsuffix{k} '_clean.tif'];
        
        CompletedJobsIdx{k}(i) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
    end
    job2resubmit{k} = find(CompletedJobsIdx{k}==0);
end



clear rt
% re-running Nov 28, cubic interpolation
rt{1} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/2500/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/5000/';
rt{4} = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/7500/';
a = 0;


for k = 1:numel(rt)
    for i = job2resubmit{k}
        a = a +1;
        
        fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) rtsuffix{k} '_clean.tif'];
        
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(a) '" -o ~/logs/Rot' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end
%%
rt = '/groups/betzig/betziglab/4Stephan/171205_SomatoSensoryCortextestconditions/0-Raw/';
fn = 'Raw_channel 0 [3196, 34491, 6488].tif';
T = [1800, 1800, 1800 1800];
MV= [246, 516, 1008 2000];
parfor i = 1:numel(T)
    GU_ExM_calcPunctaDensityVol_1ch([rt fn],...
        'FindLocalMaxima', false, 'MinThreshold', [1850],'Threshold', T(i), 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', MV(i),'ExpansionFactor', 3.62);
end

%% rename files
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed_cubicInterp_Clean_500_3200/';

for i = 1:2701
    movefile([rt 'z' num2str(i-1,'%04d') '_c0_t0.tif'], [rt num2str(i+499) '.tif'])
end
%% spectral mixing
for i = 500:3200
    
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "mouseSM' num2str(i) '" -o ~/logs/mouseSM' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_SpectralMixMouseSomato ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
    
end
%%
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/LinearYFP_cubicHomer_SM/';
tic
for i = 1:2701
    movefile( [rt num2str(i+499) '.tif'],[rt 'z' num2str(i-1,'%04d') '_c0_t0.tif'])
end
toc

%% anaylsis of good homer basson datasets
%% anaylyze v5 and nc82
rt = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/stitchedata/';
fnv1 = [rt 'ch0' filesep 'ch0.tif'];
fnv2 = [rt 'ch1' filesep 'ch1.tif'];
fnv3 = [rt 'ch2' filesep 'ch2.tif'];
MV = [1008,1008,1008];
fn = {fnv1, fnv2, fnv3};
T = [1000,500,4453];
GU_ExM_calcPunctaDensityVol(fnv2, fnv1,...
    'FindLocalMaxima', false, 'Threshold', [T(2), T(1)],'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [1008, 1008],'ExpansionFactor', 4,'minAZdia', 0.1);
% Syd1 + PN
GU_ExM_calcPunctaDensityVol(fnv3, fnv1,...
    'FindLocalMaxima', false, 'Threshold', [T(3), T(1)],'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [1008, 1008],'ExpansionFactor', 4,'minAZdia', 0.1);

%% redo
rt2 = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\';
fnv2 = [rt2 'ch1_v2' filesep 'ch1_1008_clean_v2.tif'];


GU_ExM_calcPunctaDensityVol_1ch(fnv2,...
    'ExpansionFactor', 4.04,'FindLocalMaxima', false, 'Threshold', 0,'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', 246,'MinVolumeOccupancy',0);


%% get thresholds from OTSU - crop
fn = {'ch0_crop.tif', 'ch1_crop.tif', 'ch2_crop.tif'};
for i = 1:3
    tic
    im = readtiff(fn{i});
    toc
    
    tic
    imG = filterGauss3D(im, 0.5);toc
    
    tic
    T(i) = thresholdOtsu(imG(im>0 & im<prctile(im(:),99)));toc
    
end

%%
s32rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\segment32\tiff\';
% im32 = zeros(6762, 6947, 1820, 'uint16');
fn32 = dir([s32rt '*.tif']);
% fn32 = {fn32.name}';
[fn32,~,~] = natsortfiles({fn32.name}');
% parfor i = 1:numel(fn32)
%     im32(:,:,i) = readtiff([s32rt fn32{i}]);
% end


HomerAss32 = zeros(6762, 6947, 1820, 'logical');
HomerAll32 = zeros(6762, 6947, 1820, 'logical');
YFP32 = zeros(6762, 6947, 1820, 'logical');

% yfp
tic
parfor_progress(numel(fn32));
parfor i = 1:numel(fn32)
    tmp = readtiff([s32rt fn32{i}]);
    yfp = zeros(6762, 6947, 1, 'logical');
    yfp(tmp<32768) = logical(tmp(tmp<32768));
    YFP32(:,:,i) = yfp;
    parfor_progress;
end
parfor_progress(0);
toc

writetiff(uint8(YFP32), [s32rt 'YFP32.tif']);

% yfp ass homer
tic
parfor_progress(numel(fn32));
parfor i = 1:numel(fn32)
    tmp = readtiff([s32rt fn32{i}]);
    yasshomer =  zeros(6762, 6947, 1, 'logical');
    yasshomer(tmp>=32768 & tmp<49152) = logical(tmp(tmp>=32768 & tmp<49152));
    HomerAss32(:,:,i) = yasshomer;
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(uint8(HomerAss32), ['HomerAss32.tif']);
clear HomerAss32

% all homer
tic
parfor_progress(numel(fn32));
parfor i = 1:numel(fn32)
    tmp = readtiff([s32rt fn32{i}]);
    allhomer  = zeros(6762, 6947, 1, 'logical');
    allhomer(tmp>49152) = logical(tmp(tmp>49152));
    HomerAll32(:,:,i) = allhomer;
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(uint8(HomerAll32), ['HomerAll32.tif']);

% segment 32 coordiantes
s32_x = [2816:9763];
s32_y = [29352:36114];
s32_z = [1:1820];

HomerAll32 = logical(readtiff(['HomerAll32.tif']));
HomerAll32_rp3 = regionprops3(HomerAll32);
HomerAss32 = logical(readtiff(['HomerAss32.tif']));
HomerAss32_rp3 = regionprops3(HomerAss32);
clear HomerAss32 HomerAll32
% MinVoxelVolume = 1008;
% YFP32 = logical(readtiff(['YFP32.tif']));
%  CC = bwconncomp(YFP32, 26);
%         csize = cellfun(@numel, CC.PixelIdxList);
%         idx = csize>=MinVoxelVolume;
%         CC.NumObjects = sum(idx);
%         CC.PixelIdxList = CC.PixelIdxList(idx);
%         YFP32_clean = labelmatrix(CC)~=0;
% [y,x,z] = find(YFP32_clean);
% writetiff(uint8(YFP32_clean), [ 'YFP32_clean.tif']);
% YFP32 = logical(readtiff(['YFP32_clean.tif']));
%
% yfpPerim = bwperim(YFP32_clean);
% writetiff(uint8(yfpPerim), [ 'yfpPerim.tif']);
% [py,px,pz] = find(yfpPerim);
%
% YFP32_Perim = logical(readtiff(['yfpPerim.tif']));

% yfp - read to clean (info gathered from the mip)
YFP32 = zeros(6762, 6947, 1820, 'uint16');
tic
parfor_progress(numel(fn32));
parfor i = 1:numel(fn32)
    tmp = readtiff([s32rt fn32{i}]);
    yfp = zeros(6762, 6947,1, 'uint16');
    yfp(tmp<32768) = tmp(tmp<32768);
    YFP32(:,:,i) = yfp;
    parfor_progress;
end
parfor_progress(0);
% P = prctile(YFP32(YFP32>0),99.5);
% T = thresholdOtsu(YFP48(YFP32>0 & YFP32<P)); %see T for 48, using the
T = 2500;
YFP32 = YFP32-T;
writetiff(YFP32, ['YFP32_OTSU.tif']);

YFP32 = logical(YFP32);
MinVoxelVolume = 1008;
CC = bwconncomp(YFP32, 26);
csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinVoxelVolume;
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
YFP32 = labelmatrix(CC)~=0;

% tic ; YFP32 = logical(readtiff([ 'YFP32_OTSU.tif'])); toc
nVoxSeg32 = sum(YFP32(:));
yfpPerim = bwperim(YFP32);
nVoxSeg32_Perim = sum(yfpPerim(:));

tic
[y,x,z] = ind2sub(size(YFP32), find(YFP32==1));
YFP32_VolcoordXYZ(:,1) = x';
YFP32_VolcoordXYZ(:,2) = y';
YFP32_VolcoordXYZ(:,3) = z';
toc

tic
[y,x,z] = ind2sub(size(yfpPerim), find(yfpPerim==1));
YFP32_PerimcoordXYZ(:,1) = x';
YFP32_PerimcoordXYZ(:,2) = y';
YFP32_PerimcoordXYZ(:,3) = z';
toc

writetiff(uint8(yfpPerim), [ 'yfp32Perim_OTSU.tif']);
writetiff(uint8(YFP32), [ 'YFP32_clean_OTSU_logical.tif']);

clear YFP32_clean yfpPerim YFP32
%%
clear all
s48rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\segment48\tiff\';
cd(s48rt);
% im32 = zeros(6762, 6947, 1820, 'uint16');
fn48 = dir([s48rt '*.tif']);
% fn48 = {fn48.name}';
[fn48,~,~] = natsortfiles({fn48.name}');
% parfor i = 1:numel(fn48)
%     im32(:,:,i) = readtiff([s48rt fn48{i}]);
% end




YFP48 = zeros(9006, 5260, 2075, 'logical');

% yfp
tic
parfor_progress(numel(fn48));
parfor i = 1:numel(fn48)
    tmp = readtiff([s48rt fn48{i}]);
    yfp = zeros(9006, 5260, 1, 'logical');
    yfp(tmp<32768) = logical(tmp(tmp<32768));
    YFP48(:,:,i) = yfp;
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(uint8(YFP48), [s48rt 'YFP48.tif']);


% yfp - read to clean (info gathered from the mip)
YFP48 = zeros(9006, 5260, 2075, 'uint16');
tic
parfor_progress(numel(fn48));
parfor i = 1:numel(fn48)
    tmp = readtiff([s48rt fn48{i}]);
    yfp = zeros(9006, 5260, 1, 'uint16');
    yfp(tmp<32768) = tmp(tmp<32768);
    YFP48(:,:,i) = yfp;
    parfor_progress;
end
parfor_progress(0);
% P = prctile(YFP48(YFP48>0),99.5);
% T = thresholdOtsu(YFP48(YFP48>0 & YFP48<P)); %2536
T = 2500;
YFP48 = YFP48-T;
writetiff(YFP48, [ 'YFP48_OTSU.tif']);

YFP48 = logical(YFP48);
MinVoxelVolume = 1008;
CC = bwconncomp(YFP48, 26);
csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinVoxelVolume;
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
YFP48 = labelmatrix(CC)~=0;


nVoxSeg48 = sum(YFP48(:));
yfpPerim = bwperim(YFP48);
nVoxSeg48_Perim = sum(yfpPerim(:));

tic
[y,x,z] = ind2sub(size(YFP48), find(YFP48));
YFP48_VolcoordXYZ(:,1) = x';
YFP48_VolcoordXYZ(:,2) = y';
YFP48_VolcoordXYZ(:,3) = z';
toc

tic
[y,x,z] = ind2sub(size(yfpPerim), find(yfpPerim));
YFP48_PerimcoordXYZ(:,1) = x';
YFP48_PerimcoordXYZ(:,2) = y';
YFP48_PerimcoordXYZ(:,3) = z';
toc

writetiff(uint8(yfpPerim), [ 'yfp48Perim_OTSU.tif']);
writetiff(uint8(YFP48), [ 'YFP48_clean_OTSU_logical.tif']);

clear YFP32_clean yfpPerim YFP48

tic
yfpPerim =readtiff('yfp48Perim.tif'); toc
% yfp ass homer
tic
HomerAss48 = zeros(9006, 5260, 2075, 'logical');
parfor_progress(numel(fn48));
parfor i = 1:numel(fn48)
    tmp = readtiff([s48rt fn48{i}]);
    yasshomer =  zeros(9006, 5260, 1, 'logical');
    yasshomer(tmp>=32768 & tmp<49152) = logical(tmp(tmp>=32768 & tmp<49152));
    HomerAss48(:,:,i) = yasshomer;
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(uint8(HomerAss48), ['HomerAss48.tif']);
HomerAss48 = readtiff(['HomerAss48.tif']);
HomerAss48_rp3 = regionprops3(HomerAss48);
clear HomerAss48

% all homer
tic
HomerAll48 = zeros(9006, 5260, 2075, 'logical');
parfor_progress(numel(fn48));
parfor i = 1:numel(fn48)
    tmp = readtiff([s48rt fn48{i}]);
    allhomer  = zeros(9006, 5260, 1, 'logical');
    allhomer(tmp>49152) = logical(tmp(tmp>49152));
    HomerAll48(:,:,i) = allhomer;
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(uint8(HomerAll48), ['HomerAll48.tif']);
HomerAll48_rp3 = regionprops3(HomerAll48);

%% vol calc
tic
YFP32 = readtiff('D:\GoogleDrive_DataTransfer\exllsmdata_sharedwithgokulrui\SomatoSensoryCortex_data_homer\segment32\YFP32_clean.tif'); toc

tic
YFP48 = readtiff('D:\GoogleDrive_DataTransfer\exllsmdata_sharedwithgokulrui\SomatoSensoryCortex_data_homer\segment48\YFP48_clean.tif'); toc

tic
[y,x,z] = ind2sub(size(YFP32), find(YFP32));
YFP32_VolcoordXYZ(:,1) = x';
YFP32_VolcoordXYZ(:,2) = y';
YFP32_VolcoordXYZ(:,3) = z';
toc


%% rotated voxel dimensions
rx = 114.38;
ry = 97;
rz = 155.23;

[y,x,z] = ind2sub(size(YFP32), find(YFP32_Perim));
YFP32_perimcoordXYZ(:,1) = x';
YFP32_perimcoordXYZ(:,2) = y';
YFP32_perimcoordXYZ(:,3) = z';

% read seg48 perim volume before running snipet below
[y,x,z] = ind2sub(size(yfpPerim), find(yfpPerim));
YFP48_perimcoordXYZ(:,1) = x';
YFP48_perimcoordXYZ(:,2) = y';
YFP48_perimcoordXYZ(:,3) = z';

%% crop cleaned data
% tkhpc36c
% segment 32 coordiantes
s32_x = [2816:9763];
s32_y = [29352:36114];
s32_z = [746:2565];
% s32_z = [1:1820];

% segment 48 coordiantes
s48_x = [898:6158];
s48_y = [28563:37569];
s48_z = [403:2477];

cx = [898:9762];
cy = [28563:37568];
cz = [403:2565]+500;

allHomerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\allHomer_cubic\ch0\';
YFP_Homerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\yfpAssHomer\ch0\';
YFPrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\YFP\ch0\';
srt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\';
rt = {allHomerrt, YFP_Homerrt, YFPrt};
srtsd = {'CropAllHomer', 'CropYFPassHomer', 'CropYFP'};
for j = 1:numel(rt)
    mkdir([srt srtsd{j}]);
end

parfor_progress(numel(cz));
parfor i = cz
    for j = 1:numel(rt)
        im = readtiff([rt{j} num2str(i) '.tif']);
        im = im(cy, cx);
        writetiff(im, [srt srtsd{j} filesep num2str(i) '.tif'])
    end
    parfor_progress;
end
parfor_progress(0);

%% bin yfp data
YFPrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\YFP\ch0\';
fn = dir([YFPrt '*.tif']);
fn = {fn.name}';
srt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\YFPbin_10percent\';
mkdir(srt)
parfor_progress(numel(fn));
parfor i = 1:numel(fn)
    
    im = readtiff([YFPrt fn{i}]);
    im = imresize(im, 0.1);
    im = scaleContrast(im, [0 35000], [0 255]);
    writetiff(uint8(im), [srt fn{i}]);
    
    parfor_progress;
end
parfor_progress(0);

%% mask cropped clean data
% tkhpc48
% %
% % Figure list
% % 1. raw data --> already done
% % 2. two segmented neurons
% % --> generate masks
% % --> generate masked neurons
% % --> generate unmasked data
% % 3. generate PM vs internal Homer
% % --> PM: errode 21px (500nm) on the masked regions
% % --> Internal: anything remaining after the errosion --> need to map based on centroids

% Crop_allHomerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropAllHomer\';
% Crop_YFP_Homerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropYFPassHomer\';
% Crop_YFPrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropYFP\';
% Crop_maskedNeurons = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\segmets32_48_compressed\';
%
% Crop_mallHomerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropAllHomer\';
% Crop_mYFP_Homerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropYFPassHomer\';
% Crop_mYFPrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropYFP\';
% Crop_mrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\logicalNeuronMasks\';
% Crop_imallHomerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropAllHomer\';
% Crop_imYFP_Homerrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropYFPassHomer\';
% Crop_imYFPrt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropYFP\';

Crop_allHomerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropAllHomer\';
Crop_YFP_Homerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropYFPassHomer\';
Crop_YFPrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\CropYFP\';
Crop_maskedNeurons = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\segmets32_48_compressed\';

Crop_mallHomerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropAllHomer\';
Crop_mYFP_Homerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropYFPassHomer\';
Crop_mYFPrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\maskedCropYFP\';
Crop_mrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\logicalNeuronMasks\';
Crop_imallHomerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropAllHomer\';
Crop_imYFP_Homerrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropYFPassHomer\';
Crop_imYFPrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\InvMaskedCropYFP\';

rt = {Crop_mallHomerrt, Crop_mYFP_Homerrt, Crop_mYFPrt, Crop_mrt, Crop_imallHomerrt, Crop_imYFP_Homerrt, Crop_imYFPrt};
for j = 1:numel(rt)
    mkdir(rt{j})
end

cfn_allHomer = dir([Crop_allHomerrt '*.tif']);
[cfn_allHomer,~,~] = natsortfiles({cfn_allHomer.name}');
cfn_YFP_Homer = dir([Crop_YFP_Homerrt '*.tif']);
[cfn_YFP_Homer,~,~] = natsortfiles({cfn_YFP_Homer.name}');
cfn_YFP = dir([Crop_YFPrt '*.tif']);
[cfn_YFP,~,~] = natsortfiles({cfn_YFP.name}');
cfn_makedNeurons = dir([Crop_maskedNeurons '*.tif']);
[cfn_makedNeurons,~,~] = natsortfiles({cfn_makedNeurons.name}');

parfor_progress(numel(cfn_allHomer));
parfor i = 1:numel(cfn_allHomer)
    im_m = readtiff([Crop_maskedNeurons cfn_makedNeurons{i}]);
    im_aH = readtiff([Crop_allHomerrt cfn_allHomer{i}]);
    im_YH = readtiff([Crop_YFP_Homerrt cfn_YFP_Homer{i}]);
    im_Y = readtiff([Crop_YFPrt cfn_YFP{i}]);
    
    m = logical(im_m); %mask
    
    aim_aH = im_aH;
    aim_aH(m(:)) = 0;
    im_aH(~m(:)) = 0;
    
    aim_YH = im_YH;
    aim_YH(m(:)) = 0;
    im_YH(~m(:)) = 0;
    
    aim_Y = im_Y;
    aim_Y(m(:)) = 0;
    im_Y(~m(:)) = 0;
    
    writetiff(uint8(m), [Crop_mrt num2str(i) '.tif']);
    writetiff(aim_aH, [Crop_imallHomerrt num2str(i) '.tif']);
    writetiff(im_aH, [Crop_mallHomerrt num2str(i) '.tif']);
    writetiff(aim_YH, [Crop_imYFP_Homerrt num2str(i) '.tif']);
    writetiff(im_YH, [Crop_mYFP_Homerrt num2str(i) '.tif']);
    writetiff(aim_Y, [Crop_imYFPrt num2str(i) '.tif']);
    writetiff(im_Y, [Crop_mYFPrt num2str(i) '.tif']);
    parfor_progress;
end
parfor_progress(0);

%% write files names for Imaris File conversion

srt_merge =  'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\ImarisFiles\';
mkdir(srt_merge);

parfor_progress(numel(cfn_allHomer));
parfor i = 1:numel(cfn_allHomer)
    
    
    aim_aH = readtiff([Crop_imallHomerrt num2str(i) '.tif']);
    im_aH  = readtiff([Crop_mallHomerrt num2str(i) '.tif']);
    aim_YH = readtiff([Crop_imYFP_Homerrt num2str(i) '.tif']);
    im_YH = readtiff([Crop_mYFP_Homerrt num2str(i) '.tif']);
    aim_Y = readtiff([Crop_imYFPrt num2str(i) '.tif']);
    im_Y = readtiff([Crop_mYFPrt num2str(i) '.tif']);
    
    writetiff(aim_aH, [srt_merge 'c0_t0_z' num2str(i,'%04d') '.tif']);
    writetiff(im_aH, [srt_merge 'c1_t0_z' num2str(i,'%04d') '.tif']);
    writetiff(aim_YH, [srt_merge 'c2_t0_z' num2str(i,'%04d') '.tif']);
    writetiff(im_YH, [srt_merge 'c3_t0_z' num2str(i,'%04d') '.tif']);
    writetiff(aim_Y, [srt_merge 'c4_t0_z' num2str(i,'%04d') '.tif']);
    writetiff(im_Y, [srt_merge 'c5_t0_z' num2str(i,'%04d') '.tif']);
    parfor_progress;
end
parfor_progress(0);

%% generate max projections

% srt_merge =  'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\ImarisFiles\';
% mkdir(srt_merge);
% aim_aH = zeros([9006, 8865, numel(cfn_allHomer)], 'uint16');
parfor_progress(numel(cfn_allHomer));
for i = 1:numel(cfn_allHomer)
    
    if i == 1
        aim_aH = readtiff([Crop_imallHomerrt num2str(i) '.tif']);
        im_aH  = readtiff([Crop_mallHomerrt num2str(i) '.tif']);
        aim_YH = readtiff([Crop_imYFP_Homerrt num2str(i) '.tif']);
        im_YH = readtiff([Crop_mYFP_Homerrt num2str(i) '.tif']);
        aim_Y = readtiff([Crop_imYFPrt num2str(i) '.tif']);
        im_Y = readtiff([Crop_mYFPrt num2str(i) '.tif']);
    else
        aim_aH = max(aim_aH, readtiff([Crop_imallHomerrt num2str(i) '.tif']));
        im_aH  = max(im_aH, readtiff([Crop_mallHomerrt num2str(i) '.tif']));
        aim_YH = max(aim_YH, readtiff([Crop_imYFP_Homerrt num2str(i) '.tif']));
        im_YH = max(im_YH, readtiff([Crop_mYFP_Homerrt num2str(i) '.tif']));
        aim_Y = max(aim_Y, readtiff([Crop_imYFPrt num2str(i) '.tif']));
        im_Y = max(im_Y, readtiff([Crop_mYFPrt num2str(i) '.tif']));
    end
    % aim_aH(:,:,i) = readtiff([Crop_imallHomerrt num2str(i) '.tif']);
    parfor_progress;
end
parfor_progress(0);

srt_merge =  'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\SomatoSensoryCortex_data_homer\MIPs\';
mkdir(srt_merge);
writetiff(aim_aH, [srt_merge 'inv_allHomer.tif']);
writetiff(im_aH, [srt_merge 'mask_allHomer.tif']);
writetiff(aim_YH, [srt_merge 'inv_YFPassHomer.tif']);
writetiff(im_YH, [srt_merge 'mask_YFPassHomer.tif']);
writetiff(aim_Y, [srt_merge 'inv_YFP.tif']);
writetiff(im_Y, [srt_merge 'mask_YFP.tif']);

%% supplementary figure to calculate the otsu values
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerCleanupSupp/rawstacks/';
cd(rt);
fn = dir([rt '*.tif']);
fn = {fn.name};
T = zeros(numel(fn),1);
ha = setupFigure(3,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.5 1 1.5], 'YSpace', [1.5 1 1.5]);
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
t = {'Ex1', 'Ex4', 'Ex5'};
r = {451:500, 340:389, 51:100};
for i = 1:numel(fn)
    im = readtiff([rt fn{i}]);
    im = im(:,:,r{i});
    ahs = im(:);
    T(i) = thresholdOtsu(ahs(ahs<prctile(ahs,99.9)));
    % figure
    
    axes(ha(i))
    
    histogram(ahs(ahs<T(i)), 0:10:(T(i)-1), 'FaceColor', cfO, 'EdgeColor',cfO);
    histogram(ahs(ahs>=T(i)), T(i):10:10000, 'FaceColor', cfP, 'EdgeColor',cfP);
    % line([515,517], [0,10^5])
    set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
    ylabel('# events');
    xlabel('pixel intensity');
    xlim([0 10000])
    xticks([0:2000:10000]);
        ylim([10 10^8])
%     legend(['n= <' num2str(T(i)) '-' num2str(numel(ahs(ahs<T(i)))) '--' 'n= >=' num2str(T(i)) '-' num2str(numel(ahs(ahs>=T(i))))])
    title(t{i});
end


f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'OtsuCutoff.eps']);

 
