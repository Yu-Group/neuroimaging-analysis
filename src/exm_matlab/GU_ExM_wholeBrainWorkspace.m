clc
tic
% parfor_progress(12);
parfor i = 1%1:12
    GU_ExM_calcSynapseDensityVol(data(i).framePaths{1}{2}, data(i).framePaths{1}{1},'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);
    %     parfor_progress;
end
% parfor_progress(0);
toc

%%
for i = 1:10
    [~, fn, ~] = fileparts(data(i).framePaths{1}{2});
    load([data(i).source 'Analysis' filesep fn '_densityData.mat'], 'T');
    Thresh(i, 1:2) = T(1:2);
end

%%  calculate the number of x,y,z tiles

mx = 15055;
my = 27964;
mz = 6647;
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
            %             if i <nx
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            %             else
            %                 xmax(nn) = mx;
            %             end
            
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            %             if j <ny
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            %             else
            %                 ymax(nn) = my;
            %             end
            
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            %             if k < nz
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            %             else
            %                 zmax(nn) = mz;
            %             end
            
            nn = nn + 1;
        end
    end
end

%%
%
for ex = 1:nn-1
    % generate the ch0 ch1 sub-volumes
    cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py' ...
        ' --input "/nrs/saalfeld/igor/illumination-correction/Sample1_C1/stitching/decon-warped-export/export.n5" '...
        ' --output "/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch1" '...
        ' --channel 1' ...
        ' --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ...
        ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
    
    cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py' ...
        ' --input "/nrs/saalfeld/igor/illumination-correction/Sample1_C1/stitching/decon-warped-export/export.n5" '...
        ' --output "/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch0" '...
        ' --channel 0' ...
        ' --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ...
        ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
    
    tic
    [~, commandOut1] = system(cmd1);
    toc
    
    tic
    [~, commandOut2] = system(cmd0);
    toc
    
    % run script
    ch0rt = '/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch0/';
    fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
    ch1rt = '/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch1/';
    fn1 = [ch0rt fn];
    fn2 = [ch1rt fn];
    
    GU_ExM_calcSynapseDensityVol(fn2, fn1,'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);
    
    % delete extracted tiff
    delete(char(fn1))
    delete (char(fn2))
    
    fprintf('Completed %d of %d \n', ex,numel(1:nn-1));
end

%% janelia cluster submit

for i = 1:8800
    submit = ['bsub -n 3 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SynapseCalc' num2str(i) '" -o ~/logs/SynCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_DensityCalc ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% check the status of jobs

rt = '/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch1/Analysis/DataStructures/';
% rt = 'U:\igor\Rui\Wholebrainsynapse\ch1\Analysis\DataStructures\';
fnlist = dir([rt '*.mat']);
fnlist = {fnlist.name};
for ex = 1:8800
    fn = [filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_densityData.mat'];
    CompletedJobsIdx(ex) = ismember(fn(2:end),fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);

for i = job2resubmit
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SynapseCalc' num2str(i) '" -o ~/logs/SynCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_DensityCalc ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%%

%% janelia cluster submit

for i = 1:8800
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SynRECalc' num2str(i) '" -o ~/logs/SynCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_CorrectDANDensityCalc ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%%
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DataStructures_Recalc/';
fnlist = dir([rt '*.mat']);
fnlist = {fnlist.name};
for ex = 1:8800
    fn = [filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_clean_densityDataRC.mat'];
    CompletedJobsIdx(ex) = ismember(fn(2:end),fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);

for i = job2resubmit
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SynRECalc' num2str(i) '" -o ~/logs/SynCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_CorrectDANDensityCalc ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% crop stitched z-slices
xmax = 15055;
ymax = 27964;
clear ppp
% rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch0\Analysis\clean\slice-tiff\ch0';
ppp{4} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0';
ppp{2} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/Density/slice-tiff/ch0';
ppp{3} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0';
ppp{1} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DANDensityRecalc/slice-tiff/ch0';
for k = 1:numel(ppp)
    rt = ppp{k};
    fnlist = dir([rt filesep '*.tif']);
    fnlist = {fnlist.name};
    
    for i = 1:numel(fnlist)
        fp = [rt filesep fnlist{i}];
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "Crop' num2str(i) '" -o ~/logs/Crop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_CropWholeFlyBrainSlice ' fp ' ' num2str(xmax) ' ' num2str(ymax) ''''];
        [~, ~] = system(submit,'-echo');
        %     GU_ExM_CropWholeFlyBrainSlice(fp, xmax, ymax)
    end
end

%% get min and max x y coordinates
ppp{4} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0/SliceProperties';
ppp{2} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/Density/slice-tiff/ch0/SliceProperties';
ppp{3} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0/SliceProperties';
ppp{1} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DANDensityRecalc/slice-tiff/ch0/SliceProperties';

k = 4;
rt = ppp{k};
fnlist = dir([rt filesep '*.mat']);
fnlist = {fnlist.name};
parfor_progress(numel(fnlist))
for i = 1:numel(fnlist)
    sliceProperties_nc82(i) = load([ppp{4} filesep fnlist{i}]);
    sliceProperties_DAN(i) = load([ppp{3} filesep fnlist{i}]);
    parfor_progress;
end
parfor_progress(0);
%% generate ldm list

rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch1\Analysis\clean\slice-tiff\ch0\';
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};
frameFilename=['ch0.lst'];
fid = fopen(frameFilename, 'w');

fid = fopen(frameFilename, 'a');
fprintf(fid,['Parameters {' char(13) '\n' ...
    'Raw 0' char(13) '\n' ...
    'Dims 15055 27964 6647' char(13) '\n' ...
    'Size 0 356.8035 0 662.7468 0 292.468' char(13) '\n'...
    'Channel 1' char(13) '\n'...
    '}' char(13) '\n ' char(13) '\n']);
fclose(fid);
for i = 1:numel(fnlist)
    dlmwrite(frameFilename, [num2str(i-1) '.tif' char(13)], '-append', 'delimiter','', 'precision', 16);
end

%% combine coordinates and density
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DataStructures_Recalc';
fnlist = dir([rt filesep '*.mat']);
fnlist = {fnlist.name};

nT = zeros(numel(fnlist),1, 'logical');
parfor_progress(numel(fnlist));
for i = 1:numel(fnlist)
    t = load([rt filesep fnlist{i}],'ProceedForward');
    if t.ProceedForward
        nT(i) = 1;
    end
    parfor_progress;
end
parfor_progress(0);



%% combine coordinates and density
% get file names

mx = 15055;
my = 27964;
mz = 6647;
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
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            nn = nn + 1;
        end
    end
end

%% combine all the data - from the recalculated set
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DataStructures_Recalc';
rt2 = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/DataStructures';
allData = [];
allDensity=[];
allDANDensity = [];
allDANidx = [];
allDANSynapseDensity = [];
allDis_DANsyn2DANsyn = [];
allDis_DANsyn2NONDANsyn = [];
allDis_NONDANsyn2NONDANsyn = [];
allDis = [];
allDANsynIDX_includePeriph = [];

parfor_progress(8800);
for i = 1:8800
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_clean_densityDataRC.mat'];
    fn2 = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_densityData.mat'];
    
    variableInfo = who('-file', [rt filesep fn]);
    PF = ismember('filteredDistances', variableInfo) & ismember('DANsynapseDensity', variableInfo) & ismember('Distance_DANSyn2DANSyn', variableInfo) ...
        & ismember('Distance_DANSyn2NONDANSyn', variableInfo) & ismember('Distance_NONDANSyn2NONDANSyn', variableInfo) & ...
        ismember('filteredsynapseDensity', variableInfo) & ismember('filteredData', variableInfo) & ismember('DANSynapseIdx', variableInfo);
    if PF
        t = load([rt filesep fn], 'filteredsynapseDensity', 'filteredData','DANSynapseIdx', 'filteredDistances', 'DANsynapseDensity','Distance_DANSyn2DANSyn', 'Distance_DANSyn2NONDANSyn', 'Distance_NONDANSyn2NONDANSyn');
        idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
            t.filteredData(:,1)<=725 &  t.filteredData(:,2)<=725 & t.filteredData(:,3)<=725;
        
        t.filteredData = t.filteredData(idx,:);
        t.filteredsynapseDensity = t.filteredsynapseDensity(idx);
        t.DANSynapseIdx = t.DANSynapseIdx(idx);
        t.filteredDistances = t.filteredDistances(idx);
        t.DANsynapseDensity = t.DANsynapseDensity;
        t.Distance_DANSyn2DANSyn = t.Distance_DANSyn2DANSyn;
        t.Distance_DANSyn2NONDANSyn = t.Distance_DANSyn2NONDANSyn;
        t.Distance_NONDANSyn2NONDANSyn = t.Distance_NONDANSyn2NONDANSyn;
        
        
        t.filteredData(:,1) = t.filteredData(:,1) + xmin(i)-1;
        t.filteredData(:,2) = t.filteredData(:,2) + ymin(i)-1;
        t.filteredData(:,3) = t.filteredData(:,3) + zmin(i)-1;
        allData = [allData; t.filteredData];
        allDensity = uint8([allDensity; t.filteredsynapseDensity]);
        allDANidx = logical([allDANidx; t.DANSynapseIdx]);
        
        allDis = [allDis;  t.filteredDistances];
        allDANSynapseDensity = [allDANSynapseDensity; t.DANsynapseDensity];
        allDis_DANsyn2DANsyn = [allDis_DANsyn2DANsyn; t.Distance_DANSyn2DANSyn];
        allDis_DANsyn2NONDANsyn = [allDis_DANsyn2NONDANsyn; t.Distance_DANSyn2NONDANSyn];
        allDis_NONDANsyn2NONDANsyn = [allDis_NONDANsyn2NONDANsyn; t.Distance_NONDANSyn2NONDANSyn];
        %     end
        %     % combine only the DAN syn idx - balls on + in vicinity
        %      variableInfo = who('-file', [rt2 filesep fn2]);
        %     PF = ismember('filteredData', variableInfo) & ismember('DANSynapseIdx', variableInfo);
        %     if PF
        tt = load([rt2 filesep fn2], 'filteredData','DANSynapseIdx');
        idx = tt.filteredData(:,1)>OL/2 &  tt.filteredData(:,2)>OL/2 & tt.filteredData(:,3)>OL/2 & ...
            tt.filteredData(:,1)<=725 &  tt.filteredData(:,2)<=725 & tt.filteredData(:,3)<=725;
        tt.DANSynapseIdx = tt.DANSynapseIdx(idx);
        allDANsynIDX_includePeriph = [allDANsynIDX_includePeriph; tt.DANSynapseIdx];
    end
    parfor_progress;
end
parfor_progress(0);


%% combine only the DAN syn idx - balls on + in vicinity
rt2 = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/DataStructures/';
allDANsynIDX_includePeriph = [];

parfor_progress(8800)
for i = 1:8800
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_densityData.mat'];
    variableInfo = who('-file', [rt2 filesep fn]);
    PF = ismember('filteredData', variableInfo) & ismember('DANSynapseIdx', variableInfo);
    if PF
        t = load([rt2 filesep fn], 'filteredData','DANSynapseIdx');
        idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
            t.filteredData(:,1)<=725 &  t.filteredData(:,2)<=725 & t.filteredData(:,3)<=725;
        
        t.filteredData = t.filteredData(idx,:);
        allDANsynIDX_includePeriph = [allDANsynIDX_includePeriph; t.filteredData];
    end
    parfor_progress;
end
parfor_progress(0);

%% generate amira point cloud file
% load('20171025_CombinedSynapses_labelidx.mat');
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/';
cd(rt)
load('20171024_CombinedSynapses.mat')
load('20171025_CombinedSynapses_labelidx.mat')
load('20171101_allDensityRecalc.mat')

tic
GU_amiraWriteSpots('/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/WholeBraineSynpaseDensity',allData, allDensity, allDANidx, labelidx,'scales', [0.0237,0.0237,0.044]);
toc

tic
GU_amiraWriteSpots('/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/WholeBraineDANSynpaseDensity',allData(allDANidx,:), allDensity(allDANidx), allDANidx(allDANidx),labelidx(allDANidx), 'scales', [0.0237,0.0237,0.044]);
toc

%% generate masked DAN and nc82
%%%%%% check to see if the compliled file is the same!!! -- paths are
%%%%%% built in
clear ppp
% rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch0\Analysis\clean\slice-tiff\ch0';
ppp{4} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0';
ppp{2} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/Density/slice-tiff/ch0';
ppp{3} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0';
ppp{1} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DANDensityRecalc/slice-tiff/ch0';
k = 4%:numel(ppp)
rt = ppp{k};
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};

for i = 1:numel(fnlist)
    %         fp = [rt filesep fnlist{i}];
    submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "Mask' num2str(i) '" -o ~/logs/Crop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_GenerateMaskedVolumes ' num2str(i-1) '.tif' ' ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
    %     GU_ExM_CropWholeFlyBrainSlice(fp, xmax, ymax)
end

%% calculate ONLY DAN-Associated nc82 puncta -- on each tiled volume 750^3


for i = 1:8800
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "DANass' num2str(i) '" -o ~/logs/DANass' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_ExtractDANassociatednc82 ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% crop stitched z-slices
xmax = 15055;
ymax = 27964;
clear ppp
% rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch0\Analysis\clean\slice-tiff\ch0';
ppp{1} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated/slice-tiff/ch0';
ppp{2} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated/slice-tiff/ch0';
for k = 1:numel(ppp)
    rt = ppp{k};
    fnlist = dir([rt filesep '*.tif']);
    fnlist = {fnlist.name};
    
    for i = 1:numel(fnlist)
        fp = [rt filesep fnlist{i}];
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "Crop' num2str(i) '" -o ~/logs/Crop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_CropWholeFlyBrainSlice ' fp ' ' num2str(xmax) ' ' num2str(ymax) ''''];
        [~, ~] = system(submit,'-echo');
        %     GU_ExM_CropWholeFlyBrainSlice(fp, xmax, ymax)
    end
end

%% generate masked DAN associated and non-associated nc82 data
ppp{1} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated/slice-tiff/ch0';
ppp{2} = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated/slice-tiff/ch0';
for k = 1%:numel(ppp)
    rt = ppp{k};
    fnlist = dir([rt filesep '*.tif']);
    fnlist = {fnlist.name};
    
    for i = 1:numel(fnlist)
        submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "Mask' num2str(i) '" -o ~/logs/Mask' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_GenerateMaskedVolumes ' num2str(i-1) '.tif' ' ' num2str(i) ''''];
        [~, ~] = system(submit,'-echo');
        
    end
end
%% DAN - histogram check
srt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/DANclean_masked_histograms';
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/DANclean_masked';
cd(srt)
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};
parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    %      im = readtiff([rt filesep num2str(i) '.tif']);
    %      vals = uint16(im(im>0));
    %      save([srt filesep num2str(i) '.mat'], 'vals');
    % if ~exist([srt filesep num2str(i-1) '.mat'], 'file')
    GU_extractImageValues(rt, srt, [num2str(i-1) '.tif']);
    % end
    parfor_progress;
end
parfor_progress(0);
%
clear params
parms(1:numel(fnlist)) = struct('BinCounts',[],'BinEdges',[]);
parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    parms(i) = load([srt filesep num2str(i-1) '.mat']);
    parfor_progress;
end
parfor_progress(0);


HC = zeros(1, 2^16, 'single');
parfor_progress(numel(fnlist));
for i = 1:numel(fnlist)
    if parms(i).BinEdges ~= 0
        HC(single(parms(i).BinEdges)) = HC(single(parms(i).BinEdges)) + single(parms(i).BinCounts);
    end
    parfor_progress;
end
parfor_progress(0);

H = 1:2^16;
figure, plot(H, log10(HC+1), 'r');
hold on
for i = 1:16
    line([2^12*(i),2^12*(i)], [0 10])
end

%% check the min max of the DAN channe
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0/';
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};
parfor_progress(numel(fnlist));
imMin = zeros(1, numel(fnlist), 'uint16');
imMax = zeros(1, numel(fnlist), 'uint16');
nLargerthan10k = zeros(1, numel(fnlist), 'uint16');

parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    im = readtiff([rt filesep num2str(i-1) '.tif']);
    imMin(i) = min(im(:));
    imMax(i) = max(im(:));
    nLargerthan10k(i) = sum(im(:)>10000);
    
    parfor_progress;
end
parfor_progress(0);

%% calculate ONLY DAN-Associated nc82 puncta -- on stitched slices

for i = 6:6643
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "DANass' num2str(i) '" -o ~/logs/DANass' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_ExtractDANassociatednc82_stitchedData ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%% DAN - histogram check
srt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/histCounts';
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0';

mkdir(srt)
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};

for i = 1:numel(fnlist)
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "hist' num2str(i) '" -o ~/logs/hist' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_extractImageValues ' rt ' ' srt ' ' fnlist{i} ''''];
    [stat, res] = system(submit,'-echo');
    
end
%     GU_extractImageValues(rt, srt, [num2str(i-1) '.tif']);


%
% nc82 - histogram check
srt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0/histCounts';
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0';

mkdir(srt)
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};

for i = 1:numel(fnlist)
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "hist' num2str(i) '" -o ~/logs/hist' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_extractImageValues ' rt ' ' srt ' ' fnlist{i} ''''];
    [stat, res] = system(submit,'-echo');
    
end

%% generate histograms
% dan channel
srt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/histCounts';
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0';

clear parms
parms = zeros(numel(fnlist),2^16-1, 'uint16');
parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    d = load([srt filesep num2str(i-1) '.mat'], 'BinCounts');
    parms(i,:) = d.BinCounts;
    parfor_progress;
end
parfor_progress(0);
DANBins = parms;
DANhist = sum(DANBins);
DANtotal = sum(DANhist);
DANnorm = DANhist/DANtotal;
DANcdf = cumsum(DANnorm);


srt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0/histCounts';
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0';

clear parms
parms = zeros(numel(fnlist),2^16-1, 'uint16');
parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    d = load([srt filesep num2str(i-1) '.mat'], 'BinCounts');
    parms(i,:) = d.BinCounts;
    parfor_progress;
end
parfor_progress(0);
nc82Bins = parms;
nc82hist = sum(nc82Bins);
nc82total = sum(nc82hist);
nc82norm = nc82hist/nc82total;
nc82cdf = cumsum(nc82norm);

figure, plot(1:2^16-1, DANcdf), hold on;pause;  plot(1:2^16-1, nc82cdf)

%% re-generate masked data sets

rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0';
fnlist = dir([rt filesep '*.tif']);
fnlist = {fnlist.name};

for i = 1:numel(fnlist)
    submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "Mask' num2str(i) '" -o ~/logs/Mask' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_GenerateMaskedVolumes ' num2str(i-1) '.tif' ' ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
    
end

%% anaylyze v5 and nc82
parfor i = 1:2
    GU_ExM_calcPunctaDensityVol_1ch(data.framePaths{i}{1}, ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [246],'ExpansionFactor', 3.99)
end
%%
mask = readtiff('D:\Gokul\Wholebrainsynapse\Ex15_Calyxpart1_V5+nc82_z0p18\Mask\nc82-a3-BrpV5 samples_calyx-mask.tif');
idx = sub2ind(size(mask), v5(:,2),v5(:,1),v5(:,3));
maskidx_v5 = mask(idx);
idx = sub2ind(size(mask), nc82(:,2),nc82(:,1),nc82(:,3));
maskidx_nc82 = mask(idx);

ef = 3.99;
px = 0.097;
epx = px/ef;
epz = 0.18/ef;
%% combine masks
rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/mask2/';
d = dir(rt);

fn = dir([rt d(3).name filesep '*.tif']);
fn = fn.name;
imt = readtiff([rt d(3).name filesep fn]);
maskComb = zeros(size(imt), 'uint8');
clear imt fn


for i = 3:numel(d)
    fn = dir([rt d(i).name filesep '*.tif']);
    fn = fn.name;
    imt = logical(readtiff([rt d(i).name filesep fn]));
    maskComb(imt) = i-2;
    i
end
writetiff(maskComb, [rt 'CombinedMask_37labels_fixedv2.tif']);

%% janelia calculate surface area and volume

for i = 1:8800
    submit = ['bsub -n 5 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SAV' num2str(i) '" -o ~/logs/SAV' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_WholeBrainCalcSurfaceAreaVolume ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%%
% combine all the data - surface area and volume

% get file names

mx = 15055;
my = 27964;
mz = 6647;
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
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            nn = nn + 1;
        end
    end
end

rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/dataStructures';

skelCoordXYZ = cell(1, 8800);
iSAcoordXYZ=cell(1, 8800);
VolcoordXYZ = cell(1, 8800);
zA = 0.18/0.097;
parfor_progress(8800);
parfor i = 1:8800
    %     fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '.mat'];
    %     t = load([rt filesep fn]);
    %
    %     if ~isempty(t.iVolcoordXYZ)
    %         idx = t.iVolcoordXYZ(:,1)>OL/2 &  t.iVolcoordXYZ(:,2)>OL/2 & t.iVolcoordXYZ(:,3)>OL/2*zA & ...
    %             t.iVolcoordXYZ(:,1)<=725 &  t.iVolcoordXYZ(:,2)<=725 & t.iVolcoordXYZ(:,3)<=725*zA;
    %         t.iVolcoordXYZ = t.iVolcoordXYZ(idx,:);
    %
    %         idx = t.SAcoordXYZ(:,1)>OL/2 &  t.SAcoordXYZ(:,2)>OL/2 & t.SAcoordXYZ(:,3)>OL/2*zA & ...
    %             t.SAcoordXYZ(:,1)<=725 &  t.SAcoordXYZ(:,2)<=725 & t.SAcoordXYZ(:,3)<=725*zA;
    %         t.SAcoordXYZ = t.SAcoordXYZ(idx,:);
    %
    %         t.iVolcoordXYZ(:,1) = t.iVolcoordXYZ(:,1) + xmin(i)-1;
    %         t.iVolcoordXYZ(:,2) = t.iVolcoordXYZ(:,2) + ymin(i)-1;
    %         t.iVolcoordXYZ(:,3) = t.iVolcoordXYZ(:,3) + (zmin(i)-1)*zA;
    %
    %         t.SAcoordXYZ(:,1) = t.SAcoordXYZ(:,1) + xmin(i)-1;
    %         t.SAcoordXYZ(:,2) = t.SAcoordXYZ(:,2) + ymin(i)-1;
    %         t.SAcoordXYZ(:,3) = t.SAcoordXYZ(:,3) + (zmin(i)-1)*zA;
    %
    %         iVolcoordXYZ{i} = t.iVolcoordXYZ;
    %         iSAcoordXYZ{i} = t.SAcoordXYZ;
    %     end
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '.mat'];
    t = load([rt filesep fn], 'VolcoordXYZ');
    
    if ~isempty(t.VolcoordXYZ)
        idx = t.VolcoordXYZ(:,1)>OL/2 &  t.VolcoordXYZ(:,2)>OL/2 & t.VolcoordXYZ(:,3)>OL/2*zA & ...
            t.VolcoordXYZ(:,1)<=725 &  t.VolcoordXYZ(:,2)<=725 & t.VolcoordXYZ(:,3)<=725*zA;
        t.VolcoordXYZ = t.VolcoordXYZ(idx,:);
        
        
        t.VolcoordXYZ(:,1) = t.VolcoordXYZ(:,1) + xmin(i)-1;
        t.VolcoordXYZ(:,2) = t.VolcoordXYZ(:,2) + ymin(i)-1;
        t.VolcoordXYZ(:,3) = t.VolcoordXYZ(:,3) + (zmin(i)-1)*zA;
        
        
        
        VolcoordXYZ{i} = uint16(t.VolcoordXYZ);
        
    end
    parfor_progress;
end
parfor_progress(0);

% alliSACoordXYZ = uint16(vertcat(iSAcoordXYZ{:}));
% alliVolCoordXYZ = uint16(vertcat(iVolcoordXYZ{:}));
allnonInterpVolCoordXYZ = uint16(vertcat(VolcoordXYZ{:}));

%% generate point cloud

load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/20180214_interpSurfAreaCoord_uint16.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/SurfaceAreaLabels.mat')

%%

idx = SA_labelidx>=17 & SA_labelidx<=22;
s = find(idx);
idx = s(1:10:end);
tic
GU_amiraWriteSpots('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/MBs.am',alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
toc

idx = SA_labelidx>=27 & SA_labelidx<=28;
s = find(idx);
idx = s(1:10:end);
tic
GU_amiraWriteSpots('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/ATLs.am',alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
toc

idx = SA_labelidx==29;
s = find(idx);
idx = s(1:10:end);
tic
GU_amiraWriteSpots('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/PBs.am',alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
toc

idx = SA_labelidx==30;
s = find(idx);
idx = s(1:10:end);
tic
GU_amiraWriteSpots('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/EBs.am',alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
toc

idx = SA_labelidx==31;
s = find(idx);
idx = s(1:10:end);
tic
GU_amiraWriteSpots('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/FBs.am',alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
toc
%%
n = [1:16 24:26,32:37];
for i = n
    idx = SA_labelidx==i;
    s = find(idx);
    idx = s(1:10:end);
    tic
    GU_amiraWriteSpots(['/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/' num2str(i) '.am'],alliSACoordXYZ(idx,:), ones(numel(idx),1), [], SA_labelidx(idx),'scales', [0.0237,0.0237,0.0237]);
    toc
end

%% janelia calculate skeleton

for i = 1:8800
    submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Sk' num2str(i) '" -o ~/logs/Sk' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_WholeBrainCalcSkeletonDistanceMap ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end

%%
% combine all the data - surface area and volume

% get file names

mx = 15055;
my = 27964;
mz = 6647;
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
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            nn = nn + 1;
        end
    end
end

rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/dataStructures/skel/';

SkelcoordXYZ = cell(1, 8800);
SkelRadius=cell(1, 8800);

parfor_progress(8800);
parfor i = 1:8800
    
    fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '.mat'];
    t = load([rt filesep fn], 'SkelcoordXYZ', 'SkelRadius');
    
    if ~isempty(t.SkelcoordXYZ)
        idx = t.SkelcoordXYZ(:,1)>OL/2 &  t.SkelcoordXYZ(:,2)>OL/2 & t.SkelcoordXYZ(:,3)>OL/2 & ...
            t.SkelcoordXYZ(:,1)<=725 &  t.SkelcoordXYZ(:,2)<=725 & t.SkelcoordXYZ(:,3)<=725;
        t.SkelcoordXYZ = t.SkelcoordXYZ(idx,:);
        t.SkelcoordXYZ(:,1) = t.SkelcoordXYZ(:,1) + xmin(i)-1;
        t.SkelcoordXYZ(:,2) = t.SkelcoordXYZ(:,2) + ymin(i)-1;
        t.SkelcoordXYZ(:,3) = t.SkelcoordXYZ(:,3) + (zmin(i)-1);
        SkelcoordXYZ{i} = uint16(t.SkelcoordXYZ);
        SkelRadius{i} = t.SkelRadius;
    end
    parfor_progress;
end
parfor_progress(0);


allSkelcoordXYZ = uint16(vertcat(SkelcoordXYZ{:}));
allSkelRadii = (vertcat(SkelRadius{:}));

%% label skel coordinates
cd('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain')
mask = uint8(readtiff('CombinedMask_37labels_fixedv2.tif'));

mx = 15055;
my = 27964;
mz = 6647;

% mask dimensions
sy = 1747;
sx = 940;
sz = 1661;

% mask downsampled
fy = my/sy; %16.0155;%
fx = mx/sx; %16.032;
fz = mz/sz;

x = double(allSkelcoordXYZ(:,1)/fx);
y = double(allSkelcoordXYZ(:,2)/fy);
z = double(allSkelcoordXYZ(:,3)/fz);
[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
labelidxRadii = interp3(X,Y,Z,mask,x,y,z,'nearest');
toc

%% back
% GU_amiraWriteSpots(filename,xyzCoord, Density, DensityIdx, labelidx, varargin)
n = 1:37;
for i = n
    idx = labelidxRadii==i;
    tic
    GU_amiraWriteSpots(['/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/skelAmiraPts/' num2str(i) '.am'],allSkelcoordXYZ(idx,:), allSkelRadii(idx), [], labelidxRadii(idx),'scales', [0.0237,0.0237,0.044]);
    toc
end

%% write
load('20180320_Combined_SkeletionizedFlyBrainCoordinates.mat')
tic
GU_amiraWriteSpots(['/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/skelAmiraPts/All_zCoordisDensity.am'],allSkelcoordXYZ, allSkelRadii, [], allSkelcoordXYZ(:,3),'scales', [0.0237,0.0237,0.044]);
toc

%% analysis below is to address reviewer comments
% there are artifacts in figure 6E, possibly from Arivis rendering
% generate cropped tifs
    % processing this section from ww8 workstation, as nearline is not accessible from the janelia cluster
    
clear rt
rt{1} = 'Z:\RG\ExLLSM\4Stephan-04-2018\171016_Flybrainsynapse\171016_100nmcluster\ch0\DANclean_masked_rescaled_v3\';
rt{2} = 'Z:\RG\ExLLSM\4Stephan-04-2018\171016_Flybrainsynapse\171016_100nmcluster\ch1\nc82DANassociated_masked_rescaled_v3\';
srt{1} = [rt{1} 'crop' filesep];
srt{2} = [rt{2} 'crop' filesep];
mkdir(srt{1});
mkdir(srt{2});

% crop coordinates
xrmin = 3232;
xrmax = 3232+9280;
yrmin = 640;
yrmax = 640+5280;
z = 0:6646;
for k = 1:numel(rt)
    parfor_progress(numel(z));
    parfor i = z
        fn = [rt{k} num2str(i) '.tif'];
         imCrop = GU_CropTiffSlices({fn}, xrmin, xrmax, yrmin, yrmax);
%         submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
%         [stat, res] = system(submit,'-echo');
        writetiff(imCrop,[srt{k}  num2str(i) '_crop.tif']);
        parfor_progress;
    end

    parfor_progress(0);
end


% % %% downsample fly brain ; generated a perspective downsampled z-stack
% % danrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3\';
% % danbinrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3_binned\';
% % mkdir(danbinrt);
% % fn = dir([danrt '*.tif']);
% % fn = {fn.name};
% %
% % parfor_progress(numel(fn));
% % parfor i = 1:numel(fn)
% %     im = readtiff([danrt fn{i}]);
% %     J = imresize(im, 0.25);
% %     writetiff(J, [danbinrt fn{i}], 'Compression', 'lzw');
% %     parfor_progress;
% % end
% % parfor_progress(0);
% %
% % %
% % % stitched data dimensions
% % danrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3\';
% % danbinrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3_binned\';
% % danbinMRrt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3_binned_maskedRegions\';
% % for i = 1:14
% %     mkdir([danbinMRrt num2str(i)])
% % end
% %
% % danbinMRrt_perspective = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\DANclean_masked_rescaled_v3_binned_maskedRegionsPerspective\';
% % for i = 1:14
% %     mkdir([danbinMRrt_perspective num2str(i)])
% % end
% %
% % rt = danbinrt; %danrt;%
% % fn = dir([rt '*.tif']);
% % fn = {fn.name};
% %
% % % read mask
% % % maskpath = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\mask\CombinedMask.tif';
% % % fullmask = readtiff(maskpath);
% % maskFrames = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\171016_Flybrainsynapse\mask\frames\';
% % % mkdir(maskFrames);
% % %
% % % for i = 1:size(fullmask,3)
% % %     writetiff(fullmask(:,:,i), [maskFrames num2str(i) '.tif'], 'Compression', 'none')
% % % end
% % dsf = 4; % downsample factor
% %
% % mx = 15055/dsf;
% % my = 27964/dsf;
% % mz = 6647;
% %
% % % mask dimensions
% % sx = 940;
% % sy = 1747;
% % sz = 1661;
% %
% % % mask downsampled
% % fx = 16.032/dsf;%my/sy;
% % fy = 16.0155/dsf;%my/sy;
% % fz = mz/sz;
% %
% % MPR = 0.3;% max perspective reduction over depth
% % z = 6646:-1:1;
% % parfor_progress(numel(z));
% % PRS = MPR/numel(z); % perspective reduction step size
% % parfor i = 1:numel(z)
% %     n = z(i);
% % %     if ~exist([danbinMRrt_perspective filesep num2str(14) filesep  num2str(n) '.tif'], 'file')
% %         im = readtiff([rt num2str(n) '.tif']);
% %         maskslice = readtiff([maskFrames num2str(ceil(n/fz)) '.tif']);
% %         so = size(im);
% %         mask = imresize(maskslice, so, 'nearest');
% %
% %         for j = 1:14
% %             mim = zeros(size(im), 'uint16');
% %             mim(mask==j) = im(mask==j);
% % %             writetiff(mim, [danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'lzw');
% %             pim = zeros(so, 'uint16');
% %             mim = imresize(mim, 1-PRS*i);
% %             sr = size(mim);
% %             dr = so-sr;
% %             dr2 = round((dr)/2);
% %             pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
% %             writetiff(pim, [danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'lzw');
% %         end
% % %     end
% %     parfor_progress;
% % end
% % parfor_progress(0);
% %
% %
% % %% check for missing files
% % missingF = zeros(numel(z),14);
% % missingF_P = zeros(numel(z),14);
% % parfor_progress(numel(z));
% % parfor i = 1:numel(z)
% %     n = z(i);
% %     for j = 1:14
% %         if ~exist([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif'], 'file')
% %             missingF(i,j) = 1;
% %         end
% %
% %         if ~exist([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'file')
% %             missingF_P(i,j) = 1;
% %         end
% %     end
% %     parfor_progress;
% % end
% % parfor_progress(0);
% %
% %
% %
% % %% make max projection of the data
% % parfor_progress(numel(z)*14);
% % % qqq = [11];
% % parfor j = 1:14
% %     n = 1;
% % %     j = qqq(nnn);
% %     maxPIM = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
% %     maxPIM(maxPIM>2^12*j) = 0;
% %     maxIM = readtiff([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif']);
% %     maxIM(maxIM>2^12*j) = 0;
% %     for i = 2:numel(z)
% %         n = z(i);
% %         imtmp = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
% %         imtmp(imtmp>2^12*j) = 0;
% %         maxPIM = max(maxPIM, imtmp);
% %         imtmp = readtiff([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif']);
% %          imtmp(imtmp>2^12*j) = 0;
% %         maxIM = max(maxIM, imtmp);
% %         parfor_progress;
% %     end
% %     writetiff(maxPIM, [danbinMRrt_perspective filesep  num2str(j) '_max.tif'])
% %     writetiff(maxIM, [danbinMRrt filesep  num2str(j) '_max.tif'])
% % end
% % parfor_progress(0);
%% downsample fly brain ; generated a perspective downsampled z-stack
danrt = 'G:\GokulWorkspace\171016_Flybrainsynapse\DAN\';
danbinrt = 'G:\GokulWorkspace\171016_Flybrainsynapse\DAN_binned\';
mkdir(danbinrt);

nc82rt = 'G:\GokulWorkspace\171016_Flybrainsynapse\DANAssociated_nc82_rad10\';
nc82binrt = 'G:\GokulWorkspace\171016_Flybrainsynapse\DANAssociated_nc82_rad10_binned\';
mkdir(nc82binrt);

fn = dir([danrt '*.tif']);
fn = {fn.name};
nofiles = zeros(numel(fn),1);
parfor_progress(numel(fn));
parfor i = 1:numel(fn)
    %read dan signals
    im = readtiff([danrt fn{i}]);
    J = imresize(im, 0.25);
    writetiff(J, [danbinrt fn{i}], 'Compression', 'none');
    
    %read nc82 signals
    if exist([nc82rt fn{i}], 'file')
        im = readtiff([nc82rt fn{i}]);
        J = imresize(im, 0.25);
    else
        nofiles(i) =1;
        J = zeros(size(J),'uint16');
    end
    writetiff(J, [nc82binrt fn{i}], 'Compression', 'none');
    
    parfor_progress;
end
parfor_progress(0);


% stitched data dimensions
danbinMRrt_perspective = 'G:\GokulWorkspace\171016_Flybrainsynapse\DAN_binned_maskedRegionsPerspective\';
for i = 1:14
    mkdir([danbinMRrt_perspective num2str(i)])
end

nc82binMRrt_perspective = 'G:\GokulWorkspace\171016_Flybrainsynapse\DANAssociated_nc82_rad10_binned_maskedRegionsPerspective\';
for i = 1:14
    mkdir([nc82binMRrt_perspective num2str(i)])
end

maskFrames = 'G:\GokulWorkspace\171016_Flybrainsynapse\frames\';
dsf = 4; % downsample factor

mx = 15055/dsf;
my = 27964/dsf;
mz = 6647;

% mask dimensions
sx = 940;
sy = 1747;
sz = 1661;

% mask downsampled
fx = 16.032/dsf;%my/sy;
fy = 16.0155/dsf;%my/sy;
fz = mz/sz;

MPR = 0.4;% max perspective reduction over depth
z = 6646:-1:1;
parfor_progress(numel(z));
PRS = MPR/numel(z); % perspective reduction step size
parfor i = 1:numel(z)
    n = z(i);
    %     im = readtiff([danbinrt num2str(n) '.tif']); %dan
    im2 = readtiff([nc82binrt num2str(n) '.tif']); %nc82
    maskslice = readtiff([maskFrames num2str(ceil(n/fz)) '.tif']);
    so = size(im2);
    mask = imresize(maskslice, so, 'nearest');
    
    for j = 1:14
        
        %         mim = zeros(size(im), 'uint16');
        %         mim(mask==j) = im(mask==j);
        %         mim = imresize(mim, 1-PRS*i);
        %         pim = zeros(so, 'uint16');
        %         sr = size(mim);
        %         dr = so-sr;
        %         dr2 = round((dr)/2);
        %         pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
        %         writetiff(pim, [danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'none');
        % nc82
        mim = zeros(size(im2), 'uint16');
        mim(mask==j) = im2(mask==j);
        mim = imresize(mim, 1-PRS*i);
        pim = zeros(so, 'uint16');
        sr = size(mim);
        dr = so-sr;
        dr2 = round((dr)/2);
        pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
        writetiff(pim, [nc82binMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'none');
    end
    
    parfor_progress;
end
parfor_progress(0);

% make max projection of the data
parfor_progress(numel(z)*14);
% % qqq = [11];
parfor j = 1:14
    n = 1;
    % %     j = qqq(nnn);
    %     maxPIM = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
    %     maxPIM(maxPIM>2^12*j) = 0;
    maxIM = readtiff([nc82binMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
    %     maxIM(maxIM>2^12*j) = 0;
    for i = 2:numel(z)
        n = z(i);
        %         imtmp = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
        % %         imtmp(imtmp>2^12*j) = 0;
        %         maxPIM = max(maxPIM, imtmp);
        imtmp = readtiff([nc82binMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
        %          imtmp(imtmp>2^12*j) = 0;
        maxIM = max(maxIM, imtmp);
        parfor_progress;
    end
    %     writetiff(maxPIM, [danbinMRrt_perspective filesep  num2str(j) '_max.tif'], 'Compression', 'none');
    writetiff(maxIM, [nc82binMRrt_perspective filesep  num2str(j) '_max.tif'], 'Compression', 'none');
end
parfor_progress(0);
%% check for missing files
missingF = zeros(numel(z),14);
missingF_P = zeros(numel(z),14);
parfor_progress(numel(z));
parfor i = 1:numel(z)
    n = z(i);
    for j = 1:14
        if ~exist([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif'], 'file')
            missingF(i,j) = 1;
        end
        
        if ~exist([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'file')
            missingF_P(i,j) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);



% %% make max projection of the data
% parfor_progress(numel(z)*14);
% % qqq = [11];
% parfor j = 1:14
%     n = 1;
% %     j = qqq(nnn);
%     maxPIM = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
%     maxPIM(maxPIM>2^12*j) = 0;
%     maxIM = readtiff([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif']);
%     maxIM(maxIM>2^12*j) = 0;
%     for i = 2:numel(z)
%         n = z(i);
%         imtmp = readtiff([danbinMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif']);
%         imtmp(imtmp>2^12*j) = 0;
%         maxPIM = max(maxPIM, imtmp);
%         imtmp = readtiff([danbinMRrt filesep num2str(j) filesep  num2str(n) '.tif']);
%          imtmp(imtmp>2^12*j) = 0;
%         maxIM = max(maxIM, imtmp);
%         parfor_progress;
%     end
%     writetiff(maxPIM, [danbinMRrt_perspective filesep  num2str(j) '_max.tif'])
%     writetiff(maxIM, [danbinMRrt filesep  num2str(j) '_max.tif'])
% end
% parfor_progress(0);
%% rescale custom perspective mips
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/Figures_ExLLSM/Figure5+6_Fly brain/Figure 6_Synapse/customPerspectiveMIPs_25percentRes/rawPersMIPs/';
srt = '/Users/GU/Dropbox/Manuscript_ExLLSM/Figures_ExLLSM/Figure5+6_Fly brain/Figure 6_Synapse/customPerspectiveMIPs_25percentRes/rescaledPresMIPs/';
cd(rt)
fn = dir([rt '*.tif']);
fn = {fn.name}';

for i = 1:numel(fn)
    im = readtiff([rt fn{i}]);
    %     rangeIn = [2^12*(i-1), (2^12*i)-1];
    rangeIn = [0, (2^12*i)-1];
    rangeOut = [0, 2^16];
    im = uint16(scaleContrast(im, rangeIn, rangeOut));
    writetiff(im, [srt fn{i}]);
end

%% run region separation for full res perspective mip on Janelia cluster

for i = 5:6641
    
    %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "mippers' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_ExtractDANRegions ' num2str(i) ''''];
    [stat, res] = system(submit,'-echo');
end
%% generate max projections

rt{1} = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/';
rt{2} = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/nc82_DANAssociated_masked_rad7_maskedPers/';
n = 4:50:6641;

for k = 1:numel(rt)
    for j = 1:15
        srt = [rt{k} num2str(j) 'maxproj' filesep];
        mkdir(srt);
        for i = 1:numel(n)-1
            minfn = n(i)+1;
            maxfn = n(i+1);
            sfn = [num2str(i) '.tif'];
            inrt = [rt{k} num2str(j) filesep];
            % GU_ExM_JaneliaCluster_MaxProject(num2str(minfn), num2str(maxfn), inrt, srt, sfn);
            submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "mippers' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MaxProject ' num2str(minfn) ' ' num2str(maxfn) ' ' inrt ' ' srt ' ' sfn ''''];
            [stat, res] = system(submit,'-echo');
        end
    end
end

%% round two of max proj

rt{1} = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/';
rt{2} = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/nc82_DANAssociated_masked_rad7_maskedPers/';

minfn = 1;
maxfn = numel(n)-1;
for k = 1:numel(rt)
    for i = 1:15
        srt = [rt{k}];
        
        sfn = [num2str(i) '_maxFullRes.tif'];
        inrt = [rt{k} num2str(i) 'maxproj' filesep];
        % GU_ExM_JaneliaCluster_MaxProject(num2str(minfn), num2str(maxfn), inrt, srt, sfn);
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "mippers' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MaxProject ' num2str(minfn) ' ' num2str(maxfn) ' ' inrt ' ' srt ' ' sfn ''''];
        [stat, res] = system(submit,'-echo');
        
    end
end


%% crop out EB

clear rt
rt{1} = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/1/';

xrmin = 5670;
xrmax = xrmin+2142;
yrmin = 13400;
yrmax = yrmin+2217;
z = 5:6641;
for k = 1:numel(rt)
    for i = z
        fn = [rt{k} num2str(i) '.tif'];
        %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

% check when slices are occupied by EB
rt = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/1_x_5670_7812_y_13400_15617/';
idxO = zeros(1,numel(z));
parfor_progress(numel(z));
parfor i = z
    im = readtiff([rt num2str(i) '.tif']);
    idxO(i) = sum(im(:));
    parfor_progress;
end
parfor_progress(0);

idxEB = find(idxO);


% make pespective EB resizing
MPR = 0.7;% max perspective reduction over depth

rt_pers = ['/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/EB_perspective/' num2str(MPR) '/'];
mkdir(rt_pers);

parfor_progress(numel(idxEB));
PRS = MPR/numel(idxEB); % perspective reduction step size

parfor i = 1:numel(idxEB)
    n = idxEB(i);
    im2 = readtiff([rt num2str(n) '.tif']); %nc82
    so = size(im2);
    mim = imresize(im2, 1-PRS*(i-1));
    pim = zeros(so, 'uint16');
    sr = size(mim);
    dr = so-sr;
    dr2 = round((dr)/2);
    pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
    writetiff(pim, [rt_pers filesep  num2str(n) '.tif'], 'Compression', 'none');
    
    parfor_progress;
end
parfor_progress(0);


%% process soma cleanup to check for dan signals
%% processing of stereotypy data - for paper revisions
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/Flybrain_forRevision/SomaNc82visualization/Raw/';
% rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/Flybrain_forRevision/SomaNc82visualization/Raw/';
cd(rt)
data = GU_loadConditionData;
%%
parfor i = 1:numel(data)
    fn1 = data(i).framePaths{1};
    fn2 = data(i).framePaths{2};
    GU_ExM_calcPunctaDensityVol(fn2, fn1, ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192, 192],'ExpansionFactor', 4.09, 'MinThreshold', [550,720])
    %     GU_ExM_calcSynapseDensityVol(fn2, fn1,'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);
end


%% mb data for stereotypy
rt{1} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\F1_L\Stitch_Igor\Stitch_decon_blend\slice-tiff\ch0\';
rt{2} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\F1_R\Stitch_Igor\Stitch_decon_blend\slice-tiff\ch0\';
rt{3} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\F2_L\Stitch_Igor\Stitch_decon_blend\slice-tiff\ch0\';
rt{4} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\F3_L\Stitch_Igor\Stitch_decon_blend\slice-tiff\ch0\';
rt{5} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\F3_R\Stitch_Igor\Stitch_decon_blend\slice-tiff\ch0\';
rt{6} = '\\dm11\betziglab\4Stephan\\181009_MBalpha3\RBPF1_R\Stitch_Igor\Stitch_decon_blend\slice-tiff_stage\ch1\';

srt{1} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex1_F1L_z0p18\';
srt{2} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex2_F1R_z0p18\';
srt{3} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex3_F2L_z0p18\';
srt{4} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex4_F3L_z0p18\';
srt{5} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex5_F3R_z0p18\';
srt{6} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex6_RBPF1R_z0p18\';

for i = 1:numel(rt)
    fn = dir([rt{i} '*.tif']);
    fn = {fn.name};
    tmp = readtiff([rt{i} fn{1}]);
    s = size(tmp);
    im = zeros(s(1), s(2), numel(fn), 'uint16');
    im(:,:,1) = tmp;
    parfor_progress(numel(fn));
    tic
    parfor j = 2:numel(fn)
        im(:,:,j) = readtiff([rt{i} fn{j}]);
        parfor_progress;
    end
    parfor_progress(0);
    toc
    %     tic
    %     writetiff(im, [srt{i} 'nc82.tif']); toc
    %     clear im tmp fn s
    tic 
        imtmp = im(:,:,10:50:end);
        PP = prctile(imtmp(:),99.9);
            T(1) = thresholdOtsu(imtmp(imtmp > 0 & imtmp < PP));
    
    fn = [srt{i},'part1.tif'];
    
    GU_ExM_calcPunctaDensityVol_1chVol(fn, im(:,:,1:1000), ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192],'ExpansionFactor', 3.99, 'Threshold', T,...
        'PixelSize', 0.102,'PixelSizeZ', 0.18)
    
    fn = [srt{i},'part2.tif'];
    
    GU_ExM_calcPunctaDensityVol_1chVol(fn, im(:,:,950:2000), ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192],'ExpansionFactor', 3.99, 'Threshold', T,...
        'PixelSize', 0.102,'PixelSizeZ', 0.18)
    
    fn = [srt{i},'part3.tif'];
    
    GU_ExM_calcPunctaDensityVol_1chVol(fn, im(:,:,1950:3000), ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192],'ExpansionFactor', 3.99, 'Threshold', T,...
        'PixelSize', 0.102,'PixelSizeZ', 0.18)
    
    fn = [srt{i},'part4.tif'];
    
    GU_ExM_calcPunctaDensityVol_1chVol(fn, im(:,:,2950:end), ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192],'ExpansionFactor', 3.99, 'Threshold', T,...
        'PixelSize', 0.102,'PixelSizeZ', 0.18)
end

%% extract and merge local maxia
srt{1} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex1_F1L_z0p18\';
srt{2} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex2_F1R_z0p18\';
srt{3} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex3_F2L_z0p18\';
srt{4} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex4_F3L_z0p18\';
srt{5} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex5_F3R_z0p18\';
srt{6} = '\\dm11\betziglab\4Stephan\181009_MBalpha3\GokulWorkspace\Ex6_RBPF1R_z0p18\';
suf = 'Analysis_1ch\192\DataStructures\';
zoff = [25];
zr = [0, 950, 1950,2950];
for i = 1:numel(srt)
    if i==6
        k=4;
    else
         k=3;
    end
    allFilteredData = [];
    for j = 1:k
        fd = load([srt{i} suf 'part' num2str(j) '_192_densityData.mat'], 'filteredData');
        filteredData = fd.filteredData;
        idx = filteredData(:,3) > zoff & filteredData(:,3) <= zr(2)+zoff;
        filteredData = filteredData(idx,:);
        filteredData(:,3) = filteredData(:,3) + zr(j);
        allFilteredData = [allFilteredData; filteredData];
    end
    save([srt{i} 'MergedFilteredData' num2str(i) '.mat'], 'allFilteredData')
    i
end

