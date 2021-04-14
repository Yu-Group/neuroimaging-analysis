%% mouse visual cortex workspace

% rt = '/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/';
s = GU_calcGaussianIntegral(1, [3.2 3.2 3.2]);

GU_ExM_calcPunctaDensityVol(data.framePaths{1}{3}, data.framePaths{1}{1},...
    'FindLocalMaxima', false, 'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [round(s), 1008],'ExpansionFactor', 3.62); % sigma 3.9
% T = 1400 for casper; T = 1100 for yfp; T = 200 for myelin

GU_ExM_calcPunctaDensityVol_1ch(data.framePaths{1}{2}, ...
    'FindLocalMaxima', false, 'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [1008],'ExpansionFactor', 3.62)



%% submit for sub-vol extraction,
rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/';
mkdir([rt 'ch0'])
mkdir([rt 'ch1'])
mkdir([rt 'ch2'])

mx = 10727;
my = 39519;
mz = 5413;
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
    submit = ['bsub -n 6 -R"affinity[core(1)]" -J' ' "VC' num2str(i) '" -o ~/logs/VC' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MouseVizCortex_Chromatic_cleanup ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end
%%
rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1/Analysis/516/';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};
for ex = 1:nn-1
    fn = [filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_516_clean.tif'];
    CompletedJobsIdx(ex) = ismember(fn(2:end),fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);
for i = job2resubmit
    submit = ['bsub -n 3 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "VC' num2str(i) '" -o ~/logs/VC' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_MouseVizCortex_Chromatic_cleanup ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end

%%
mx = 10727;
my = 39519;
mz = 5413;
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
rt{1} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0/Analysis/1008/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1/Analysis_1ch/1008/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch2/Analysis/516/';
a = 0;
for k = 1:3
    for i = 1:nn-1
        a = a +1;
        if k == 1
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_clean.tif'];
        end
        if k==2
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_1008_clean.tif'];
        end
        if k ==3
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean.tif'];
        end
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Rot' num2str(a) '" -o ~/logs/Rot' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end


%% find failed jobs
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0/Analysis/1008/RotatedStacks/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1/Analysis_1ch/1008/RotatedStacks/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch2/Analysis/516/RotatedStacks/';
for k = 1:3
    fnlist = dir([rt{k} '*.tif']);
    fnlist = {fnlist.name};
    CompletedJobsIdx = cell(1, 3);
    for ex = 1:nn-1
        if k == 1
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_clean.tif'];
        end
        if k==2
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_1008_clean.tif'];
        end
        if k ==3
            fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_516_clean.tif'];
        end
        CompletedJobsIdx{k}(ex) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
    end
    job2resubmit{k} = find(CompletedJobsIdx{k}==0);
end

%% change coordinates
mx = 10727;
my = 39519;
mz = 5413;
zxRatio = 0.18/0.097;
theta = 32;
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0/Analysis/1008/RotatedStacks/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1/Analysis_1ch/1008/RotatedStacks/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch2/Analysis/516/RotatedStacks/';
fileSufix = {'_clean', '_1008_clean', '_516_clean'};
for k = 1:3
    GU_ExM_WriteRotatedCoordinates(rt{k}, mx, my, mz, zxRatio, theta, fileSufix{k})
end


%% DAN - histogram check
srt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0/Analysis/1008/RotatedStacks/params';

% cd(srt)
fnlist = dir([srt filesep '*.mat']);
fnlist = {fnlist.name};

%
clear params
params(1:numel(fnlist)) = struct('BinCounts',[],'edges',[]);
parfor_progress(numel(fnlist));
parfor i = 1:numel(fnlist)
    params(i) = load([srt filesep fnlist{i}], 'BinCounts', 'edges');
    parfor_progress;
end
parfor_progress(0);


HC = zeros(1, 2^16, 'single');
parfor_progress(numel(fnlist));
for i = 1:numel(fnlist)
    if params(i).edges ~= 0
        HC(single(params(i).edges)) = HC(single(params(i).edges)) + single(params(i).BinCounts);
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

%% submit masking with yfp and spectral mixing
df = 9;
for i = 3000:6000
    submit = ['bsub -n 2 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Viz' num2str(i) '" -o ~/logs/Viz' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_Janelia_MouseVisualCortex_Mask_SpectralMix ' num2str(i) ' ' num2str(df)  ''''];
    [stat, res] = system(submit,'-echo');
end

%% resubmit missing jobs
clear rt
rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch1/YFPMaskedSpectralMixed/';
df = 7;
for i = 3000:6000
    if ~exist([rt num2str(i) '.tif'], 'file')
        submit = ['bsub -n 3 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "Viz' num2str(i) '" -o ~/logs/Viz' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_Janelia_MouseVisualCortex_Mask_SpectralMix ' num2str(i) ' ' num2str(df)  ''''];
        [stat, res] = system(submit,'-echo');
        %         FDE(i) = 1;
    end
    
end
%% generate cropped tifs
rt{1} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch0/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch1/YFPMaskedSpectralMixed/';
rt{3} = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch2/YFPMaskedSpectralMixed/';
xrmin = 4000;
xrmax = 5000;
yrmin = 11000;
yrmax = 12000;
for k = 1:numel(rt)
    for i = 4500:5250
        fn = [rt{k} num2str(i) '.tif'];
        %         GU_CropTiffSlices({fn}, xrmin, xymax, yrmin, yrmax);
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SliceCrop' num2str(i) '" -o ~/logs/SliceCrop' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_CropTiffSlices ' fn ' ' num2str(xrmin) ' ' num2str(xrmax) ' ' num2str(yrmin) ' ' num2str(yrmax)  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

%%
% generate figures for the cover
% generate perective rescale of each of the files
rt = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/Position2_tile/Stitch_Igor/Converted/slice-tiff_fusion/ch0/';
rt_perspective = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/Position2_tile/Stitch_Igor/Converted/slice-tiff_fusion/ch0/pers';
mkdir(rt_perspective);

MPR = 0.4;% max perspective reduction over depth
z = 840:-1:0;
parfor_progress(numel(z));
PRS = MPR/numel(z); % perspective reduction step size
parfor i = 1:numel(z)
    n = z(i);
    im2 = readtiff([rt num2str(n) '.tif']);
    so = size(im2);
    mim = imresize(im2, 1-PRS*i);
    pim = zeros(so, 'uint16');
    sr = size(mim);
    dr = so-sr;
    dr2 = round((dr)/2);
    pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
    writetiff(pim, [rt_perspective filesep  num2str(n) '.tif'], 'Compression', 'none');
    
    parfor_progress;
end
parfor_progress(0);

%% spine detection

% read the distance transform valued skeleton
rt = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/ch0_Analysis_1ch/1008/DistMapSkel/';
s = [4113, 3170, 700];
im = zeros(s, 'uint16');

fn = dir([rt '*.tif']);
fn = {fn.name}';
for i = 1:7
    zr = (100*(i-1)+1):100*i;
    im(:,:,zr) = readtiff([rt fn{i}]);
    i
end

% first filter: to remove the backbones of the soma, dendrites
T_den = 8;
tmpim = im;
tmpim(im>T_den) = 0; % set part of the dendridic backbone to zero
CC = bwconncomp(tmpim, 26);
MaxVoxelVolume = 150; % this # accounts for ~3-5um lenght
MinVoxelVolume = 10; % accounts of debris

csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=MinVoxelVolume & csize<MaxVoxelVolume; % keeping the short objects with large distance transforms
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
imGL = labelmatrix(CC)~=0;
imGL(im<=T_den & im>0) = 1; % recovering all short distance transforms such that we only have short distance transform data and 
tmpim = im;
tmpim(~imGL) = 0;
writetiff(uint8(tmpim(:,:,1:100)), [rt 'test.tif'], 'Compression', 'none');

%% temp workspace
enodes_idx = [node.ep];
enodes = node(logical(enodes_idx));

s = [4113, 3170, 100];
imnode = zeros(s, 'logical');



np = arrayfun(@(x) numel(x.idx), enodes);
idx = arrayfun(@(x) x.idx, enodes, 'unif',0);
idx = vertcat(idx{:});

linkidx = [enodes.links];
nllen = arrayfun(@(x) numel(x.point), link);
nllen = nllen(linkidx);
linkpts = arrayfun(@(x) x.point, link, 'unif',0);
linkpts = linkpts(linkidx);
linkpts = horzcat(linkpts{:});


imnode(linkpts) = 1;
writetiff(uint8(imnode(:,:,1:100)), [rt 'test_nodes.tif'], 'Compression', 'none');


%% yfp associated homer - create dialated masked regions of spines

% read the distance transform valued skeleton
rt = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/ch0_Analysis_1ch/';
rtYFP = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/Position2_tile/Stitch_Igor/Converted/slice-tiff_fusion/ch0/';

s = [4113, 3170, 700];
im = zeros(s, 'uint16');

Homer = readtiff([rt 'YFPass_HomerJuxB.tif']);
Homer = logical(Homer);
Homer = Homer(:,:,1:750);

parfor_progress(750);
parfor i = 1:750
   im(:,:,i) = readtiff([rtYFP num2str(i-1) '.tif']);
    parfor_progress;
end
parfor_progress(0);

% dilate Homer postion
SE = strel('sphere',15);
dH = imdilate(Homer, SE);

imdH = zeros(size(dH), 'uint16');
imdH(dH) = im(dH);

writetiff(imdH, [rt 'homerDilated15YFPSignal.tif']);
writetiff(uint8(dH), [rt 'dilated15_YFPassHomer.tif']);

