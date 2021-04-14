%% DPM3 

MinHoleVolume = 1008;% for membrane
MinHoleVolume = 246; % for puncta

for i = 1%1:numel(data)
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
        
        
        fprintf('writing image...')
        tic
        writetiff(imG, [data(i).framePathsDS{j}{1}]);
        
        toc
        
        clear im imG imGL;
    end
end
% check signal intensity of the clean mask in 1px eroded region

se = strel('sphere',1);
cleanDAN_E = imerode(logical(imG), se);
figure, histogram(imG(imG>0 & ~cleanDAN_E ),1000);
muCut = prctile(imG(imG>0 & ~cleanDAN_E),99.9);

% ch0 - PN threshold = 4542 -->otsu (2980, 2021) +99.9 (2134, 1949) prctle of the edge voxels (from 1);
% ch1 - nc82 - 960: mean([1000,917])


%% mouse visual cortex workspace

% rt = '/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/';
s = GU_calcGaussianIntegral(1, [2.5 2.5 2.5]);
parfor k = 1:numel(data)
GU_ExM_calcPunctaDensityVol(data(k).framePaths{2}{1}, data(k).framePaths{1}{1},...
    'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [round(s), 1008],'ExpansionFactor', 3.99,...
    'MinThreshold', [960,2500]); % sigma 3.9
end
%%

mx = 10617;
my = 10712;
mz = 6638;
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

%% calculate the local maxia of the puncta and clean the data
for i = 1:nn-1
    submit = ['bsub -n 3 -R"affinity[core(1)]" -J' ' "DPM3' num2str(i) '" -o ~/logs/DPM3' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_FlyDPM3_Puncta_cleanup ' num2str(i) ''''];
    [~, ~] = system(submit,'-echo');
end

%% submit 
clear rt
% rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0/Analysis/1008/';
rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2nonAss/';
rt{2} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2Ass/';
a = 0;
% rtsuffix = {'_clean', '_246_clean_V2nonAss', '_246_clean_V2Ass'};
rtsuffix = {'_246_clean_V2nonAss', '_246_clean_V2Ass'};
for k = 1:numel(rt)
    for i = 1:nn-1
        a = a +1;
        
        fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) rtsuffix{k} '.tif'];
        
        submit = ['bsub -n 10 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "RotDPM' num2str(a) '" -o ~/logs/RotDPM' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_JaneliaCluster_rotateTiles ' [rt{k} fn]  ''''];
        [stat, res] = system(submit,'-echo');
    end
end

%% change coordinates
mx = 10617;
my = 10712;
mz = 6638;
zxRatio = 0.18/0.097;
theta = 32;
clear rt
% rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0/Analysis/1008/RotatedStacks/';
rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2nonAss/RotatedStacks/';
rt{2} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2Ass/RotatedStacks/';
fileSufix = {'_246_clean_V2nonAss', '_246_clean_V2Ass'};
for k = 1:numel(rt)
    GU_ExM_WriteRotatedCoordinates(rt{k}, mx, my, mz, zxRatio, theta, fileSufix{k})
end

%% spectrally mix

   for i = 2000:7500
        submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "SMDPM' num2str(i) '" -o ~/logs/RotDPM' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_SpectralMixFlyDPM ' num2str(i)  ''''];
        [stat, res] = system(submit,'-echo');
   end
    
   %%
   
   rt = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0/Analysis/1008/RotatedStacks/allRenamedTiles/slice-tiff/ch0/';
 rtsave = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/SpectrallyMixed/';
tic
for i = 2000:7500
   movefile( [rt num2str(i) '.tif'],[rtsave 'z' num2str(i-2000,'%04d') '_c0_t0.tif']) 
end
toc