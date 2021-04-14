function GU_ExM_calcPunctaDensityVol_1chVol(vol,im, varargin)
% calculated the synapse density in a given volume
% input: vol = punta/discrete structures; input for this function is a
% volume
% updated compared to the original -- integrated code from GU_ExM_correctDANSynapseDensityCalc
% also integrated GU_ExM_JaneliaCluster_ExtractDANassociatednc82
% Gokul Upadhyayula, Nov, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isstr);
ip.addRequired('im');
ip.addParameter('ExpansionFactor', 4.09, @isnumeric);
ip.addParameter('PixelSize', 0.097, @isnumeric); % in um
ip.addParameter('PixelSizeZ', 0.18, @isnumeric); % in um
ip.addParameter('Sigmas', [2.5, 2.5], @isnumeric); % necessary for 3D Gauss Fitting
ip.addParameter('DensityVolume', 1, @isnumeric); % pre-expanded surveyed volume; in um^3
ip.addParameter('minAZdia', 0.1, @isnumeric); % in nm, pre-expanded - to remove active zone duplicates; only when FindLocalMaxima is true
ip.addParameter('FindLocalMaxima', true, @islogical); % minAZdia not used to remove duplicate if set to false
% image filtering options
ip.addParameter('Threshold', [], @isvector); % if left empty, will use OTSU to calculate intensity to threshold; [ch1, ch2]
ip.addParameter('MinThreshold', [0,0], @isvector); % if left empty, will use OTSU to calculate intensity to threshold; [ch1, ch2]
ip.addParameter('OTSUGaussKernel', 1, @isnumeric); % 3D Gaussian kernel size to smooth image
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('MinVoxelVolume', [516, 1008], @isvector); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise; second value is for the DAN channel
ip.addParameter('StrelSphereSize', 2, @isnumeric); % to make spheres from the local max
ip.addParameter('Verbose', true, @islogical); % to make spheres from the local max
ip.addParameter('MinVolumeOccupancy', 0.1, @isnumeric); % the loaded volume must have fraction of voxels with values > 0

ip.parse(vol,im, varargin{:});
pr = ip.Results;

%% options
minAZdia = pr.minAZdia;
sigmas = pr.Sigmas;
ef = pr.ExpansionFactor;
px = pr.PixelSize; %in um
bandwidth = minAZdia/2/px*ef; % search radius
V = pr.DensityVolume;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % radius in expanded samples; in um
bandwidth2 = re/px;
T = pr.Threshold;
MinVoxelVolume = pr.MinVoxelVolume;
se = strel('sphere', pr.StrelSphereSize);
Vrb = pr.Verbose;
MVO = pr.MinVolumeOccupancy;
zAniso = pr.PixelSizeZ/px;
tic
%% read volumes

% if Vrb
%     fprintf('reading volume...')
%     
% end
% 
% fileInfo = dir(vol);
% fileSize = fileInfo.bytes/1000^3;
% if fileSize < 4
%     im = readtiff(vol);
% else
%     [im,~] = imread_big(vol);
% end

ProceedForward = logical(sum(im(:)) > 0) & (sum(im(:)>0) >= numel(im)*MVO);
if ProceedForward
    if Vrb
        toc
    end
    
    %% filtering volumes
    
    if Vrb
        fprintf('filtering volumes...')
        tic
    end
    
    if isempty(T)
        tic
        imG = filterGauss3D(im,pr.OTSUGaussKernel); toc
        PP = prctile(im(:),pr.OTSUMaxPer);
        if PP > 0
            T(1) = thresholdOtsu(im(im > 0 & im < PP));
        else
            T(1) = 0;
            ProceedForward = 0;
        end
        OT(1) = T(1);
        if T(1) < pr.MinThreshold(1)
            T(1) = pr.MinThreshold(1);
        end
        imG = imG-T(1);
        imG(imG<0) = 0;
    else
        imG = im-T(1);
    end
    
    if Vrb
        toc
    end
end

if ProceedForward
    %% clean volumes and discard random voxels
    
    if MinVoxelVolume(1) > 0
        if Vrb
            fprintf('discarding small noise voxels for VOL1 channel...')
            tic
        end
        
        CC = bwconncomp(logical(imG), 26);
        csize = cellfun(@numel, CC.PixelIdxList);
        idx = csize>=MinVoxelVolume(1);
        CC.NumObjects = sum(idx);
        CC.PixelIdxList = CC.PixelIdxList(idx);
        imGL = labelmatrix(CC)~=0;
        imG = uint16(imG) .* uint16(imGL);
        
        if Vrb
            toc
        end
    end
    if sum(imG(:) > 0) ==0
        ProceedForward = 0;
    end
end

if ProceedForward
    %% calculate local maxima and centroids of the synapses
    
    if Vrb
        fprintf('calculating local maxima volume...')
        tic
    end
    if pr.FindLocalMaxima
        allmax = locmax3d(imG, 2*ceil(sigmas([1 1 2]))+1, 'ClearBorder', false);
        [y,x,z] = ind2sub(size(allmax), find(allmax));
        ptData(:,1) = x';
        ptData(:,2) = y';
        ptData(:,3) = z';
    else
        s  = regionprops(logical(imG),'centroid');
        ptData = cat(1, s.Centroid);
    end
    
    if Vrb
        toc
    end
    
    %% filter detections closer than pre-set values (Default: 100nm)
    % take the mean of two points closer than the pre-set value
    if pr.FindLocalMaxima
        if Vrb
            fprintf('kd tree search to remove duplicates (diameter pre expanded ~ see options)*expansion factor...')
            tic
        end
        
        Mdl = KDTreeSearcher(ptData);
        [~,D] = knnsearch(Mdl, ptData, 'k',2);
        if isempty(D) || size(D,2) == 0 || size(ptData,1) ==1
            ProceedForward = 0;
        end
    end
end


if ProceedForward
    if pr.FindLocalMaxima % only remove duplicates if using local maxima
        ptDistances = D(:,2);
        rmidx = ptDistances <= bandwidth*2; % index of duplicates
        
        [idx,~] = knnsearch(Mdl, ptData(rmidx,:), 'k',2); % find duplicate pairs
        for nn = 1:length(idx)
            idx(idx(:,1) == idx(nn,2),1) = idx(nn,1);
        end
        
        % calculate the average of duplicates and remove one from the pairs
        [~,iaF,~] = unique(idx(:,1), 'first');
        [~,iaL,~] = unique(idx(:,1), 'last');
        
        dupPts(:,:,1) = ptData(idx(iaF,2),:);
        dupPts(:,:,2) = ptData(idx(iaL,1),:);
        meanptData = mean(dupPts,3);
        filteredData = [ptData(~rmidx,:); round(meanptData)];
        
        if Vrb
            toc
        end
    else
        filteredData = ptData;
    end
    %% density calculation with kdtree range search
    
    if Vrb
        fprintf('range search of every synapse and calculating density...')
        tic
    end
    
    filteredData_Zcorr(:,1) = filteredData(:,1);
    filteredData_Zcorr(:,2) = filteredData(:,2);
    filteredData_Zcorr(:,3) = filteredData(:,3) * zAniso;
    
    MdlKDT = KDTreeSearcher(filteredData_Zcorr);
    IdxKDT = rangesearch(MdlKDT,filteredData_Zcorr,bandwidth2);% Search radius = bandwidth2;
    filteredsynapseDensity = cellfun(@numel, IdxKDT)/V;
    
    if Vrb
        toc
    end
    
    %% Generate intensity coded synpase density
    
    if Vrb
        fprintf('generating intenisty coded raw mask data ...')
        tic
    end
end
mask = zeros(size(im), 'uint8');
fprintf('initializing new volumes...'); tic
toc

if ProceedForward
    for ii = 1:numel(filteredData(:,2))
        mask(round(filteredData(ii,2)),round(filteredData(ii,1)),round(filteredData(ii,3))) = filteredsynapseDensity(ii); % #synapses/volume
    end
    mask = imdilate(mask, se);
    
    if Vrb
        toc
    end
    
    %% write labeled image mask
    
    if Vrb
        fprintf('clearing workspace and writing intensity coded image and saving structure info...')
        tic
    end
    
end

cleanVol1 = zeros(size(im), 'uint16');

if ProceedForward
    cleanVol1(logical(imG)) = im(logical(imG));
end

[p1,fn1,~] = fileparts(vol);

if ~exist([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep ], 'dir')
    mkdir([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep ]);
end

if ~exist([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density'], 'dir')
    mkdir([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density']);
end

writetiff(cleanVol1, [p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean.tif']);

clear im im2 imG imG2 imGL allmax CC D dupPts iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn
clear rmidx varargin x y z idx

writetiff(mask, [p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_density.tif']);


clear mask DANmask

if Vrb
    toc
end

%% integrating the GU_ExM_correctDANSynapseDensityCalc script below
% GU_ExM_correctDANSynapseDensityCalc(matfile, cleanDANpath, varargin)

if Vrb
    fprintf('clearing workspace and writing intensity coded image and saving structure info...')
    tic
end

if ~exist([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures'], 'dir')
    mkdir([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures']);
end

if ~exist([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density'], 'dir')
    mkdir([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density']);
end


clear DANmask OcleanDAN ocleanDAN cleanDAN rmidx se varargin x y z idx cleanVol1 cleanVol2 CC im im2 imG imG2 imGL allmax CC D dupPts iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn
clear cleanVol1 imDAN imNONDAN filteredData_DAN ip

save([p1 filesep 'Analysis_1ch' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_densityData.mat']);

if Vrb
    toc
end
end

