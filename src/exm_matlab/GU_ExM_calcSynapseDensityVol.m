function GU_ExM_calcSynapseDensityVol(vol, vol2, varargin)
% calculated the synapse density in a given volume
% input: volume

% Gokul Upadhyayula, Oct, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isstr);
ip.addRequired('vol2', @isstr); % second channel
ip.addParameter('ExpansionFactor', 4.09, @isnumeric);
ip.addParameter('PixelSize', 0.097, @isnumeric); % in um
ip.addParameter('Sigmas', [2.5, 2.5], @isnumeric); % necessary for 3D Gauss Fitting
ip.addParameter('DensityVolume', 1, @isnumeric); % pre-expanded surveyed volume; in um^3
ip.addParameter('minAZdia', 0.1, @isnumeric); % in um, pre-expanded - to remove active zone duplicates
% image filtering options
ip.addParameter('Threshold', [], @isvector); % if left empty, will use OTSU to calculate intensity to threshold; [ch1, ch2]
ip.addParameter('MinThreshold', [0,0], @isvector); % if left empty, will use OTSU to calculate intensity to threshold; [ch1, ch2]
ip.addParameter('OTSUGaussKernel', 1, @isnumeric); % 3D Gaussian kernel size to smooth image
ip.addParameter('OTSUMaxPer', 99, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('MinVoxelVolume', 192, @isnumeric); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise (168 comes from sig = 2.3)
ip.addParameter('StrelSphereSize', 2, @isnumeric); % to make spheres from the local max
ip.addParameter('Verbose', false, @islogical); % to make spheres from the local max
ip.addParameter('MinVolumeOccupancy', 0.1, @isnumeric); % the loaded volume must have fraction of voxels with values > 0
ip.parse(vol, vol2, varargin{:});
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
%% read volumes

if Vrb
    fprintf('reading volumes...')
    tic
end

im = readtiff(vol);
im2 = readtiff(vol2);
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
        imG = filterGauss3D(im,pr.OTSUGaussKernel);
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
        
        % for vol 2 (DAN)
        imG2 = filterGauss3D(im2,pr.OTSUGaussKernel);
        PP = prctile(im2(:),pr.OTSUMaxPer);
        if PP > 0
            T(2) = thresholdOtsu(im2(im2 > 0 & im2 < PP));
        else
            T(2) = 0;
        end
        OT(2) = T(2);
        if T(2) < pr.MinThreshold(2)
            T(2) = pr.MinThreshold(2);
        end
        imG2 = imG2-T(2);
        imG2(imG2<0) = 0;
    else
        imG = im-T(1);
        imG2 = im2-T(2);
    end
    
    if Vrb
        toc
    end
end

if ProceedForward
    %% clean volumes and discard random voxels
    
    if MinVoxelVolume > 0
        if Vrb
            fprintf('discarding small noise voxels for synapse channel...')
            tic
        end
        
        CC = bwconncomp(logical(imG), 26);
        csize = cellfun(@numel, CC.PixelIdxList);
        idx = csize>=MinVoxelVolume;
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
    allmax = locmax3d(imG, 2*ceil(sigmas([1 1 2]))+1, 'ClearBorder', false);
    [y,x,z] = ind2sub(size(allmax), find(allmax));
    ptData(:,1) = x';
    ptData(:,2) = y';
    ptData(:,3) = z';
    
    if Vrb
        toc
    end
    
    %% filter detections closer than pre-set values (Default: 100nm)
    % take the mean of two points closer than the pre-set value
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


if ProceedForward
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
    
    clear Mdl D
    Mdl = KDTreeSearcher(filteredData);
    [~,D] = knnsearch(Mdl, filteredData, 'k',2);
    filteredDistances = D(:,2);
    
    if Vrb
        toc
    end
    
    %% density calculation with kdtree range search
    
    if Vrb
        fprintf('range search of every synapse and calculating density...')
        tic
    end
    
    MdlKDT = KDTreeSearcher(filteredData);
    IdxKDT = rangesearch(MdlKDT,filteredData,bandwidth2);% Search radius = bandwidth2;
    filteredsynapseDensity = cellfun(@numel, IdxKDT)/V;
    
    if Vrb
        toc
    end
    % %% check figures
    % figure,
    % scatter3(filteredData(:,1),filteredData(:,2),filteredData(:,3),20,synapseDensity, 'filled')
    % colormap jet
    
    %% Generate intensity coded synpase density
    
    if Vrb
        fprintf('generating intenisty coded raw mask data & read out DAN signal...')
        tic
    end
end
mask = zeros(size(im), 'uint8');
DANmask = zeros(size(im), 'uint8');

if ProceedForward
    DANSynapseIdx = zeros(numel(filteredData(:,2)),1, 'logical');
    DANchSignal = zeros(numel(filteredData(:,2)),1, 'uint16');
    for ii = 1:numel(filteredData(:,2))
        mask(filteredData(ii,2),filteredData(ii,1),filteredData(ii,3)) = filteredsynapseDensity(ii); % #synapses/volume
        DANchSignal(ii) = imG2(filteredData(ii,2),filteredData(ii,1),filteredData(ii,3));
        if DANchSignal(ii) > 0
            DANmask(filteredData(ii,2),filteredData(ii,1),filteredData(ii,3)) = filteredsynapseDensity(ii); % #synapses/volume
            DANSynapseIdx(ii) = 1;
        end
    end
    mask = imdilate(mask, se);
    DANmask = imdilate(DANmask, se);
    
    if Vrb
        toc
    end
    
    %% calculate DAN synapse distance to non DAN
    if Vrb
        tic
    end
    
    Mdl = KDTreeSearcher(filteredData(~DANSynapseIdx,:)); % all nonDAN synapse positions
    [~,D] = knnsearch(Mdl, filteredData(DANSynapseIdx,:), 'k',1);
    if ~isempty(D) && size(D,2) > 1
        Distance_DANSyn2NONDANSyn = D(:,1);
    end
    
    Mdl = KDTreeSearcher(filteredData(DANSynapseIdx,:)); % all DAN synapse positions
    [~,D] = knnsearch(Mdl, filteredData(DANSynapseIdx,:), 'k',2);
    if ~isempty(D) && size(D,2) > 1
        Distance_DANSyn2DANSyn = D(:,2);
    end
    
    Mdl = KDTreeSearcher(filteredData(~DANSynapseIdx,:)); % all DAN synapse positions
    [~,D] = knnsearch(Mdl, filteredData(~DANSynapseIdx,:), 'k',2);
    if ~isempty(D) && size(D,2) > 1
        Distance_NONDANSyn2NONDANSyn = D(:,2);
    end
    
    %% misc calculations
    
    nSynCloserthanSetValue = numel(ptData(:,2))-numel(filteredData(:,2)); % number of synapses closer than pre-set value
    
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
cleanVol2 = zeros(size(im), 'uint16');

if ProceedForward
    cleanVol1(logical(imG)) = im(logical(imG));
    cleanVol2(logical(imG2)) = im2(logical(imG2));
end
[p1,fn1,~] = fileparts(vol);
[p2,fn2,~] = fileparts(vol2);

if ~exist([p1 filesep 'Analysis'], 'dir')
    mkdir([p1 filesep 'Analysis']);
end

if ~exist([p2 filesep 'Analysis'], 'dir')
    mkdir([p2 filesep 'Analysis']);
end
if ~exist([p1 filesep 'Analysis' filesep 'DataStructures'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep 'DataStructures']);
end
if ~exist([p1 filesep 'Analysis' filesep 'Density'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep 'Density']);
end
if ~exist([p1 filesep 'Analysis' filesep 'DANDensity'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep 'DANDensity']);
end

writetiff(cleanVol1, [p1 filesep 'Analysis' filesep fn1 '_clean.tif']);
writetiff(cleanVol2, [p2 filesep 'Analysis' filesep fn2 '_clean.tif']);

clear im im2 imG imG2 imGL allmax CC D dupPts iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn
clear rmidx se varargin x y z idx pr cleanVol1 cleanVol2

writetiff(mask, [p1 filesep 'Analysis' filesep 'Density' filesep fn1 '_density.tif']);
writetiff(DANmask, [p1 filesep 'Analysis' filesep 'DANDensity' filesep fn1 '_DANdensity.tif']);

clear mask DANmask
save([p1 filesep 'Analysis' filesep 'DataStructures' filesep fn1 '_densityData.mat']);

if Vrb
    toc
end

