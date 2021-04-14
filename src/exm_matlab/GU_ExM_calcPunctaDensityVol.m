function GU_ExM_calcPunctaDensityVol(vol, vol2, varargin)
% calculated the synapse density in a given volume
% input: vol = punta/discrete structures; vol2 = membrane/cytosolic filled
% updated compared to the original -- integrated code from GU_ExM_correctDANSynapseDensityCalc
% also integrated GU_ExM_JaneliaCluster_ExtractDANassociatednc82
% Gokul Upadhyayula, Nov, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isstr);
ip.addRequired('vol2', @isstr); % second channel
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
ip.addParameter('AssSphereSize', 10, @isnumeric); % to make spheres from the local max
ip.addParameter('Verbose', true, @islogical); % to make spheres from the local max
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
zAniso = px/pr.PixelSizeZ;
%% read volumes

if Vrb
    fprintf('reading volumes...')
    tic
end

fileInfo = dir(vol);
fileSize = fileInfo.bytes/1000^3;s
if fileSize < 4
    im = readtiff(vol);
else
    [im,~] = imread_big(vol);
end

fileInfo = dir(vol2);
fileSize = fileInfo.bytes/1000^3;
if fileSize < 4
    im2 = readtiff(vol2);
else
    [im2,~] = imread_big(vol2);
end

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
    allptData = ptData;
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
        fprintf('generating intenisty coded raw mask data & read out VOL2 signal...')
        tic
    end
end
mask = zeros(size(im), 'uint8');
% DANmask = zeros(size(im), 'uint8');
fprintf('initializing new volumes...'); tic
imDAN = zeros(size(im), 'uint16');
imNONDAN = zeros(size(im), 'uint16');
toc

if ProceedForward
    DANchSignal = zeros(numel(filteredData(:,2)),1, 'uint16');
    for ii = 1:numel(filteredData(:,2))
        mask(round(filteredData(ii,2)),round(filteredData(ii,1)),round(filteredData(ii,3))) = filteredsynapseDensity(ii); % #synapses/volume
        DANchSignal(ii) = imG2(round(filteredData(ii,2)),round(filteredData(ii,1)),round(filteredData(ii,3)));
        
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
cleanVol2 = zeros(size(im), 'uint16');

if ProceedForward
    cleanVol1(logical(imG)) = im(logical(imG));
    cleanVol2(logical(imG2)) = im2(logical(imG2));
end


[p1,fn1,~] = fileparts(vol);
[p2,fn2,~] = fileparts(vol2);

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep ], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep ]);
end

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2Ass' filesep ], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2Ass' filesep ]);
end

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2nonAss' filesep ], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2nonAss' filesep ]);
end

% if ~exist([p2 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep ], 'dir')
%     mkdir([p2 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep ]);
% end

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density']);
end

writetiff(cleanVol1, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean.tif']);

clear im im2 imG imG2 imGL allmax CC D dupPts iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn
clear rmidx varargin x y z idx

writetiff(mask, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Density' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_density.tif']);


clear mask DANmask

if Vrb
    toc
end

%% integrating the GU_ExM_correctDANSynapseDensityCalc script below
% GU_ExM_correctDANSynapseDensityCalc(matfile, cleanDANpath, varargin)

% directories
cleanDANpath = [p2 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(2)) filesep fn2 '_clean.tif'];

if ~exist([p2 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(2))], 'dir')
    mkdir([p2 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(2))])
end


if Vrb
    fprintf('clearing workspace and writing intensity coded image and saving structure info...')
    tic
end
% [p2,fn2,~] = fileparts(cleanDANpath);

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures']);
end

if ~exist([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density'], 'dir')
    mkdir([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density']);
end

cleanDAN = cleanVol2;
clear cleanVol2
OcleanDAN = cleanDAN;
%% remove small debris from cleanDAN file
se2 = strel('sphere', pr.AssSphereSize);
if sum(cleanDAN(:)>0) > 0
    if Vrb
        fprintf('removing debris for VOL2 Data...')
        tic
    end
    CC = bwconncomp(cleanDAN, 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    
    idx = csize>=pr.MinVoxelVolume(2);
    CC.NumObjects = sum(idx);
    CC.PixelIdxList = CC.PixelIdxList(idx);
    cleanDAN = labelmatrix(CC)~=0;
    OcleanDAN(~cleanDAN) = 0;
    cleanDAN = imdilate(cleanDAN, se);
    cleanDAN = imdilate(cleanDAN, se);
    for i = 1:numel(filteredData(:,1))
        idx_smalldebris(i) = cleanDAN(round(filteredData(i,2)),round(filteredData(i,1)),round(filteredData(i,3)));
    end
    if Vrb
        toc
    end
    
    if Vrb
        fprintf('generating dilated centroid tiff stack...')
        tic
    end
    
    Vol2PunctaIdx = zeros(1,numel(filteredData(:,2)))';
    Vol2mask = zeros(size(cleanDAN), 'uint8');
    for ii = 1:numel(filteredData(:,2))
        if DANchSignal(ii) > 0 && idx_smalldebris(ii) == 1
            Vol2mask(round(filteredData(ii,2)),round(filteredData(ii,1)),round(filteredData(ii,3))) = filteredsynapseDensity(ii); % #synapses/volume
            Vol2PunctaIdx(ii) = 1;
        end
    end
    Vol2PunctaIdx = logical(Vol2PunctaIdx);
    Vol2mask = imdilate(Vol2mask, se);
    Vol2Puncta2Total = sum(Vol2PunctaIdx)/numel(filteredData(:,1));
    if Vrb
        toc
    end
    
    % generate VOL2 associated and non-associated VOL1 object file
    fprintf('generating Vol2 associated Vo1 objects file...'); tic
    
    filteredData_DAN = round(filteredData(Vol2PunctaIdx,:));
    toc
    
%     fprintf('geneating labels...'); tic
%     lholes = bwlabeln(logical(cleanVol1));
%     toc
    
    fprintf('looping through coordinates...'); tic
    DM = zeros(size(cleanVol1), 'logical');
    DM(sub2ind(size(cleanVol1),filteredData_DAN(:,2),filteredData_DAN(:,1),filteredData_DAN(:,3))) = 1;
    DM = imdilate(DM,se2);
    imDAN(DM(:,:,:)) = cleanVol1(DM(:,:,:));
    imNONDAN(~DM(:,:,:)) = cleanVol1(~DM(:,:,:));

%     danPuncatIdx = zeros(numel(filteredData_DAN(:,1)),1);
%     
%     for n = 1:numel(filteredData_DAN(:,1))
%         danPuncatIdx(n) = lholes(filteredData_DAN(n,2),filteredData_DAN(n,1),filteredData_DAN(n,3));
%     end
%     idx = ismember(lholes,danPuncatIdx);
%     imDAN(idx) = cleanVol1(idx);
%     imNONDAN(~idx) = cleanVol1(~idx);
    toc
    
    %% save data
    if ~exist(cleanDANpath, 'file')
        writetiff(OcleanDAN, cleanDANpath);
    end
    writetiff(Vol2mask, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density' filesep fn2 '_' num2str(pr.MinVoxelVolume(1)) '_Density.tif']);
    writetiff(imDAN, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2Ass' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean_V2Ass.tif']);
    writetiff(imNONDAN, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2nonAss' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean_V2nonAss.tif']);
    
else
    if Vrb
        fprintf('nothing there...saving blank volumes...')
        tic
    end
    
    if ~exist(cleanDANpath, 'file')
        writetiff(OcleanDAN, cleanDANpath);
    end
    writetiff(zeros(size(OcleanDAN),'uint8'), [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'Vol2Density' filesep fn2 '_' num2str(pr.MinVoxelVolume(1)) '_Density.tif']);
    writetiff(imDAN, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2Ass' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean_V2Ass.tif']);
    writetiff(imNONDAN, [p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) '_V2nonAss' filesep fn1 '_' num2str(pr.MinVoxelVolume(1)) '_clean_V2nonAss.tif']);
    
    if Vrb
        toc
    end
end

clear DANmask OcleanDAN ocleanDAN cleanDAN rmidx se varargin x y z idx cleanVol1 cleanVol2 CC im im2 imG imG2 imGL allmax CC D iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn
clear cleanVol1 imDAN imNONDAN

save([p1 filesep 'Analysis' filesep num2str(pr.MinVoxelVolume(1)) filesep 'DataStructures' filesep fn2 '_' num2str(pr.MinVoxelVolume(1)) '_densityData.mat']);

if Vrb
    toc
end
end

