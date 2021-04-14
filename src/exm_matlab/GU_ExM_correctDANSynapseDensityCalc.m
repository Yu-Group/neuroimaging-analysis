function GU_ExM_correctDANSynapseDensityCalc(matfile, cleanDANpath, varargin)

% calculated the synapse density in a given volume
% input: volume
% Gokul Upadhyayula, Oct, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('matfile', @isstr);
ip.addRequired('cleanDANpath', @isstr); % second channel
ip.addParameter('ExpansionFactor', 4.09, @isnumeric);
ip.addParameter('PixelSize', 0.097, @isnumeric); % in um
ip.addParameter('DensityVolume', 1, @isnumeric); % pre-expanded surveyed volume; in um^3
ip.addParameter('StrelSphereSize', 2, @isnumeric); % to make spheres from the local max
% image filtering options
ip.addParameter('Threshold', 600, @isnumeric); % determined in GU_Expansion_nc83SynapseDensity
ip.addParameter('MinVoxelVolume', 246, @isnumeric); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise (168 comes from sig = 2.3)
ip.addParameter('Verbose', true, @islogical); % to make spheres from the local max
ip.parse(matfile, cleanDANpath, varargin{:});
pr = ip.Results;



%% load files
if pr.Verbose
    fprintf('reading volumes...')
    tic
end
load(matfile);
cleanDAN = readtiff(cleanDANpath);
if pr.Verbose
    toc
end
OcleanDAN = cleanDAN;

ef = pr.ExpansionFactor;
px = pr.PixelSize; %in um
V = pr.DensityVolume;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % radius in expanded samples; in um
bandwidth2 = re/px;
Vrb = pr.Verbose;
se = strel('sphere', pr.StrelSphereSize);

%% directories

if Vrb
    fprintf('clearing workspace and writing intensity coded image and saving structure info...')
    tic
end
[p2,fn2,~] = fileparts(cleanDANpath);


if ~exist([p2(1:end-6) filesep 'DataStructures_Recalc'], 'dir')
    mkdir([p2(1:end-6) filesep 'DataStructures_Recalc']);
end

if ~exist([p2(1:end-6) filesep 'DANDensityRecalc'], 'dir')
    mkdir([p2(1:end-6) filesep 'DANDensityRecalc']);
end


%% remove small debris from cleanDAN file

if sum(cleanDAN(:)>0) > 0
    if Vrb
        fprintf('removing debris for DAN Data...')
        tic
    end
    CC = bwconncomp(cleanDAN, 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    idx = csize>=pr.MinVoxelVolume;
    CC.NumObjects = sum(idx);
    CC.PixelIdxList = CC.PixelIdxList(idx);
    cleanDAN = labelmatrix(CC)~=0;
    OcleanDAN(~cleanDAN) = 0;
    
    for i = 1:numel(filteredData(:,1))
        idx_smalldebris(i) = cleanDAN(filteredData(i,2), filteredData(i,1), filteredData(i,3));
    end
    if Vrb
        toc
    end
    
    if Vrb
        fprintf('generating dilated DAN synapse tiff stack...')
        tic
    end
    DANSynapseIdx_original_preThreshold = DANSynapseIdx;
    DANSynapseIdx = zeros(1,numel(filteredData(:,2)))';
    DANmask = zeros(size(cleanDAN), 'uint8');
    for ii = 1:numel(filteredData(:,2))
        if DANchSignal(ii) > pr.Threshold && idx_smalldebris(ii) == 1
            DANmask(filteredData(ii,2),filteredData(ii,1),filteredData(ii,3)) = filteredsynapseDensity(ii); % #synapses/volume
            DANSynapseIdx(ii) = 1;
        end
    end
    DANSynapseIdx = logical(DANSynapseIdx);
    DANmask = imdilate(DANmask, se);
    DANSyn2Total = sum(DANSynapseIdx)/numel(filteredData(:,1));
    if Vrb
        toc
    end
    %%
    if Vrb
        fprintf('kd tree search for distances and density of synapses...')
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
    
    MdlKDT = KDTreeSearcher(filteredData(DANSynapseIdx,:));
    IdxKDT = rangesearch(MdlKDT,filteredData(DANSynapseIdx,:),bandwidth2);% Search radius = bandwidth2;
    DANsynapseDensity = cellfun(@numel, IdxKDT)/V;
    
    if Vrb
        toc
    end
    %% save data
    
    writetiff(OcleanDAN, [p2 filesep fn2 '.tif']);
    writetiff(DANmask, [p2(1:end-6) filesep 'DANDensityRecalc' filesep fn2 '_DANdensityRC.tif']);
else
    if Vrb
        fprintf('nothing there...saving blank volumes...')
        tic
    end
    writetiff(OcleanDAN, [p2 filesep fn2 '.tif']);
    writetiff(zeros(size(OcleanDAN),'uint8'), [p2(1:end-6) filesep 'DANDensityRecalc' filesep fn2 '_DANdensityRC.tif']);
    if Vrb
        toc
    end
end
clear DANmask ocleanDAN cleanDAN rmidx se varargin x y z idx pr cleanVol1 cleanVol2 CC im im2 imG imG2 imGL allmax CC D dupPts iaF iaL IdxKDT ii Mdl MdlKDT meanptData MinVoxelVolume nn

save([p2(1:end-6) filesep 'DataStructures_Recalc' filesep fn2 '_densityDataRC.mat']);

if Vrb
    toc
end
end