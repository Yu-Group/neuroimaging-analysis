function GU_ExpansionMapSynapticDensity(impath, varargin)
% calculate synaptic density map calcuation script
% Gokul Upadhyayula, 2017

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('impath', @isstr);

% general options
ip.addParameter('overwrite', false , @islogical); %
ip.addParameter('AnalysisPath', 'Analysis', @isstr); %
ip.addParameter('DensityMapPath', 'Density', @isstr); %
ip.addParameter('FiguresPath', 'Figures', @isstr); %
ip.addParameter('GenerateFigures', true , @islogical); %

% image preprocessing
ip.addParameter('PixelSizeXY', 0.097 , @isnumeric); % in um
ip.addParameter('PixelSizeZ', 0.18 , @isnumeric); % in um
ip.addParameter('sigma', [1.3, 1.4] , @isvector); %
ip.addParameter('MinHoleVolume', 50 , @islogical); %
ip.addParameter('GaussSigma', 1 , @islogical); %
ip.addParameter('SphereSigma', 2 , @islogical); %

% clustering points
ip.addParameter('ActiveZoneDiameter', 0.2 , @isnumeric); % in um and real space
ip.addParameter('ExpansionFactor', 4 , @isnumeric); % x fold expanded
ip.addParameter('SurveyedVolume', 1 , @isnumeric); % in um3 and real space
ip.addParameter('Threshold', 700 , @isnumeric);
ip.parse(impath, varargin{:});
pr = ip.Results;

% file name
[p, fn, ext] = fileparts(impath);


if ~exist([p filesep pr.DensityMapPath fn ext], 'file') || pr.overwrite
    % check folders exist
    mkdir([p filesep pr.AnalysisPath]);
    mkdir([p filesep pr.DensityMapPath]);
    mkdir([p filesep pr.FiguresPath]);
    
    % options
    sigma = pr.sigma;
    MinHoleVolume = pr.MinHoleVolume;
    AZdia = pr.ActiveZoneDiameter;% diameter of activezone
    ef = pr.ExpansionFactor;
    px = pr.PixelSizeXY; %in um
    pxZ = pr.PixelSizeZ;
    bandwidth = AZdia/2/px*ef;
    V = pr.SurveyedVolume;% in um3 surveyed volume
    r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
    re = r*ef; % expanded radius in um
    bandwidth2 = re/px;
    se = strel('sphere', pr.SphereSigma);
  
    % read image volume
    im = readtiff(impath);
    
    % filter image volume
    imG = filterGauss3D(im,pr.GaussSigma);
    if isempty(pr.Threshold)
        T = thresholdOtsu(im(im > 0 & im < prctile(im(:),99)));
    else
        T = pr.Threshold;
    end
    imG = imG-T;
    imG(imG<0) = 0;
    
    % discard small objects
    CC = bwconncomp(logical(imG), 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    idx = csize>=MinHoleVolume;
    CC.NumObjects = sum(idx);
    CC.PixelIdxList = CC.PixelIdxList(idx);
    imGL = labelmatrix(CC)~=0;
    imG = uint16(imG) .* uint16(imGL);
    clear imGL;
    
    % calculate local max
    im = logical(locmax3d(imG, 2*ceil(sigma([1 1 2]))+1, 'ClearBorder', false)); % im now just points for all local max
    clear imG
    
    % %%%%%%cluster analysis%%%%%%% %
    % get centroids
    s  = regionprops(im,'centroid');
    ptData = cat(1, s.Centroid);
    clear s imcentroids
    
    % first round of clustering to remove duplicate detections
    [clusterInfo, ~, ~] = MeanShiftClustering(ptData, bandwidth);
    ClusCenters = cat(1, clusterInfo.ptClusterCenter);
    
    % second round of clustering to group
    [clusterInfo2, ~, ~] = MeanShiftClustering(ClusCenters, bandwidth2);
    ClusCenters2 = cat(1, clusterInfo2.ptClusterCenter);
    nSynapses = numel(clusterInfo);
    
    % save matlab variables
    save([p filesep pr.AnalysisPath filesep fn '.mat'], 'clusterInfo', 'ClusCenters','ClusCenters2', 'clusterInfo2', 'nSynapses', 'ptData');
    
    % sort clusters based on size - used for color coding
    [~,idx] = sort([clusterInfo2.numPoints], 'descend');
    
    % generating intenisty coded raw mask data
    [lholes, ~] = bwlabeln(imdilate(uint8(im), se));
    lholes = (lholes);
    mask = zeros(size(im), 'single');
    parfor_progress(numel(idx));
    for n = 1:numel(idx)
        clear lidx
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
    
    % write Density map tiff 
    writetiff(mask, [p filesep pr.DensityMapPath filesep fn ext]);
    
    % make figures
    if pr.GenerateFigures
        cm = jet(numel([clusterInfo2.numPoints]));
        cm = flip(cm,2);
        ha = setupFigure(1,1, 'AxesWidth', 10, 'AxesHeight', 10,'SameAxes', false,...
            'XSpace', [1.5 2 1], 'YSpace', [1.25 1 1]);
        set(ha, 'FontSize', 8);
        for n = 1:numel(idx)
            x2= ClusCenters(clusterInfo2(idx(n)).ptIdData,1);
            y2= ClusCenters(clusterInfo2(idx(n)).ptIdData,2);
            z2= ClusCenters(clusterInfo2(idx(n)).ptIdData,3);
            scatter3(x2*px, y2*px, z2*px, 10, cm(n,:), 'filled','MarkerFaceAlpha',1); hold on
        end
        axis equal
        view(3)
        xlabel('x-axis (um^3)')
        ylabel('y-axis (um^3)')
        zlabel('z-axis (um^3)')
        toc
        f0 = gcf();
        print(f0, '-dsvg', [p filesep pr.FiguresPath filesep fn '_surveypreexpandedVol' num2str(V) 'um3.svg']);
    end
end