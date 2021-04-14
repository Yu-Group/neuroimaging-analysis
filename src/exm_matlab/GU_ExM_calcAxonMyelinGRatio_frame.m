function GU_ExM_calcAxonMyelinGRatio_frame(data, varargin)

% Gokul Upadhyayula, April 2018

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('Sigma', [1,1], @isvector); % two values for 3D Gaussian smoothing
ip.addParameter('TValue', [0.1, 0.9], @isvector);
ip.addParameter('Mask', []);
ip.addParameter('MinBranchLength', 50, @isnumeric);
ip.addParameter('zAniso', 1.8557, @isnumeric);
ip.addParameter('OTSUMaxPer', 99.5, @isnumeric);
ip.addParameter('ExpansionFactor', 3.95, @isnumeric);
ip.addParameter('MyelinOffset', 0, @isnumeric); % # frames
ip.addParameter('LargestObject', false, @islogical);
ip.addParameter('UseCanny', false, @islogical);
ip.addParameter('Overwrite', true, @islogical);
ip.parse(data, varargin{:});
pr = ip.Results;

fnA = data.framePaths{1}{1};
fnM = data.framePaths{1}{2};
ef = pr.ExpansionFactor;
nF = pr.MyelinOffset;

rt = data.source;
cd(rt)
if isempty(pr.zAniso)
    pr.zAniso = data.dz/data.pixelSize;
end
%% read stacks
imA = readtiff(fnA);
imM = readtiff(fnM);

%% offset correction, if necessary
if nF > 0
    imMO = zeros([size(imM)]+[0 0 nF]);
    imAO = zeros([size(imA)]+[0 0 nF]);
    
    imMO(:,:,nF+1:end) = imM;
    imAO(:,:,1:end-nF) = imA;
    
    imA = imAO;
    imM = imMO;
end

im = {imA, imM};
%% 3d canny

for i = 1:2  % loop for axon and myelin
    tic
    fgsim{i} = filterGauss3D(im{i}, pr.Sigma);
    toc
    
    if pr.UseCanny
        tic
        [e{i}, thresh{i}] = canny(fgsim{i}, [0 0 0], ...
            'TMethod', 'relMax', 'TValue', pr.TValue, ...
            'SubPixel', false);
        toc
    else
        tim = fgsim{i};
        per99 = prctile(tim(:), pr.OTSUMaxPer);
        t = thresholdOtsu(tim(im{i}>0 & tim<per99));
        tim = tim-t;
        tim(tim<0) = 0;
        e{i} = logical(tim);
    end
    
    % mask data
    if ~isempty(pr.Mask)
        if nF > 0
            mask = zeros([size(pr.Mask)]+[0 0 nF]);
            if i == 1
                mask(:,:,1:end-nF) = pr.Mask;
            elseif i ==2
                mask(:,:,nF+1:end) = pr.Mask;
            end
            e{i}(~mask) = 0;
        else
            e{i}(~logical(pr.Mask)) = 0;
        end
    end
end

%% help patch mbp by adding axon
e{2} = logical(e{1}+e{2});

%% select largest segment
if pr.LargestObject
    for i = 1:2
        
        CC = bwconncomp(e{i}, 26);
                % discard small components (assumed to be noise or debris on glass slide)
                csize = cellfun(@numel, CC.PixelIdxList);
                CC.NumObjects = sum(csize == max(csize));
                CC.PixelIdxList = CC.PixelIdxList(csize == max(csize));
                fmask{i} = labelmatrix(CC)~=0;
    end
else
    fmask = e;
end


%% fill holes

HSIZE = 3;

se = strel('sphere',3);

for j = 1:numel(fmask)
    
    cm = fmask{j};
    
    for k = 1:3
        for z = 1:size(cm,3)
            cm(:,:,z) = imfill(cm(:,:,z), 8, 'holes');
        end
        for x = 1:size(cm,2)
            cm(:,x,:) = imfill(squeeze(cm(:,x,:)), 8, 'holes');
        end
        for y = 1:size(cm,1)
            cm(y,:,:) = imfill(squeeze(cm(y,:,:)), 8, 'holes');
        end
        %             cm = imclose(cm, binarySphere(HSIZE));
    end
    %         cm = imerode(cm,se);
    
    im{j} = cm;
end

imMyelin = im{2};
imAxon = im{1};

writetiff(uint8(imMyelin), [rt 'filled_3d_Myelin.tif']);
writetiff(uint8(imAxon), [rt 'filled_3d_Axon.tif']);

%% mex function
%pad edges
pimAxon = padarray(imAxon,[3 3 3],0,'both');
tic
skel = Skeleton3D(logical(pimAxon));
toc
skel_Axon = skel(4:end-3,4:end-3,4:end-3);
writetiff(uint8(skel), [rt 'Skel_filled_3d_Axon.tif']);

% pimMyelin = padarray(imMyelin,[3 3 3],0,'both');
% tic
% skel = Skeleton3D(logical(pimMyelin));
% toc
% % skel_Myelin = skel(4:end-3,4:end-3,4:end-3);
% writetiff(uint8(skel), [rt 'Skel_filled_3d_Myelin.tif']);
%%
% skel_Axon = readtiff([rt 'Skel_filled_3d_Axon.tif']);
% skel_Myelin = readtiff([rt 'Skel_filled_3d_Myelin.tif']);

skel = skel_Axon;
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
% for i = 1:3
[~,node2,link2] = Skel2Graph3D(skel,pr.MinBranchLength);
skel = Graph2Skel3D(node2,link2,w,l,h);
% end
skel_Axon_clean = skel;

% skel = skel_Myelin;
% syxz = size(skel);
% w = size(skel,1);
% l = size(skel,2);
% h = size(skel,3);
% for i = 1:3
%     [~,node2,link2] = Skel2Graph3D(skel,pr.MinBranchLength);
%     skel = Graph2Skel3D(node2,link2,w,l,h);
% end
%
% skel_Myelin_clean = skel;

writetiff(uint8(skel_Axon_clean), [rt 'Skel_filled_3d_Axon_clean.tif']);
% writetiff(uint8(skel_Myelin_clean), [rt 'Skel_filled_3d_Myelin_clean.tif']);
%%

tic

% interpolate
[ny,nx,nz] = size(skel_Axon_clean);
[y,x,z] = ndgrid(1:ny,1:nx,1:nz);
[Y,X,Z] = ndgrid(1:ny,1:nx,1:1/pr.zAniso:nz);
imask = interp3(x,y,z,double(skel_Axon_clean),X,Y,Z,'nearest');

iperim_imAxon = bwperim(interp3(x,y,z,double(imAxon),X,Y,Z,'nearest'));
iperim_imMyelin = bwperim(interp3(x,y,z,double(imMyelin),X,Y,Z,'nearest'));

dt_AxonCenter = bwdist(imask,'euclidean'); toc
writetiff(single(dt_AxonCenter), [rt 'Skel_dt_AxonCenter.tif']);
% dt_AxonCenter = readtiff([rt 'Skel_dt_AxonCenter.tif']);

% radius of Axon and Myelin Perim from Axon center
rad_imAxon = dt_AxonCenter .* single(iperim_imAxon);
rad_imMyelin = dt_AxonCenter .* single(iperim_imMyelin);
writetiff(uint8(rad_imAxon), [rt 'radius_imAxon.tif']);
writetiff(uint8(rad_imMyelin), [rt 'radius_imMyelin.tif']);

% rad_imAxon = readtiff([rt 'radius_imAxon.tif']);
% rad_imMyelin = readtiff( [rt 'radius_imMyelin.tif']);
% % distance from the Axon Center
% % tic
% % dt_AxonCenter = bwdist(skel_Axon_clean,'euclidean'); toc
% % writetiff(single(skel_Axon_clean), [rt 'DistanceTransform_from_Skel_filled_3d_Axon_clean.tif']);
% dt_AxonCenter = readtiff([rt 'DistanceTransform_from_Skel_filled_3d_Axon_clean.tif']);


%%
[radAxon_xyz(:,2),radAxon_xyz(:,1),radAxon_xyz(:,3)] = ind2sub(size(rad_imAxon),find(rad_imAxon(:)>0));
[radMyelin_xyz(:,2),radMyelin_xyz(:,1),radMyelin_xyz(:,3)] = ind2sub(size(rad_imMyelin),find(rad_imMyelin(:)>0));
tic
Mdl = KDTreeSearcher(radAxon_xyz);
toc
tic
[idx_Myelin,~] = knnsearch(Mdl, radMyelin_xyz, 'k',1);
toc

MyelinRadius = rad_imMyelin(rad_imMyelin(:)>0);
AxonRadius = rad_imAxon(rad_imAxon(:)>0);
gRatio = AxonRadius(idx_Myelin)./MyelinRadius;

gRatioVol = zeros(size(rad_imAxon), 'single');


corrGRatio = gRatio;
corrGRatio(corrGRatio>1) = 1;

for i = 1:numel(gRatio)
    gRatioVol(radMyelin_xyz(i,2),radMyelin_xyz(i,1),radMyelin_xyz(i,3)) = gRatio(i);
end

writetiff(single(gRatioVol), [rt 'GRATIO.tif']);
save([data.source 'corrGRatio.mat'], 'corrGRatio');

%% calculate myelin thickness

myelinThickness = MyelinRadius - AxonRadius(idx_Myelin);
myelinThickness(myelinThickness<0) = 0;
myelinThickness = myelinThickness*.097/ef;

% figure, scatter3(radMyelin_xyz(:,1),radMyelin_xyz(:,2), radMyelin_xyz(:,3),5,myelinThickness)
% colormap jet
% % caxis([0 1])
% xlabel('x')
% axis equal

%% write rotated data
%
rotatedCoord = radMyelin_xyz;
sx = ceil(max(rotatedCoord(:,1))-min(rotatedCoord(:,1))+1);
sy = ceil(max(rotatedCoord(:,2))-min(rotatedCoord(:,2))+1);
sz = ceil(max(rotatedCoord(:,3))-min(rotatedCoord(:,3))+1);

rotatedCoord(:,1) = rotatedCoord(:,1)-min(rotatedCoord(:,1))+1;
rotatedCoord(:,2) = rotatedCoord(:,2)-min(rotatedCoord(:,2))+1;
rotatedCoord(:,3) = rotatedCoord(:,3)-min(rotatedCoord(:,3))+1;

filename = [rt date 'GRATIO_Rotated.am'];
GU_amiraWriteSpots(filename,rotatedCoord, corrGRatio, [], [])


filename = [rt date 'myelinaxonDistance_um_Rotated.am'];
GU_amiraWriteSpots(filename,rotatedCoord, myelinThickness, [], [])
%% figure
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
ha = setupFigure(1,2, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [2 2.0 1.5], 'YSpace', [0.85 1 1]);
set(ha, 'FontSize', 6);

xr = [20];
AR = AxonRadius(idx_Myelin);
[~,unique_idx]= unique(idx_Myelin);
AR = AR(rotatedCoord(unique_idx,1)>xr & rotatedCoord(unique_idx,1)<(sx-xr) & rotatedCoord(unique_idx,2)>xr & rotatedCoord(unique_idx,2)<(sy-xr) & rotatedCoord(unique_idx,3)>xr & rotatedCoord(unique_idx,3)<(sz-xr))*0.097/ef*1000;
MR = MyelinRadius(rotatedCoord(unique_idx,1)>xr & rotatedCoord(unique_idx,1)<(sx-xr) & rotatedCoord(unique_idx,2)>xr & rotatedCoord(unique_idx,2)<(sy-xr) & rotatedCoord(unique_idx,3)>xr & rotatedCoord(unique_idx,3)<(sz-xr))*0.097/ef*1000;

axes(ha(1))

GU_stairsHistPlot({[],AR,MR}, 'ShowECDF', true, 'BinSize',25,'CreateFigure', false,...
    'PercentileEnd', 100, 'LineWidth',0.75,'ylabel', 'Relative Frequency','xlabel', 'Distance to Axon Center (um)')
set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')

legend('',['Axon n=' num2str(numel(AR))], '',['Myelin n=' num2str(numel(MR))])
xlim([0 800]);
xticks([0:250:1000])
yyaxis left
ylim([0 0.015])

axes(ha(2))
% GU_stairsHistPlot({[],[],[], corrGRatio}, 'ShowECDF', true, 'BinSize',0.05,'CreateFigure', false,...
%         'PercentileEnd', 100, 'LineWidth',0.75,'ylabel', 'Relative Frequency','xlabel', 'G-Ratio')
histogram(corrGRatio(rotatedCoord(unique_idx,1)>xr & rotatedCoord(unique_idx,1)<(sx-xr) & rotatedCoord(unique_idx,2)>xr & rotatedCoord(unique_idx,2)<(sy-xr) & rotatedCoord(unique_idx,3)>xr & rotatedCoord(unique_idx,3)<(sz-xr)), 25, 'EdgeColor', ceG, 'FaceColor', cfG, 'FaceAlpha',1)
xlabel('G-Ratio')
% yyaxis left
ylabel('Counts')
set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
xlim([0.2 1]);
xticks([0:0.2:1])

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' PLots_distancetoAxon_GRatio.eps']);