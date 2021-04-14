55%% March 27, 2018 ||||| reprocessing the data using the new stitched data
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/Raw/Axon1_Newstitch/';
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/Raw/Axon2_newstitch/';
cd(rt)

fn = dir('*.tif');
fn = {fn.name};
fnA = fn{1};
fnM = fn{2};

%%
imA = readtiff(fnA);
imM = readtiff(fnM);
imA = imA(:,:,450:850);
imM = imM(:,:,450:850);

%% 3d canny
tic
fgsimA = filterGauss3D(imA, [1,1]);
toc
tic
[e, thresh] = canny(fgsimA, [0 0 0], ...
          'TMethod', 'relMax', 'TValue', [0.1 0.9], ...
          'SubPixel', false);
      toc
%       imtool3D(uint16(e*10000)+imA)
spy3d(e)
imA_e = e;

%%
tic
fgsimM = filterGauss3D(imM, [1,1]);
toc
tic
[e, thresh] = canny(fgsimM, [0 0 0], ...
          'TMethod', 'relMax', 'TValue', [0.1 0.7], ...
          'SubPixel', false);
      toc
imM_e = e;
%% clean up, segment and discard debris
imAM = {imA, imM};
fn = {fnA, fnM};
MinHoleVolume = 1008;% for membrane
% MinHoleVolume = 516; % for puncta

for i = 1%1:numel(data)
    for j = 1:numel(imAM)
        im = imAM{j};
        fprintf('filtering image...')
        
        tic
        imG = filterGauss3D(im,1);
        T = thresholdOtsu(im(im > 0));
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
        im(imG==0) =0;
        
        fprintf('writing image...')
        tic
        writetiff(im, ['segmented' fn{j}]);
        toc
        clear im imG imGL;
    end
end
%%  read segmented masks
simA = readtiff(['segmented' fnA]);
simM = readtiff(['segmented' fnM]);

simA = simA(:,:,450:850);
simM = simM(:,:,450:850);
se = strel('sphere',3);
simA_e = imerode(simA, se);
bwl_A = bwlabeln(simA_e);

mask_A = zeros(size(bwl_A),'logical');
mask_A(bwl_A==1) = 1;
mask_A = imdilate(mask_A,se);

se = strel('sphere',15);
dmask_A = imdilate(mask_A, se);
mask_M = zeros(size(bwl_A),'logical');
mask_M(dmask_A) = logical(simM(dmask_A));

writetiff(uint8(mask_M), ['Masked_15px' fnM]);
writetiff(uint8(mask_A), ['Masked' fnA]);

% %% cleanCurated
% curatedM = readtiff([date 'curated_Masked_15pxAxon1_channel 1 [3519, 522, 5240]_newstitch.tif']);
% curatedM = bwlabeln(curatedM);
% cleanCuratedM = zeros(size(curatedM),'logical');
% cleanCuratedM(curatedM==1) = 1;
% writetiff(uint8(cleanCuratedM), [date 'curated_Masked_15pxAxon1_channel 1 [3519, 522, 5240]_newstitch.tif']);

%%

imM = readtiff(['Masked_15px' fnM]);
imA = readtiff(['Masked' fnA]);

% imM = imM(:,:,1:1104);
% imA = imA(:,:,1:1104);
HSIZE = 3;
im = {imM, imA};
se = strel('sphere',3);

for j = 1:numel(im)
    
    cm = im{j};
    cm = imdilate(cm,se);
    
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
        cm = imclose(cm, binarySphere(HSIZE));
    end
    cm = imerode(cm,se);
    
    im{j} = cm;
end
perim_imAxon = bwperim(im{2});
perim_imMyelin = bwperim(im{1});

% perim_imAxon = perim_imAxon(:,:,1:1085);
% perim_imMyelin = perim_imMyelin(:,:,1:1085);

imMyelin = im{1};
imAxon = im{2};

writetiff(uint8(imMyelin), [rt 'filled_3d_Myelin.tif']);
writetiff(uint8(imAxon), [rt 'filled_3d_Axon.tif']);

%% mex function
%pad edges
pimAxon = padarray(imAxon,[3 3 3],0,'both');
tic
skel = Skeleton3D(logical(pimAxon));
toc
skel = skel(4:end-3,4:end-3,4:end-3);
writetiff(uint8(skel), [rt 'Skel_filled_3d_Axon.tif']);

pimMyelin = padarray(imMyelin,[3 3 3],0,'both');
tic
skel = Skeleton3D(logical(pimMyelin));
toc
skel = skel(4:end-3,4:end-3,4:end-3);
writetiff(uint8(skel), [rt 'Skel_filled_3d_Myelin.tif']);
%%
skel_Axon = readtiff([rt 'Skel_filled_3d_Axon.tif']);
skel_Myelin = readtiff([rt 'Skel_filled_3d_Myelin.tif']);
% skip section
MinBranchLength = 50;

skel = skel_Axon;
syxz = size(skel);
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
for i = 1:3
    [~,node2,link2] = Skel2Graph3D(skel,MinBranchLength);
    skel = Graph2Skel3D(node2,link2,w,l,h);
end
skel_Axon_clean = skel;

skel = skel_Myelin;
syxz = size(skel);
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
for i = 1:3
    [~,node2,link2] = Skel2Graph3D(skel,MinBranchLength);
    skel = Graph2Skel3D(node2,link2,w,l,h);
end

skel_Myelin_clean = skel;

writetiff(uint8(skel_Axon_clean), [rt 'Skel_filled_3d_Axon_clean.tif']);
writetiff(uint8(skel_Myelin_clean), [rt 'Skel_filled_3d_Myelin_clean.tif']);
%%
zA = 0.18/0.097;
tic

% interpolate
[ny,nx,nz] = size(skel_Axon_clean);
[y,x,z] = ndgrid(1:ny,1:nx,1:nz);
[Y,X,Z] = ndgrid(1:ny,1:nx,1:1/zA:nz);
imask = interp3(x,y,z,double(skel_Axon_clean),X,Y,Z,'nearest');

iperim_imAxon = interp3(x,y,z,double(perim_imAxon),X,Y,Z,'nearest');
iperim_imMyelin = interp3(x,y,z,double(perim_imMyelin),X,Y,Z,'nearest');

dt_AxonCenter = bwdist(imask,'euclidean'); toc
writetiff(single(dt_AxonCenter), [rt 'Skel_dt_AxonCenter.tif']);
dt_AxonCenter = readtiff([rt 'Skel_dt_AxonCenter.tif']);

% radius of Axon and Myelin Perim from Axon center
rad_imAxon = dt_AxonCenter .* single(iperim_imAxon);
rad_imMyelin = dt_AxonCenter .* single(iperim_imMyelin);
writetiff(uint8(rad_imAxon), [rt 'radius_imAxon.tif']);
writetiff(uint8(rad_imMyelin), [rt 'radius_imMyelin.tif']);
%%
rad_imAxon = readtiff([rt 'radius_imAxon.tif']);
rad_imMyelin = readtiff( [rt 'radius_imMyelin.tif']);
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
%%
radMyelin_xyz(:,3) = radMyelin_xyz(:,3)*0.18/0.097;

%% calculate myelin thickness
ef = 3.95;
myelinThickness = MyelinRadius - AxonRadius(idx_Myelin);
myelinThickness(myelinThickness<0) = 0;
myelinThickness = myelinThickness*.097/ef;

figure, scatter3(radMyelin_xyz(:,1),radMyelin_xyz(:,2), radMyelin_xyz(:,3),5,myelinThickness)
colormap jet
% caxis([0 1])
xlabel('x')
axis equal
%% rotate the data
theta = 45;
[nc, R, t] = AxelRot(radMyelin_xyz', theta, [0,1,0], [350,400,1000]);
theta = -20;
[nc, ~, ~] = AxelRot(nc, theta, [0,1,0], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
theta = 70;
[nc, ~, ~] = AxelRot(nc, theta, [0,0,1], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
figure, scatter3(nc(1,:),nc(2,:), nc(3,:),5,corrGRatio')
colormap jet
caxis([0 1])
xlabel('x')
axis equal

%% write rotated data

rotatedCoord = radMyelin_xyz;%nc';
sx = ceil(max(rotatedCoord(:,1))-min(rotatedCoord(:,1))+1);
sy = ceil(max(rotatedCoord(:,2))-min(rotatedCoord(:,2))+1);
sz = ceil(max(rotatedCoord(:,3))-min(rotatedCoord(:,3))+1);

rotatedCoord(:,1) = rotatedCoord(:,1)-min(rotatedCoord(:,1))+1;
rotatedCoord(:,2) = rotatedCoord(:,2)-min(rotatedCoord(:,2))+1;
rotatedCoord(:,3) = rotatedCoord(:,3)-min(rotatedCoord(:,3))+1;

filename = [rt date 'GRATIO.am'];
GU_amiraWriteSpots(filename,rotatedCoord, corrGRatio, [], [])


filename = [rt date 'myelinThickness_um.am'];
GU_amiraWriteSpots(filename,rotatedCoord, myelinThickness, [], [])
