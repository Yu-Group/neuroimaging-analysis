% rt = '\\dm11\betziglab\4Stephan\160402_Sample9_YFP_MBP_S35\Stitch_Igor\Analysis\Axon1\Axon1_ch0,ch1,ch0x2+1\Axon1_mask\';
% fn_A = 'Axon1_channel 0 [3519, 522, 5240]_Axon.wrl';
% fn_M = 'Axon1_channel 0 [3519, 522, 5240]_Myelin.wrl';
% cd(rt)
%
% [~,Axon,~] = read_vrml([rt fn_A]);
% [~,Melin,~] = read_vrml([rt fn_M]);
%
%%

figure
scatter3(Axon(1).pts(:,1)*SF,Axon(1).pts(:,2)*SF,Axon(1).pts(:,3)*SF*zAniso, 1)
hold on
pause
scatter3(round(Axonpts(:,1)*SF)/SF,round(Axonpts(:,2)*SF)/SF,round(Axonpts(:,3)*SF)/SF, 1)

hold on
pause
scatter3(round(Myelinpts(:,1)*SF)/SF,round(Myelinpts(:,2)*SF)/SF,round(Myelinpts(:,3)*SF)/SF, 1)

%% skip to section below
clear all
rt = '\\dm11\betziglab\4Stephan\160402_Sample9_YFP_MBP_S35\Stitch_Igor\Analysis\Axon1\Axon1_ch0,ch1,ch0x2+1\Axon1_mask\';
% rt = '/groups/betzig/betziglab/4Stephan/160402_Sample9_YFP_MBP_S35/Stitch_Igor/Analysis/Axon1/Axon1_ch0,ch1,ch0x2+1/Axon1_mask/';
load([rt 'Axon_MyelinData.mat'])
cd(rt)
zAniso = 1; %0.18/0.097;
SF = 10;
Axonpts(:,1) = Axon(1).pts(:,1)*SF;
Axonpts(:,2) = Axon(1).pts(:,2) *SF;
Axonpts(:,3) = Axon(1).pts(:,3)*zAniso*SF;
Myelinpts(:,1) = Melin(1).pts(:,1)*SF;
Myelinpts(:,2) = Melin(1).pts(:,2)*SF;
Myelinpts(:,3) = Melin(1).pts(:,3)*zAniso*SF;
% rt = '/groups/betzig/betziglab/4Stephan/160402_Sample9_YFP_MBP_S35/Stitch_Igor/Analysis/Axon1/Axon1_ch0,ch1,ch0x2+1/Axon1_mask/';
%% skip to section below
xr = [min([Axonpts(:,1); Myelinpts(:,1)]), max([Axonpts(:,1); Myelinpts(:,1)])];
yr = [min([Axonpts(:,2); Myelinpts(:,2)]), max([Axonpts(:,2); Myelinpts(:,2)])];
zr = [min([Axonpts(:,3); Myelinpts(:,3)]), max([Axonpts(:,3); Myelinpts(:,3)])];
nx = ceil(ceil(xr(2))-floor(xr(1)))+1;
ny = ceil(ceil(yr(2))-floor(yr(1)))+1;
nz = ceil(ceil(zr(2))-floor(zr(1)))+1;

% [k,v ]= boundary(Axonpts(:,1)-xr(1)+1,Axonpts(:,2)-yr(1)+1,Axonpts(:,3)-zr(1)+1);

%% skip to below
imAxon = zeros([ny,nx,nz],'logical');
imMyelin = zeros([ny,nx,nz],'logical');

idx = sub2ind([ny,nx,nz], ceil(Axonpts(:,2)-yr(1)+1), ceil(Axonpts(:,1)-xr(1)+1), ceil(Axonpts(:,3)-zr(1)+1));
imAxon(idx) = 1;
idx = sub2ind([ny,nx,nz], floor(Axonpts(:,2)-yr(1)+1), floor(Axonpts(:,1)-xr(1)+1), floor(Axonpts(:,3)-zr(1)+1));
imAxon(idx) = 1;


idx = sub2ind([ny,nx,nz], ceil(Myelinpts(:,2)-yr(1)+1), ceil(Myelinpts(:,1)-xr(1)+1), ceil(Myelinpts(:,3)-zr(1)+1));
imMyelin(idx) = 1;
idx = sub2ind([ny,nx,nz], floor(Myelinpts(:,2)-yr(1)+1), floor(Myelinpts(:,1)-xr(1)+1), floor(Myelinpts(:,3)-zr(1)+1));
imMyelin(idx) = 1;

% % fill holes
se = strel('sphere', 1);
imAxon = imdilate(imAxon,se);
imMyelin = imdilate(imMyelin,se);

parfor_progress(size(imAxon,3));
for z = 1:size(imAxon,3)
    imAxon(:,:,z) = imfill(imAxon(:,:,z), 'holes');
    imMyelin(:,:,z) = imfill(imMyelin(:,:,z), 'holes');
    parfor_progress;
end
parfor_progress(0);

imAxon = imerode(imAxon,se);
imMyelin = imerode(imMyelin,se);

writetiff(uint8(imAxon), [rt 'Filled_Axon.tif']);
writetiff(uint8(imMyelin), [rt 'Filled_Myelin.tif']);

%
tic
dt_Axon = bwdist(~imAxon,'euclidean'); toc
tic
dt_Myelin  = bwdist(~imMyelin,'euclidean');toc
writetiff(single(dt_Axon), [rt 'DistanceTransfrom_Axon.tif']);
writetiff(single(dt_Myelin), [rt 'DistanceTransfrom_Myelin.tif']);

%
imAxon = logical(readtiff([rt 'Filled_Axon.tif']));
imMyelin = logical(readtiff([rt 'Filled_Myelin.tif']));

% mex function
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

%%  read distance transform and skel images %%%%%%%%%%%%%%%% read
clear skel
dt_Axon = readtiff([rt 'DistanceTransfrom_Axon.tif']);
dt_Myelin = readtiff([rt 'DistanceTransfrom_Myelin.tif']);
skel_Axon = readtiff([rt 'Skel_filled_3d_Axon.tif']);
skel_Myelin = readtiff([rt 'Skel_filled_3d_Myelin.tif']);
%% skip section
MinBranchLength = 20;

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

%% read cleaned data
skel_Axon_clean = readtiff([rt 'Skel_filled_3d_Axon_clean.tif']);
skel_Myelin_clean = readtiff([rt 'Skel_filled_3d_Myelin_clean.tif']);
idx_Axon = find(skel_Axon_clean);
[Ya,Xa,Za] = ind2sub(size(skel_Axon_clean),idx_Axon);
idx_Myelin = find(skel_Myelin_clean);
[Ym,Xm,Zm] = ind2sub(size(skel_Myelin_clean),idx_Myelin);
%%
figure, plot(dt_Axon(idx_Axon)*2/SF); hold on; plot(dt_Myelin(idx_Myelin)*2/SF)
figure, plot(dt_Axon(idx_Axon)*2/SF); hold on; plot(dt_Myelin(idx_Axon)*2/SF)
gRatio = (dt_Axon(idx_Axon)*2/SF)./(dt_Myelin(idx_Axon)*2/SF);

centerDist = squareform(pdist(horzcat(Xa, Ya, Za)));


figure, plot(gRatio(25:end-25))
% figure, plot(smooth(dt_Axon(idx_Axon)*2/SF)); hold on; plot(smooth(dt_Myelin(idx_Myelin)*2/SF))

%% calculate 3D g-ration ---- tkhpc36c
rt = 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\G-RatioData\';
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/';
cd(rt)

skel_Axon_clean = readtiff([rt 'Skel_filled_3d_Axon_clean.tif']); % file is isotropic
imAxon = logical(readtiff([rt 'Filled_Axon.tif']));
imMyelin = logical(readtiff([rt 'Filled_Myelin.tif']));
rad_imAxon = readtiff([rt 'radius_imAxon.tif']);
rad_imMyelin = readtiff( [rt 'radius_imMyelin.tif']);
% % distance from the Axon Center
% % tic
% % dt_AxonCenter = bwdist(skel_Axon_clean,'euclidean'); toc
% % writetiff(single(skel_Axon_clean), [rt 'DistanceTransform_from_Skel_filled_3d_Axon_clean.tif']);
% dt_AxonCenter = readtiff([rt 'DistanceTransform_from_Skel_filled_3d_Axon_clean.tif']);
%
% % perimeter of Axon and Myelin
% perim_imAxon = bwperim(imAxon);
% perim_imMyelin = bwperim(imMyelin);
% perim_imAxon(:,:,[1, end]) = 0;
% perim_imMyelin(:,:,[1, end]) = 0;
% writetiff(uint8(perim_imAxon), [rt 'perim_imAxon_clean.tif']);
% writetiff(uint8(perim_imMyelin), [rt 'perim_imMyelin_clean.tif']);
%
% % radius of Axon and Myelin Perim from Axon center
% rad_imAxon = dt_AxonCenter .* single(perim_imAxon);
% rad_imMyelin = dt_AxonCenter .* single(perim_imMyelin);
% writetiff(uint8(rad_imAxon), [rt 'radius_imAxon.tif']);
% writetiff(uint8(rad_imMyelin), [rt 'radius_imMyelin.tif']);


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
%% calculate myelin thickness
ef = 3.95;
myelinThickness = MyelinRadius - AxonRadius(idx_Myelin);
myelinThickness(myelinThickness<0) = 0;
myelinThickness = myelinThickness*.097/ef;
%% rotate the data
theta = 120;
[nc, R, t] = AxelRot(radMyelin_xyz', theta, [1,0,0], [350,400,1000]);
theta = -20;
[nc, ~, ~] = AxelRot(nc, theta, [0,1,0], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
theta = 70;
[nc, ~, ~] = AxelRot(nc, theta, [0,0,1], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
figure, scatter3(nc(1,:),nc(2,:), nc(3,:))
xlabel('x')
axis equal

%% write rotated data

rotatedCoord = nc';
sx = ceil(max(rotatedCoord(:,1))-min(rotatedCoord(:,1))+1);
sy = ceil(max(rotatedCoord(:,2))-min(rotatedCoord(:,2))+1);
sz = ceil(max(rotatedCoord(:,3))-min(rotatedCoord(:,3))+1);

rotatedCoord(:,1) = rotatedCoord(:,1)-min(rotatedCoord(:,1))+1;
rotatedCoord(:,2) = rotatedCoord(:,2)-min(rotatedCoord(:,2))+1;
rotatedCoord(:,3) = rotatedCoord(:,3)-min(rotatedCoord(:,3))+1;

filename = 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\G-RatioData\GRATIO_Rotated.am';
GU_amiraWriteSpots(filename,rotatedCoord, corrGRatio, [], [])


filename = 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\G-RatioData\myelinThickness_um_Rotated.am';
GU_amiraWriteSpots(filename,rotatedCoord, myelinThickness, [], [])
%%
nc = rotatedCoord';
figure, scatter3(nc(1,:),nc(2,:), nc(3,:))
xlabel('x')
axis equal
%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
ha = setupFigure(1,2, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [2 2.0 1.5], 'YSpace', [0.85 1 1]);
set(ha, 'FontSize', 6);

xr = [100,2100];
AR = AxonRadius(idx_Myelin);
[~,unique_idx]= unique(idx_Myelin);
AR = AR(rotatedCoord(unique_idx,1)>xr(1) & rotatedCoord(unique_idx,1)<xr(2))*0.097/ef*1000;
MR = MyelinRadius(rotatedCoord(:,1)>xr(1) & rotatedCoord(:,1)<xr(2))*0.097/ef*1000;

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
histogram(corrGRatio(rotatedCoord(:,1)>xr(1) & rotatedCoord(:,1)<xr(2)), 25, 'EdgeColor', ceG, 'FaceColor', cfG, 'FaceAlpha',1)
xlabel('G-Ratio')
% yyaxis left
ylabel('Counts')
set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
xlim([0.2 1]);
xticks([0:0.2:1])

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' PLots_distancetoAxon_GRatio.eps']);
%%
ha = setupFigure(1,1, 'AxesWidth', 7, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%%
figure

scatter3(rotatedCoord(:,1),rotatedCoord(:,2),rotatedCoord(:,3), 1, corrGRatio,'.'); axis equal; colormap jet

axis equal
view([-40,20])
% daspect([1 1 1])
%
% xlim([0, 200])
% ylim([0, 10])
% zlim([0 20]);
xlabel('x')
%%
% c = colorbar;
% c.Color = 'black';

set(gcf,'color','black')
set(gcf,'DefaultAxesColor','k')
ax = gca;
ax.GridColor = 'black';%
ax.Color = 'black';
ax.YColor = 'black';%
ax.XColor = 'black';%
ax.ZColor = 'black';%

axis vis3d

% Clean it up

axis off
l = light('Position',[-50 -15 29])
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
% colorbar EastOutside

%%
% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.
y = radMyelin_xyz(:,2)*0.097;
x = radMyelin_xyz(:,1)*0.097;
z = radMyelin_xyz(:,3)*0.097;
c = corrGRatio;

tri = delaunay(x,y,z);
plot3(x,y,z,'.')

% Plot it with TRISURF

h = trisurf(tri, x, y, z, c);
axis vis3d

% Clean it up

axis off
l = light('Position',[-50 -15 29])
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
colorbar EastOutside


%%
%%
%%
%% March 27, 2018 ||||| reprocessing the data using the new stitched data
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/Raw/Axon1_Newstitch/';
fnA = 'Axon1_channel 0 [3519, 522, 5240]_newstitch.tif';
fnM = 'Axon1_channel 1 [3519, 522, 5240]_newstitch.tif';
cd(rt)
%%
imA = readtiff(fnA);
imM = readtiff(fnM);
imA = imA(:,:,570:1770);
imM = imM(:,:,570:1770);

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

bwl_A = bwlabeln(simA);

mask_A = zeros(size(bwl_A),'logical');
mask_A(bwl_A==1) = 1;

se = strel('sphere',15);
dmask_A = imdilate(mask_A, se);
mask_M = zeros(size(bwl_A),'logical');
mask_M(dmask_A) = logical(simM(dmask_A));

writetiff(uint8(mask_M), ['Masked_15px' fnM]);
writetiff(uint8(mask_A), ['Masked' fnA]);

%% cleanCurated
curatedM = readtiff([date 'curated_Masked_15pxAxon1_channel 1 [3519, 522, 5240]_newstitch.tif']);
curatedM = bwlabeln(curatedM);
cleanCuratedM = zeros(size(curatedM),'logical');
cleanCuratedM(curatedM==1) = 1;
writetiff(uint8(cleanCuratedM), [date 'curated_Masked_15pxAxon1_channel 1 [3519, 522, 5240]_newstitch.tif']);

%%
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/Raw/Axon1_Newstitch/';
cd(rt)
imM = readtiff('28-Mar-2018curated_Masked_15pxAxon1_channel 1 [3519, 522, 5240]_newstitch.tif');
imA = readtiff('MaskedAxon1_channel 0 [3519, 522, 5240]_newstitch.tif');

imM = imM(:,:,1:1104);
imA = imA(:,:,1:1104);
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


%% calculate myelin thickness
ef = 3.95;
myelinThickness = MyelinRadius - AxonRadius(idx_Myelin);
myelinThickness(myelinThickness<0) = 0;
myelinThickness = myelinThickness*.097/ef;
%% rotate the data
theta = 120;
[nc, R, t] = AxelRot(radMyelin_xyz', theta, [1,0,0], [350,400,1000]);
theta = -20;
[nc, ~, ~] = AxelRot(nc, theta, [0,1,0], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
theta = 70;
[nc, ~, ~] = AxelRot(nc, theta, [0,0,1], [(max(nc(1,:))-min(nc(1,:)))/2,(max(nc(2,:))-min(nc(2,:)))/2,(max(nc(3,:))-min(nc(3,:)))/2]);
figure, scatter3(nc(1,:),nc(2,:), nc(3,:))
xlabel('x')
axis equal

%% write rotated data

rotatedCoord = nc';
sx = ceil(max(rotatedCoord(:,1))-min(rotatedCoord(:,1))+1);
sy = ceil(max(rotatedCoord(:,2))-min(rotatedCoord(:,2))+1);
sz = ceil(max(rotatedCoord(:,3))-min(rotatedCoord(:,3))+1);

rotatedCoord(:,1) = rotatedCoord(:,1)-min(rotatedCoord(:,1))+1;
rotatedCoord(:,2) = rotatedCoord(:,2)-min(rotatedCoord(:,2))+1;
rotatedCoord(:,3) = rotatedCoord(:,3)-min(rotatedCoord(:,3))+1;

filename = [rt date 'GRATIO_Rotated.am'];
GU_amiraWriteSpots(filename,rotatedCoord, corrGRatio, [], [])


filename = [rt date 'myelinThickness_um_Rotated.am'];
GU_amiraWriteSpots(filename,rotatedCoord, myelinThickness, [], [])
