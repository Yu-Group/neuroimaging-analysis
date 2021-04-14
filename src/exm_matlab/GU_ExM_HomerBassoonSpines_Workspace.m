rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch0_Analysis_1ch\1008\';
cd(rt)

%% load distance transform skeleton
rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch0_Analysis_1ch\1008\DistMapSkel\';
cd(rt)
fn = 'ch0_1008_clean_DTSK.tif';

[DTSkel,~] = imread_big([rt fn]);
% DTSkel = readtiff([rt fn]);
skel = logical(DTSkel);
%% clean up and genearte nodes
MinBranchLength = 10;
fprintf('Converting Skeleton to Graph, cleaning up and back to Skeleton...')
tic
[~,node,link] = Skel2Graph3D(skel,MinBranchLength);
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
skel = Graph2Skel3D(node,link,w,l,h);
toc

%% ch1 is HOMER (not Basson) - re do analysis:

%% reclean bassoon and find pairs of homer and basson

% imB = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1\ch1_1008_clean.tif']);
% imB = imB-1500;
% writetiff(imB, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1\ch1_1008_clean_v2.tif');
imH = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1_v2\Analysis_1ch\246\ch1_1008_clean_v2_246_clean.tif']);
imB = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch2\ch2_1008_clean.tif']);
%%
SE = strel('sphere', 3);
tic
labelB = bwlabeln(logical(imB));
toc
logIM = logical(imH);
dilogIM = imdilate(logIM, SE);
JuxHOMERIdx = unique(labelB(dilogIM));
JuxHOMERIdx = JuxHOMERIdx(JuxHOMERIdx>0);
JuxHOMERYFP = zeros(size(imB),'logical');
JuxHOMERYFP(ismember(labelB,JuxHOMERIdx)) = 1;

JIB = zeros(size(imB), 'uint16');
JIB(JuxHOMERYFP) = imB(JuxHOMERYFP); %% bassoon that is juxaposed with the homer - for all

labelH = bwlabeln(logical(imH));
dilogBas = imdilate(JuxHOMERYFP, SE);
JuxBASIdx = unique(labelH(dilogBas));
JuxBASOONYFP = zeros(size(imH),'logical');
JuxBASOONYFP(ismember(labelH, JuxBASIdx)) = 1;

JIH = zeros(size(imH), 'uint16');
JIH(JuxBASOONYFP) = imH(JuxBASOONYFP); %% bassoon that is juxaposed with the homer - for all

writetiff(JIB, 'BassonJuxH.tif');
writetiff(JIH, 'HomerJuxB.tif');
%% read cleaned YFP signal
rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch0_Analysis_1ch\1008\';
cd(rt)
fn = 'ch0_1008_clean.tif';
im = readtiff([rt fn]);
%% get yfp associated homer and basson (paird with homer)
SE = strel('sphere', 3);
tic
% [labelH, nh] = bwlabeln(logical(imH));
toc
tic
% logIM = logical(im);
dilogIM = imdilate(logIM, SE);
JuxHOMERIdx = unique(labelH(dilogIM));
JuxHOMERIdx = JuxHOMERIdx(JuxHOMERIdx>0);
JuxHOMERYFP = zeros(size(imH),'logical');
JuxHOMERYFP(ismember(labelH,JuxHOMERIdx)) = 1;

JIH_YFP = zeros(size(imH), 'uint16');
JIH_YFP(JuxHOMERYFP) = JIH(JuxHOMERYFP); %% bassoon that is juxaposed with the homer - for all
writetiff(JIH_YFP, 'YFPass_HomerJuxB.tif');

% labelB = bwlabeln(logical(imB));
dilogBas = imdilate(JuxHOMERYFP, SE);
JuxBASIdx = unique(labelB(dilogBas));
JuxBASOONYFP = zeros(size(JIB),'logical');
JuxBASOONYFP(ismember(labelB, JuxBASIdx)) = 1;

JIB_YFP = zeros(size(JIB), 'uint16');
JIB_YFP(JuxBASOONYFP) = JIB(JuxBASOONYFP); %% bassoon that is juxaposed with the homer - for all

writetiff(JIB_YFP, 'YFPass_BassonJuxH.tif');

toc

%%
%
JIB = readtiff('BassonJuxH.tif');
JIH = readtiff('HomerJuxB.tif');

JIB_YFP = readtiff('YFPass_BassonJuxH.tif');
JIH_YFP = readtiff('YFPass_HomerJuxB.tif');

%% get the centroids
tic
B_rp = regionprops(logical(JIB), 'Centroid'); toc
tic
H_rp = regionprops(logical(JIH), 'Centroid'); toc
tic
BY_rp = regionprops(logical(JIB_YFP), 'Centroid'); toc
tic
HY_rp = regionprops(logical(JIH_YFP), 'Centroid'); toc

save('20180529_HomerBassson_rp.mat', 'B_rp', 'H_rp','BY_rp','HY_rp');
%% get weigthed centroids
tic
B_wc = regionprops3(JIB>0, JIB, {'WeightedCentroid'}); toc
tic
H_wc = regionprops3(JIH>0, JIH, {'WeightedCentroid'}); toc

tic
BY_wc = regionprops3(JIB_YFP>0, JIB_YFP, {'WeightedCentroid'}); toc


tic
HY_wc =regionprops3(JIH_YFP>0, JIH_YFP, {'WeightedCentroid'}); toc

save('20180529_HomerBassson_rp.mat', 'B_rp', 'H_rp','BY_rp','HY_rp', 'B_wc', 'H_wc','BY_wc','HY_wc');

%%

ef = 4.04;
px = .097; %in um
pz = 0.25;
zAniso = pz/px;

rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\';
cd(rt)
load('20180529_HomerBassson_rp.mat')

% H = H_wc;
% B = B_wc;
% HY = HY_wc;
% BY = BY_wc;

% Homer_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),H_rp)]*px/ef;
% Homer_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),H_rp)]*px/ef;
% Homer_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),H_rp)]*zAniso*px/ef;
% 
% Bassoon_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),B_rp)]*px/ef;
% Bassoon_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),B_rp)]*px/ef;
% Bassoon_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),B_rp)]*zAniso*px/ef;
% 
% HomerYFP_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),HY_rp)]*px/ef;
% HomerYFP_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),HY_rp)]*px/ef;
% HomerYFP_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),HY_rp)]*zAniso*px/ef;
% 
% BassoonYFP_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),BY_rp)]*px/ef;
% BassoonYFP_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),BY_rp)]*px/ef;
% BassoonYFP_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),BY_rp)]*zAniso*px/ef;

Homer_centroid(:,1) = H_wc.WeightedCentroid(:,1)*px/ef;
Homer_centroid(:,2) = H_wc.WeightedCentroid(:,2)*px/ef;
Homer_centroid(:,3) = H_wc.WeightedCentroid(:,3)*zAniso*px/ef;

Bassoon_centroid(:,1) = B_wc.WeightedCentroid(:,1)*px/ef;
Bassoon_centroid(:,2) = B_wc.WeightedCentroid(:,2)*px/ef;
Bassoon_centroid(:,3) = B_wc.WeightedCentroid(:,3)*zAniso*px/ef;

HomerYFP_centroid(:,1) = HY_wc.WeightedCentroid(:,1)*px/ef;
HomerYFP_centroid(:,2) = HY_wc.WeightedCentroid(:,2)*px/ef;
HomerYFP_centroid(:,3) = HY_wc.WeightedCentroid(:,3)*zAniso*px/ef;

BassoonYFP_centroid(:,1) = BY_wc.WeightedCentroid(:,1)*px/ef;
BassoonYFP_centroid(:,2) = BY_wc.WeightedCentroid(:,2)*px/ef;
BassoonYFP_centroid(:,3) = BY_wc.WeightedCentroid(:,3)*zAniso*px/ef;


tic
B_Mdl = KDTreeSearcher(Bassoon_centroid);
toc

tic
H_Mdl = KDTreeSearcher(Homer_centroid);
toc

tic
BY_Mdl = KDTreeSearcher(BassoonYFP_centroid);
toc

tic
HY_Mdl = KDTreeSearcher(HomerYFP_centroid);
toc

tic
[~,D_HB] = knnsearch(H_Mdl, Bassoon_centroid, 'k',1);
toc

tic
[~,D_YHB] = knnsearch(HY_Mdl, BassoonYFP_centroid, 'k',1);
toc

%%
ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
histogram(D_HB*1000,0:10:600)
set(ha, 'FontSize', 6);
xlim([0 600])
xticks(0:150:600)
ylabel('# events');
xlabel('Distance (nm)');
legend(['n=' num2str(numel(D_HB))])
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'DistanceFromHomer2ClosestBasoon.eps']);


ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
histogram(D_YHB*1000,0:10:600)
set(ha, 'FontSize', 6);
xlim([0 600])
xticks(0:150:600)
ylabel('# events');
xlabel('Distance (nm)');
legend(['n=' num2str(numel(D_YHB))])
f0 = gcf();

%%
tic
S_B = bwskel(logical(JIB)); toc

tic
S_H = bwskel(logical(JIH)); toc

tic
S_BY = bwskel(logical(JIB_YFP)); toc

tic
S_HY = bwskel(logical(JIH_YFP)); toc

writetiff(uint8(S_B), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\skel_pairedBassoon.tif');
writetiff(uint8(S_H), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\skel_pairedHomer.tif');
writetiff(uint8(S_BY), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\skel_pairedBassoon_proxHomer_YFP.tif');
writetiff(uint8(S_HY), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\skel_pairedHomer_proxYFP.tif');


% get voxel values from the skeletonized data
tic
B_VV = regionprops3(S_B>0, JIB, {'VoxelValues'}); toc
tic
H_VV = regionprops3(S_H>0, JIH, {'VoxelValues'}); toc
tic
BY_VV = regionprops3(S_BY>0, JIB_YFP, {'VoxelValues'}); toc
tic
HY_VV =regionprops3(S_HY>0, JIH_YFP, {'VoxelValues'}); toc

save('20180529_HomerBassson_rp3_VoxelValues.mat', 'B_VV', 'H_VV','BY_VV','HY_VV');

%% interpolate skeletons

tic
JIB_YFP = readtiff('YFPass_BassonJuxH.tif'); toc

tic
iJIB_YFP = GU_interp3DlargeVolume(JIB_YFP, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
writetiff(iJIB_YFP, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\interpolated_pairedBassoon_proxHomer_YFP.tif');
clear JIB_YFP
iJIB_YFP = logical(iJIB_YFP);
tic
iS_BY = bwskel(iJIB_YFP); toc
writetiff(uint8(iS_BY), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\interp_skel_pairedBassoon_proxHomer_YFP.tif');
clear iS_BY iJIB_YFP

tic
JIH_YFP = readtiff('YFPass_HomerJuxB.tif'); toc
tic
iJIH_YFP = GU_interp3DlargeVolume(JIH_YFP, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
writetiff(iJIH_YFP, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\interpolated_pairedHomer_proxYFP.tif');
clear JIH_YFP
iJIH_YFP = logical(iJIH_YFP);
tic
iS_HY = bwskel(logical(iJIH_YFP)); toc

writetiff(uint8(iS_HY), 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\skel\interp_skel_pairedHomer_proxYFP.tif');

%% interpolate all pairs and calculate the major axis length
tic
JIB = readtiff('BassonJuxH.tif'); toc

tic
iJIB = GU_interp3DlargeVolume(JIB, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
writetiff(iJIB, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\interpolated_pairedBassoon_proxHomer.tif');
clear JIB_YFP
iJIB = logical(iJIB);

tic
JIH = readtiff('HomerJuxB.tif'); toc
tic
iJIH = GU_interp3DlargeVolume(JIH, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
writetiff(iJIH, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\interpolated_pairedHomer.tif');
clear JIH
iJIH = logical(iJIH);

% get principle axis lenght from the interpolated skeletonized data
tic
iB_PA = regionprops3(iJIB, {'PrincipalAxisLength', 'Centroid'}); toc
tic
iH_PA = regionprops3(iJIH, {'PrincipalAxisLength', 'Centroid'}); toc

save('20180601_HomerBassson_rp3_PrincipleAxisLength.mat', 'iB_PA', 'iH_PA');
%%
ef = 4.04;
px = .097; %in um
pz = 0.25;
zAniso = pz/px;

Homer_centroid(:,1) = iH_PA.Centroid(:,1)*px/ef;
Homer_centroid(:,2) = iH_PA.Centroid(:,2)*px/ef;
Homer_centroid(:,3) = iH_PA.Centroid(:,3)*px/ef;

Bassoon_centroid(:,1) = iB_PA.Centroid(:,1)*px/ef;
Bassoon_centroid(:,2) = iB_PA.Centroid(:,2)*px/ef;
Bassoon_centroid(:,3) = iB_PA.Centroid(:,3)*px/ef;


tic
B_Mdl = KDTreeSearcher(Bassoon_centroid);
toc

tic
H_Mdl = KDTreeSearcher(Homer_centroid);
toc

tic
[idx,D_HB] = knnsearch(H_Mdl, Bassoon_centroid, 'k',1);
toc

Homer_PA_all(:,1) = iH_PA.PrincipalAxisLength(:,1)*px/ef;
Bassoon_PA(:,1) = iB_PA.PrincipalAxisLength(:,1)*px/ef;
Homer_PA = Homer_PA_all(idx);
save('20180601_HomerBassson_rp3_PrincipleAxisLength.mat', 'iB_PA', 'iH_PA', 'Bassoon_centroid',...
    'Bassoon_PA', 'Homer_centroid', 'Homer_PA', 'Homer_PA_all', 'idx');
%%