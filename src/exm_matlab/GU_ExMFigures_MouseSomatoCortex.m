% 516
% load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/20171104_CombinedHomorCoordinates_YFPisDAN_516.mat')
% 246
load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/20171128_CombinedHomorCoordinates_YFPisDAN_246.mat')
% load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/20171128_CombinedHomorCoordinates_YFPisDAN_246.mat')

zxRatio = 0.18/0.097;
allData(:,3) = allData(:,3)*zxRatio;
%%

mx = 10580;
my = 66496;
mz = 4488*zxRatio;

cx = mx/2;
cy = my/2;
cz = mz/2;

theta = -32;
clear rotData
[rotData, R, t] = AxelRot(allData', theta, [0,1,0], [cx+1,cy+1,cz+1]);

rotData = rotData';


zl = [4367, 5300]; % 25um - 50um
% zlim = [4600, 4900];

idx = rotData(:,3) > zl(1) & rotData(:,3) < zl(2);
rotDataZlim = rotData(idx,:);
allDensityZlim = allDensity(idx);
%%
figure,
s = 10;
scatter3(rotDataZlim(1:s:end,1),rotDataZlim(1:s:end,2),rotDataZlim(1:s:end,3), 1, allDensityZlim(1:s:end))
% scatter3(rotData(1:s:end,1),rotData(1:s:end,2),rotData(1:s:end,3), 1, allDensity(1:s:end))
% scatter3(allData(1:s:end,1),allData(1:s:end,2),allData(1:s:end,3), 1, allDensity(1:s:end))
xlabel('x')
ylabel('y')
zlabel('z')
% view([90,0]) % yz
view([0,0]) %xz
axis equal
% z range to use : 5300:4367 5300-933

%% get the sum of the sum projected in z

rotDataZlim_round = round(rotDataZlim);
xr = [min(rotDataZlim_round(:,1)), max(rotDataZlim_round(:,1))];
yr = [min(rotDataZlim_round(:,2)), max(rotDataZlim_round(:,2))];
zr = [min(rotDataZlim_round(:,3)), max(rotDataZlim_round(:,3))];


%% hist binning - for z-limted - projected over z
EF = 3.62;
px = 0.097;
bs = 25; % in um
nx = 66414;
ny = 11501;

zl = [4367, 5300]; % 25um - 50um
% zlim = [4600, 4900];

idx = rotData(:,3) > zl(1) & rotData(:,3) < zl(2);
rotDataZlim = rotData(idx,:);
rotDataZlim(:,1) = rotDataZlim(:,1) - xr(1)+1;
rotDataZlim(:,2) = rotDataZlim(:,2) - yr(1)+1;

X = [rotDataZlim(:,1)/EF*px,rotDataZlim(:,2)/EF*px];
% X = [rotData(:,1)/EF*px,rotData(:,2)/EF*px]; % all data

bf = bs/(px/EF); % bin factor

% setup figure and save
ha = setupFigure(1,1, 'AxesWidth', 20, 'AxesHeight', 20,'SameAxes', false,...
    'XSpace', [5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))

hist3(X, [ceil(ny/bf), ceil(nx/bf)],'FaceAlpha',.55, 'EdgeAlpha',0);
% set(ha,'YDir','rev')
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
colormap jet
caxis([1000 10000])
c = colorbar;
c.Color = 'white';
grid on
view([-120,40]) %xz

set(gcf,'color','black')
set(gcf,'DefaultAxesColor','k')
ax = gca;
ax.GridColor = 'white';
ax.Color = 'black';
ax.YColor = 'white';
ax.XColor = 'white';
ax.ZColor = 'white';
daspect([1 1 bs])
zlim([10 1.5*10^4])
xlim([0 300])
ylim([0 1800])
yticks(0:100:1800)
yticklabels(1800:-100:0)
xticks(0:100:300)
ylabel('y-axis (µm)');
xlabel('x-axis (µm)');
title('Homer Counts')
legend(['\color{white} BinSize:XY ' num2str(bs) ' µm; Z: 25µm']);
f0 = gcf();
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/';
% print(f0, '-painters','-dsvg', '-loose',[rt 'HomerCounts_jet_246.svg']);
% print(f0, '-painters','-depsc', '-loose',[rt 'HomerCounts.eps']);

%%

%recalculate density

pef = 3.62;
px = .097; %in um
zAniso = 1; % rotated data
V = 10;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*pef; % radius in expanded samples; in um
bandwidth2 = re/px;


tic
Mdl = KDTreeSearcher(rotData); % all nonDAN synapse positions
toc
tic
IdxKDT = rangesearch(Mdl,rotData,bandwidth2);% Search radius = bandwidth2;
toc
allDensity_recalc = cellfun(@numel, IdxKDT)/V;

%% plot all density data
ha = setupFigure(1,1, 'AxesWidth', 4, 'AxesHeight', 4,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


axes(ha(1))
Didx = allDANidx==1;
NDidx = allDANidx==0;
GU_stairsHistPlot({double(allDensity_recalc(NDidx)),[],[],double(allDensity_recalc(Didx))}, 'ShowECDF', true, 'BinSize',0.1,'CreateFigure', false,...
    'PercentileEnd', 99.9, 'LineWidth',0.75)

title('Density of nc82')
legend(['n = ' num2str(numel(single(allDensity_recalc(NDidx))))], 'Non-DAN', ['n = ' num2str(numel(single(allDensity_recalc(Didx))))], 'DAN')
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 2]);
yyaxis left
% ylim([0 0.15]);
xticks(0:0.5:3)
grid on
ylabel('Relative Frequency');
xlabel('Density (# Homer punca per µm^3)');
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[rt 'Density_246.eps']);

%% hist binning - for all z - projected over y
rotDataZlim = rotData;
xr = [min(rotData(:,1)), max(rotData(:,1))];
yr = [min(rotData(:,2)), max(rotData(:,2))];
zr = [min(rotData(:,3)), max(rotData(:,3))];


EF = 3.62;
px = 0.097;
bs = 3; % in um
nx = 66414;
ny = 11501;

zl = [4367, 5300]; % 25um - 50um
% zlim = [4600, 4900];

% idx = rotData(:,3) > zl(1) & rotData(:,3) < zl(2);

rotDataZlim(:,1) = rotDataZlim(:,1) - xr(1)+1;
rotDataZlim(:,3) = rotDataZlim(:,3) - zr(1)+1;

X = [rotDataZlim(:,1)/EF*px,rotDataZlim(:,3)/EF*px];
% X = [rotData(:,1)/EF*px,rotData(:,2)/EF*px]; % all data

bf = bs/(px/EF); % bin factor

% setup figure and save
ha = setupFigure(1,1, 'AxesWidth', 10, 'AxesHeight', 10,'SameAxes', false,...
    'XSpace', [5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))

hist3(X, [ceil(ny/bf), ceil(nx/bf)],'FaceAlpha',.55, 'EdgeAlpha',0);
set(ha,'yDir','rev')
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
colormap jet
% caxis([1000 34000])
c = colorbar;
c.Color = 'white';
grid on
% view([-70,60]) %
view([0 90])
daspect([1 1 bs*5])
set(gcf,'color','black')
set(gcf,'DefaultAxesColor','k')
ax = gca;
ax.GridColor = 'white';
ax.Color = 'black';
ax.YColor = 'white';
ax.XColor = 'white';
ax.ZColor = 'white';

% zlim([10 3.5*10^4])
% xlim([0 300])
ylim([0 120])
% yticks(0:100:1800)
yticklabels(120:-20:0)
% xticks(0:100:300)
ylabel('z-axis (µm)');
xlabel('x-axis (µm)');
title('Homer Counts')
legend(['\color{white} BinSize:XY ' num2str(bs) ' µm; Y: complete']);
f0 = gcf();
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/';
print(f0, '-painters','-dsvg', '-loose',[rt 'HomerCounts_jet_246_XZ.svg']);

%% homer cleanup punta distribution
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
cutoff = 246;
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerCleanupSupp/';
cd(rt)
load('ex1to5_csizesAfterOTSU.mat')
ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 1 1.5], 'YSpace', [1.5 1 1.5]);
axes(ha(1))
ahs = horzcat(allHomerSizes{:});
histogram(ahs(ahs<cutoff), 0:10:(cutoff-1), 'FaceColor', cfO, 'EdgeColor',cfO);
histogram(ahs(ahs>=cutoff), cutoff:10:3000, 'FaceColor', cfP, 'EdgeColor',cfP);
% line([515,517], [0,10^5])
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
ylabel('# Homer puncta');
xlabel('Size (# Connected Voxels)');
xlim([0 3000])
ylim([1 10^4])
legend(['n= <' num2str(cutoff) '-' num2str(numel(ahs(ahs<cutoff))) '--' 'n= >=' num2str(cutoff) '-' num2str(numel(ahs(ahs>=cutoff)))])
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'HomerPuncaSizeAfterOtsu.eps']);

%% variable
%% homer
%% neurons
%% analysis

%% plot distance of homer to the two segmented neurons

rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/VariableHomerStaining2Neurons/';
cd(rt)
load('seg48_HomerInfo.mat')
load('Set32_HomerInfo.mat')
load('20180703_seg32_OTSU_volume_perim.mat')
load('20180703_seg48_OTSU_volume_perim.mat')
load('20180707_Homer_rp3.mat')
pef = 3.62;
rx = 114.38/1000/pef;
ry = 97/1000/pef;
rz = 155.23/1000/pef;

HomVolDensity_seg32 = numel(HomerAss32_rp3.Volume)/(nVoxSeg32 * rx * ry *rz);
HomVolDensity_seg48 = numel(HomerAss48_rp3.Volume)/(nVoxSeg48 * rx * ry *rz);

sa32 = nVoxSeg32_Perim * ry^2
sa48 = nVoxSeg48_Perim * ry^2


%%

% load('20180703_Seg32_YFPVol_PerimCorrdinates.mat', 'YFP32_VolcoordXYZ')
% load('20180703_Seg48_YFPVol_PerimCorrdinates.mat', 'YFP48_VolcoordXYZ')
load('20180706_Seg32_YFPVol_PerimCorrdinates_params_checked.mat')
load('20180706_Seg48_YFPVol_PerimCorrdinates_params_checked.mat')

coord = YFP32_VolcoordXYZ;
coord(:,1) = coord(:,1) * rx;
coord(:,2) = coord(:,2) * ry;
coord(:,3) = coord(:,3) * rz;

Homer = HomerAll32_rp3.Centroid;  % this is non-yfp-associated homer
Homer(:,1) = Homer(:,1) * rx;
Homer(:,2) = Homer(:,2) * ry;
Homer(:,3) = Homer(:,3) * rz;

tic
Mdl32 = KDTreeSearcher(coord);
toc
% search for the closest V5
tic
[~,D32] = knnsearch(Mdl32, Homer); % in um
toc


coord = YFP48_VolcoordXYZ;
coord(:,1) = coord(:,1) * rx;
coord(:,2) = coord(:,2) * ry;
coord(:,3) = coord(:,3) * rz;

Homer = HomerAll48_rp3.Centroid; % this is non-yfp-associated homer
Homer(:,1) = Homer(:,1) * rx;
Homer(:,2) = Homer(:,2) * ry;
Homer(:,3) = Homer(:,3) * rz;

tic
Mdl48 = KDTreeSearcher(coord);
toc
% search for the closest V5
tic
[~,D48] = knnsearch(Mdl48, Homer); % in um
toc
%%
HomerSel32 = HomerAss32_rp3.Centroid;
HomerSel48 = HomerAss48_rp3.Centroid;
figure,
s = 1100;
% scatter3(data(1:s:end,1),data(1:s:end,2),data(1:s:end,3), 1)
hold on
% scatter3(Homer(:,1),Homer(:,2),Homer(:,3), 1)
% scatter3(coord(1:s:end,1),coord(1:s:end,2),coord(1:s:end,3), 1)
scatter3(HomerSel32(:,1),HomerSel32(:,2),HomerSel32(:,3), 1)
scatter3(HomerSel48(:,1),HomerSel48(:,2),HomerSel48(:,3), 1)
% scatter3(datap32(1:s:end,1),datap32(1:s:end,2),datap32(1:s:end,3), 1)
% scatter3(rotData(1:s:end,1),rotData(1:s:end,2),rotData(1:s:end,3), 1, allDensity(1:s:end))
% scatter3(allData(1:s:end,1),allData(1:s:end,2),allData(1:s:end,3), 1, allDensity(1:s:end))
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (\mum)')
% view([90,0]) % yz
view(3) %xz
axis equal
grid on
%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
% load('20180704_HomerDistance2Neurons.mat')
% load('seg48_HomerInfo.mat', 'HomerAll48_rp3', 'HomerAss48_rp3')
% load('Set32_HomerInfo.mat', 'HomerAll32_rp3', 'HomerAss32_rp3')
% load('20180703_seg32_OTSU_volume_perim.mat')
% load('20180703_seg48_OTSU_volume_perim.mat')

pef = 3.62;
% setup figure and save
ha = setupFigure(2,1, 'AxesWidth', 5, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 1 1.5], 'YSpace', [1.5 1 1.5]);
axes(ha(1))

a32 = numel(HomerAss32_rp3.Volume);
a48 = numel(HomerAss48_rp3.Volume);
t32 = [D32; ones(a32,1)*-0.1];
t48 = [D48; ones(a48,1)*-0.1];
histogram(t32, -0.10:0.05:1, 'FaceColor', cfO, 'FaceAlpha', 0.4, 'EdgeColor', ceO);
histogram(t48, -0.10:0.05:1, 'FaceColor', cfB, 'FaceAlpha', 0.4, 'EdgeColor', ceB);
xlabel('Distance to Neuron (\mum)');
ylabel('# Homer');
set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
legend(['Neuron 1 n=' num2str(numel(t32(t32<=1)))],['Neuron 2 n=' num2str(numel(t48(t48<=1)))])
xlim([-0.1 1]);
xticks([-.1:0.1:1])
ylim([0 2000])



% axes(ha(2))
% histogram(D32(D32>0.1*pef)/pef, 0:0.1:10, 'FaceColor', cfO, 'FaceAlpha', 0.4, 'EdgeColor', ceO);
% histogram(D48(D48>0.1*pef)/pef, 0:0.1:10, 'FaceColor', cfB, 'FaceAlpha', 0.4, 'EdgeColor', ceB);
% xlabel('Distance to Neuron (\mum)');
% ylabel('# Homer');
% set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
% 
% legend(['Neuron 1 n=' num2str(numel(D32(D32>0.1*pef)))],['Neuron 2 n=' num2str(numel(D48(D48>0.1*pef)))])
% xlim([0 3]);
% xticks([0:0.5:3])
% ylim([0 8000])
% 
% 
% d = 0:0.1/pef:3;
% 
% 
% rx = 114.38/1000/pef;
% ry = 97/1000/pef;
% rz = 155.23/1000/pef;
% 
% for i = 1:numel(d)
%     n32(i) = numel(D32(D32<d(i)*pef))/(nVoxSeg32 * rx * ry *rz);
%     n48(i) = numel(D48(D48<d(i)*pef))/(nVoxSeg48 * rx * ry *rz);
%     an32(i) = numel(t32(t32<d(i)))/(nVoxSeg32 * rx * ry *rz);
%     an48(i) = numel(t48(t48<d(i)))/(nVoxSeg48 * rx * ry *rz);
% end
% 
% 
% % plot homer volume density
% axes(ha(3))
% plot(d, an32, 'color', ceO); hold on
% plot(d, an48, 'color', ceB)
% xlabel('Distance from Neuron (\mum)')
% ylabel('All Homer Vol Density')
% xlim([0 3]);
% xticks([0:0.5:3])
% set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
% ylim([0 8])
% 
% axes(ha(5))
% plot(d, an32./an48, 'color', ceR)
% xlabel('Distance from Neuron (\mum)')
% ylabel('Fold Difference')
% xlim([0 3]);
% xticks([0:0.5:3])
% set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
% ylim([0 8])
% 
% axes(ha(4))
% plot(d(d>0.1), n32(d>0.1), 'color', ceO); hold on
% plot(d(d>0.1), n48(d>0.1), 'color', ceB)
% xlabel('Distance from Neuron (\mum)')
% ylabel('non-YFP-Assoc Homer Vol Density')
% xlim([0 3]);
% xticks([0:0.5:3])
% set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
% ylim([0 8])
% 
% axes(ha(6))
% plot(d(d>0.1), n32(d>0.1)./n48(d>0.1), 'color', ceR)
% xlabel('Distance from Neuron (\mum)')
% ylabel('Fold Difference')
% xlim([0 3]);
% xticks([0:0.5:3])
% set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
% ylim([0 8])

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'nHomerDistanceFromNeuron.eps']);

%%
%%
%%
%%
%%

%% analysis of the homer around the two masked neurons
% rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/VariableHomerStaining2Neurons/';
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/VariableHomerStaining2Neurons/';
cd(rt)
load('20180706_Seg32_YFPVol_PerimCorrdinates_params_checked.mat')
load('20180706_Seg48_YFPVol_PerimCorrdinates_params_checked.mat')
% load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/20171128_CombinedHomorCoordinates_YFPisDAN_246.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/20171128_CombinedHomorCoordinates_YFPisDAN_246.mat')
pefx = .114/3.62;
pefy = .097/3.62;
pefz = .155/3.62;
zxRatio = 0.18/0.097;
allData(:,3) = allData(:,3)*zxRatio;

mx = 10580;
my = 66496;
mz = 4488*zxRatio;

cx = mx/2;
cy = my/2;
cz = mz/2;

theta = -32;
clear rotData
[rotData, R, t] = AxelRot(allData', theta, [0,1,0], [cx+1,cy+1,cz+1]);

rotData = rotData';

datar = rotData;
datar(:,1) = (datar(:,1))*pefy;
datar(:,2) = (datar(:,2))*pefy;
datar(:,3) = (datar(:,3)-500*0.155/0.097)*pefy;


% segment 32 coordiante offsets
s32_x = [2816:9763];
s32_y = [29352:36114];
s32_z = [746:2565];

% segment 48 coordiantes
s48_x = [898:6158];
s48_y = [28563:37569];
s48_z = [403:2477];

% perimeter of neurons
datap48 = [YFP48_PerimcoordXYZ];
datap48(:,1) = (datap48(:,1)+s48_x(1))*pefx;
datap48(:,2) = (datap48(:,2)+s48_y(1))*pefy;
datap48(:,3) = (datap48(:,3)+s48_z(1))*pefz;

datap32 = [YFP32_PerimcoordXYZ];
datap32(:,1) = (datap32(:,1)+s32_x(1))*pefx;
datap32(:,2) = (datap32(:,2)+s32_y(1))*pefy;
datap32(:,3) = (datap32(:,3)+s32_z(1))*pefz;

% volume of neurons
datav48 = [YFP48_VolcoordXYZ];
datav48(:,1) = (datav48(:,1)+s48_x(1))*pefx;
datav48(:,2) = (datav48(:,2)+s48_y(1))*pefy;
datav48(:,3) = (datav48(:,3)+s48_z(1))*pefz;

datav32 = [YFP32_VolcoordXYZ];
datav32(:,1) = (datav32(:,1)+s32_x(1))*pefx;
datav32(:,2) = (datav32(:,2)+s32_y(1))*pefy;
datav32(:,3) = (datav32(:,3)+s32_z(1))*pefz;

Homer = datar;  % this is non-yfp-associated homer
Homer = Homer(Homer(:,2)>750 & Homer(:,2)<1075 & Homer(:,3)<120,:);
%% test
figure,
s = 200;
% scatter3(data(1:s:end,1),data(1:s:end,2),data(1:s:end,3), 1)
hold on
scatter3(Homer(1:s:end,1),Homer(1:s:end,2),Homer(1:s:end,3), 1)
scatter3(datap48(1:s:end,1),datap48(1:s:end,2),datap48(1:s:end,3), 1)
scatter3(datap32(1:s:end,1),datap32(1:s:end,2),datap32(1:s:end,3), 1)
% scatter3(rotData(1:s:end,1),rotData(1:s:end,2),rotData(1:s:end,3), 1, allDensity(1:s:end))
% scatter3(allData(1:s:end,1),allData(1:s:end,2),allData(1:s:end,3), 1, allDensity(1:s:end))
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (\mum)')
% view([90,0]) % yz
view(3) %xz
axis equal
grid on
%% distance from the neuron volumes

tic
Mdl32 = KDTreeSearcher(datav32);
toc

tic
Mdl48 = KDTreeSearcher(datav48);
toc

r = 0.0; % Search radius
IdxKDT = rangesearch(Mdl32,Homer,r);
a = cellfun(@numel, IdxKDT);
sum(a)


tic
[~,D32] = knnsearch(Mdl32, Homer, 'k',1); % in um
toc
tic
[~,D48] = knnsearch(Mdl48, Homer, 'k',1); % in um
toc
