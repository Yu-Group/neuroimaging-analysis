%% original label calculation
cd('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain')
% cd('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain')
load('20171025_CombinedSynapses_withParams.mat')
% load('20171025_CombinedSynapses_labelidx.mat')
load('20171231_CombinedSynapses_37labelsidxfixed.mat')
load('20171101_allDensityRecalc.mat')
%%
% mask = uint16(readtiff('CombinedMask.tif'));
mask = uint8(readtiff('CombinedMask_37labels_fixedv2.tif'));

mx = 15055;
my = 27964;
mz = 6647;

% mask dimensions
sy = 1747;
sx = 940;
sz = 1661;

% mask downsampled
fy = my/sy; %16.0155;%
fx = mx/sx; %16.032;
fz = mz/sz;

x = ceil(allData(:,1)/fx);
y = ceil(allData(:,2)/fy);
z = ceil(allData(:,3)/fz);
[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
labelidx = interp3(X,Y,Z,mask,x,y,z,'nearest');
toc
%% load data
cd('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain')
load('20171025_CombinedSynapses_withParams.mat')
% load('20171025_CombinedSynapses_labelidx.mat')
load('20171231_CombinedSynapses_37labelsidxfixed.mat')
%% volume of each label
ef = 4.09;
px = .097; %in um
zAniso = 0.18/px;
clear labelVol
uidx = unique(mask(:));
n = 1;
for i = 2:numel(uidx)
    labelVol(n) = sum(mask(:)==uidx(i))*zAniso*(px/ef)^3; % BINNED volume in um3 pre-expansion corrected for zAniso
    n = n+1;
end
%labelVol = labelVol*fx*fy*fz; % accounting for the binning difference -- done below
%%
% ln = {'EB',...
%     'FB',...
%     'PB',...
%     'MB-lobe',...
%     'NO Left',...
%     'NO Right',...
%     'Calyx Left',...
%     'Calyx Right',...
%     'Lobulla Left',...
%     'Lobulla Right',...
%     'Lobulla-plate Left',...
%     'Lobulla-plate Right',...
%     'Medulla Left',...
%     'Medulla Right',...
%     'Rest'}

ln = {'L Medulla',...
    'R Medulla',...
    'L Lobulla-Plate',...
    'R Lobulla-Plate',...
    'L Lobulla',...
    'R Lobulla',...
    'L OTU',...
    'R OTU',...
    'L VLPR',...
    'R VLPR',...
    'L LH',...
    'R LH',...
    'L Calyx',...
    'R Calyx',...
    'L MB-a1',... % remove
    'R MB-a1',... % remove
    'L MB-\alpha^{\prime}2\alpha2',...
    'R MB-\alpha^{\prime}2\alpha2',...
    'L MB-\alpha3',...
    'R MB-\alpha3',...
    'L MB-\alpha^{\prime}3',...
    'R MB-\alpha^{\prime}3',...
    'L MB-b"2',...% remove
    'R MB-b"2',...% remove
    'L MB-\gamma1',...
    'R MB-\gamma1',...
    'L ATL',...
    'R ATL',...
    'PB',...
    'EB',...
    'FB',...
    'L NO',...
    'R NO',...
    'L LAL',...
    'R LAL',...
    'L SP',...
    'R SP'};
lidx = [1:14,17:22,25:37];
ln = ln(lidx);
rmidx = [15,16,23,24];
c = categorical(ln');
%% recalculate density

ef = 4.09;
px = .097; %in um
zAniso = 0.18/px;
V = 1;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % radius in expanded samples; in um
bandwidth2 = re/px;


allData_zcorr(:,1) = allData(:,1);
allData_zcorr(:,2) = allData(:,2);
allData_zcorr(:,3) = allData(:,3)*zAniso;

tic
Mdl = KDTreeSearcher(allData_zcorr); % all nonDAN synapse positions
toc
tic
IdxKDT = rangesearch(Mdl,allData_zcorr,bandwidth2);% Search radius = bandwidth2;
toc
allDensity_recalc = cellfun(@numel, IdxKDT)/V;
%%
allDensity = allDensity_recalc;
%%
ha = setupFigure(7,5, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);
n = 1;
for i = lidx%1:numel(uidx)-1
    axes(ha(n))
    Didx = labelidx == i & allDANidx==1;
    NDidx = labelidx == i & allDANidx==0;
    GU_stairsHistPlot({double(allDensity(NDidx)),[],[],double(allDensity(Didx))}, 'ShowECDF', true, 'BinSize',1,'CreateFigure', false,...
        'PercentileEnd', 99.9, 'LineWidth',0.75)
%     ax = gca;
%     ax.BoxStyle = 'full';
    box on
    title(ln{n})
    %     legend(['n = ' num2str(numel(single(allDensity(NDidx))))], 'Non-DAN', ['n = ' num2str(numel(single(allDensity(Didx))))], 'DAN')
    set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
    xlim([0 25]);
    yyaxis left
    ylim([0 0.2]);
    xticks(0:5:25)
    grid off
    %     if i == any(1:5:40)
    ylabel('Relative Frequency','color', 'k');
    %     end
    %     if ~isEven(i)
    xlabel('Density (# nc82 punca per µm^3)','color', 'k');
    %     end
    n = n+1;
end
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' DensityPlots.eps']);
%% distance - DAN 2 DAN
tic
fData(:,1) = allData(:,1)*0.0237;
fData(:,2) = allData(:,2)*0.0237;
fData(:,3) = allData(:,3)*0.044;
Mdl = KDTreeSearcher(fData(allDANidx,:));
toc
tic
[~,D] = knnsearch(Mdl, fData(allDANidx,:), 'k',2);
toc
D2D = D(:,2);

% distance - any 2 any
tic
fData(:,1) = allData(:,1)*0.0237;
fData(:,2) = allData(:,2)*0.0237;
fData(:,3) = allData(:,3)*0.044;
Mdl = KDTreeSearcher(fData);
toc
tic
[~,D] = knnsearch(Mdl, fData, 'k',2);
toc
D = D(:,2);
%% plots for the 37 different areas

ha = setupFigure(7,5, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);
n = 1;
for i = lidx%1:numel(uidx)-1
    axes(ha(n))
    Didx = labelidx == i & allDANidx==1;
    NDidx = labelidx == i & allDANidx==0;
    D2Didx = labelidx(allDANidx==1) == i;
    GU_stairsHistPlot({double(D(NDidx))*1000,[],[],single(D(Didx))*1000,double(D2D(D2Didx))*1000}, 'ShowECDF', true, 'BinSize',20,'CreateFigure', false,...
        'PercentileEnd', 99.9, 'LineWidth',0.75)
    ax = gca;
    ax.BoxStyle = 'full';
    box on
    title(ln{n})
    legend(['n = ' num2str(numel(D(NDidx)))], 'Non-DAN', ['n = ' num2str(numel(D(Didx)))], 'DAN2Any',...
        ['n = ' num2str(numel(D2D(D2Didx)))], 'DAN2DAN')
    legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
    xlim([0 1000]);
    yyaxis left
    ylim([0 0.008]);
    xticks(0:200:1000)
    grid off
    %     if i == 1
    ylabel('Relative Frequency', 'color', 'k');
    %     end
    %     if ~isEven(i)
    xlabel('Distance to closest nc82 punca (nm)', 'color', 'k');
    %     end
    n = n+1;
end
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' DistancePlots.eps']);


%% plots for all points

ha = setupFigure(2,1, 'AxesWidth', 4, 'AxesHeight', 4,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


axes(ha(1))
Didx = allDANidx==1;
NDidx = allDANidx==0;
GU_stairsHistPlot({double(allDensity(NDidx)),[],[],double(allDensity(Didx))}, 'ShowECDF', true, 'BinSize',1,'CreateFigure', false,...
    'PercentileEnd', 99.9, 'LineWidth',0.75)

title('Density of nc82')
legend(['n = ' num2str(numel(single(allDensity(NDidx))))], 'Non-DAN', ['n = ' num2str(numel(single(allDensity(Didx))))], 'DAN')
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 25]);
yyaxis left
ylim([0 0.15]);
xticks(0:10:40)
grid on
ylabel('Relative Frequency');
xlabel('Density (# nc82 punca per µm^3)');

axes(ha(2))
Didx = allDANidx==1;
NDidx = allDANidx==0;

GU_stairsHistPlot({double(D(NDidx))*1000,[],[],double(D(Didx))*1000,single(D2D)*1000}, 'ShowECDF', true, 'BinSize',20,'CreateFigure', false,...
    'PercentileEnd', 99.9, 'LineWidth',0.75)
title('Distance to closest nc82 puncta')
legend(['n = ' num2str(numel(D(NDidx)))], 'Non-DAN', ['n = ' num2str(numel(D(Didx)))], 'DAN2Any',...
    ['n = ' num2str(numel(D2D))], 'DAN2DAN')
legend('boxoff')
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000]);
yyaxis left
ylim([0 0.008]);
xticks(0:250:1000)
grid on
ylabel('Relative Frequency');
xlabel('Distance (nm)');

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' DensityPlots_forALLdata_square.eps']);

%% disptance plots for DAN2DAN

ha = setupFigure(1,15, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

for i = 1:15
    axes(ha(i))
    D2Didx = labelidx(allDANidx==1) == i;
    %     NDidx = labelidx(allDANidx==0) == i;
    GU_stairsHistPlot({[],[],[],[],single(D(Didx))*1000}, 'ShowECDF', true, 'BinSize',20,'CreateFigure', false,...
        'PercentileEnd', 99.9)
    
    title(ln{i})
    %     legend(['n = ' num2str(numel(single(allDensity(NDidx))))], 'Non-DAN', ['n = ' num2str(numel(single(allDensity(Didx))))], 'DAN')
    % legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica')
    xlim([0 1000]);
    yyaxis left
    ylim([0 0.008]);
    xticks(0:250:1000)
    grid on
    if i == 1
        ylabel('Relative Frequency');
    end
    if ~isEven(i)
        xlabel('Distance to closest nc82 punca');
    end
end
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',['DistancePlots_DAN2DAN.eps']);

%% collect n's for each of the category

for i = 1:numel(uidx)-1
    ntotalSyn(i) = sum(labelidx == i);
    nDANSyn(i) = sum(labelidx == i & allDANidx == 1);
    nDANSyn_incPeriph(i) = sum(labelidx == i & allDANsynIDX_includePeriph == 1);
end
%%
%load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/20180111_synapseCount_vol.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/20180111_synapseCount_vol.mat')
clear nSyn_normVol
labelVol = labelVol*16.016*16.0069*4.0018;
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
% c = reordercats(c, ln(:));
x = [0:32];
x = x + cumsum([repmat([1, 0],1,12), 1, 1, 1, repmat([1, 0],1,3)]);
ha = setupFigure(1,1, 'AxesWidth', 16, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

nSyn_normVol_percentDAN = nDANSyn./ntotalSyn * 100;

% nSyn_normVol(:,1) = (nDANSyn'./labelVol');
% nSyn_normVol(:,2) = (nDANSyn_incPeriph'./labelVol') - nSyn_normVol(:,1);
% nSyn_normVol(:,3) = (ntotalSyn'./labelVol') - nSyn_normVol(:,2);
% 
% b = bar(x,(nSyn_normVol(lidx,:)),0.5, 'stacked');
% ylabel('# nc82 / \mum^3')

hAx=gca;             % handles to current figure, axes
hAx.XTick = x;
hAx.XTickLabel= ln;
hAx.XTickLabelRotation = 45;

nSyn(:,1) = (nDANSyn');
nSyn(:,2) = (nDANSyn_incPeriph') - nSyn(:,1);
nSyn(:,3) = (ntotalSyn') - nSyn(:,2);
b = bar(x,(nSyn(lidx,:)),0.5, 'stacked');
ylabel('# nc82')

grid off
% ylim([0.01 30])
set(gca, 'YScale', 'log')

b(1).FaceColor = 'flat';
b(1).CData = cfG;
b(1).EdgeColor = ceG;
b(2).FaceColor = 'flat';
b(2).CData = cfB;
b(2).EdgeColor = ceB;
b(3).FaceColor = 'flat';
b(3).CData = cfR;
b(3).EdgeColor = ceR;

ax = gca;
ax.BoxStyle = 'full';
box on
grid on
ax1_pos = ax.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
yyaxis right
plot(x, nSyn_normVol_percentDAN(lidx),'-', 'color','g', 'linewidth',.1,'Markersize',10)
xlim([0 51])
ylim([0 12])
ylabel('% DAN', 'FontSize', 6)
yyaxis left
yticks(500)
xticks(500)

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' Maskednc82CountDensity_log.eps']);
% print(f0, '-dsvg', '-loose', [date ' Maskednc82CountDensity_log.svg']);
% print(f0, '-painters','-depsc', '-loose',[date ' Maskednc82Count.eps']);
%%
nAll = [nDANSyn' nDANSyn_incPeriph' ntotalSyn'];
ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);
axes(ha(1))
bar(nAll, 'stacked')
set(gca, 'YScale', 'log')

%% surface area and volume
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain';
rt = 'E:\Gokul\ExM_Flybrain_SA_Vol';
cd(rt)
mask = uint8(readtiff([rt filesep 'CombinedMask_37labels_fixedv2.tif']));
% rt2 = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol';
load([rt filesep '20180214_interpSurfAreaCoord_uint16.mat'])

mx = 15055;
my = 27964;
mz = 6647;
zA = 0.18/0.097;
% mask dimensions
sy = 1747;
sx = 940;
sz = 1661;

% mask downsampled
fy = my/sy; %16.0155;%
fx = mx/sx; %16.032;
fz = mz/sz;

x = ceil(alliSACoordXYZ(:,1)/fx);
y = ceil(alliSACoordXYZ(:,2)/fy);
z = ceil(alliSACoordXYZ(:,3)/fz/zA);
[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
SA_labelidx = interp3(X,Y,Z,mask,x,y,z,'nearest');
toc

tic
load(['20180214_NONinterpVolCoord_uint16.mat'])
toc
mx = 15055;
my = 27964;
mz = 6647;
zA = 0.18/0.097;
% mask dimensions
sy = 1747;
sx = 940;
sz = 1661;

% mask downsampled
fy = my/sy; %16.0155;%
fx = mx/sx; %16.032;
fz = mz/sz;

x = single(allnonInterpVolCoordXYZ(:,1)/fx);
y = single(allnonInterpVolCoordXYZ(:,2)/fy);
z = single(allnonInterpVolCoordXYZ(:,3)/fz/zA);


[X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
tic
Vol_labelidx = interp3(X,Y,Z,mask,x,y,z,'nearest');
toc
%% load surface area labels
cd('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/')
load('SurfaceAreaLabels.mat')
load('VolumeLabels_nonInterp.mat')

%%
uidx = unique(SA_labelidx(:));

ef = 4.09;
px = .097; %in um
zAniso = 0.18/px;
clear labelVol

label_Volume = zeros(1, numel(uidx));
label_SurfArea = zeros(1, numel(uidx));
n = 1;
for i = 1:numel(uidx)
    label_Volume(i) = sum(Vol_labelidx(:)==uidx(i))*zAniso*(px/ef)^3; % volume in um3 pre-expansion corrected for zAniso
    label_SurfArea(i) = sum(SA_labelidx(:)==uidx(i))*(px/ef)^2; % volume in um3 pre-expansion corrected for zAniso
    n = n+1;
end
%%
clear nSyn_normVol labelVol
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/20180111_synapseCount_vol.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/LabledRegions_SurfArea_Vol.mat')

% load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/20180111_synapseCount_vol.mat')
% load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/WholeFlyBrain/SurfAreaVol/LabledRegions_SurfArea_Vol.mat')

labelVol = labelVol*16.016*16.0069*4.0018; % accounting for the binning difference
DANVolume = label_Volume(2:end);
DANSurfArea = label_SurfArea(2:end);
%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
% c = reordercats(c, ln(:));
x = [0:32];
x = x + cumsum([repmat([1, 0],1,12), 1, 1, 1, repmat([1, 0],1,3)]);
ha = setupFigure(1,1, 'AxesWidth', 16, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

clear nSyn_normSA nSyn_normVol

nSyn_normVol(:,1) = (nDANSyn'./DANVolume');
nSyn_normVol(:,2) = (nDANSyn_incPeriph'./DANVolume') - nSyn_normVol(:,1);
b = bar(x,(nSyn_normVol(lidx,:)),0.5, 'stacked');
ylabel('# nc82 / DAN Volume \mum^3')

% nSyn_normSA(:,1) = (nDANSyn'./DANSurfArea') ;
% nSyn_normSA(:,2) = (nDANSyn_incPeriph'./DANSurfArea') - nSyn_normSA(:,1);
% b = bar(x,(nSyn_normSA(lidx,:)),0.5, 'stacked');
% ylabel('# nc82 / DAN Surface Area \mum^2')

hAx=gca;             % handles to current figure, axes
hAx.XTick = x;
hAx.XTickLabel= ln;
hAx.XTickLabelRotation = 45;


grid off
% ylim([0.01 30])
set(gca, 'YScale', 'log')

b(1).FaceColor = 'flat';
b(1).CData = cfG;
b(1).EdgeColor = ceG;
b(2).FaceColor = 'flat';
b(2).CData = cfG;
b(2).EdgeColor = ceG;


ax = gca;
ax.BoxStyle = 'full';
box on

ax1_pos = ax.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
yyaxis right
plot(x, DANVolume(lidx), 'color','g', 'linewidth',1)
% plot(x, DANSurfArea(lidx), 'color','g', 'linewidth',1)
set(gca, 'YScale', 'log')
% xlim([0 51])
% ylim([0 12])
ylabel('DAN Volume \mum^3', 'FontSize', 6)
% ylabel('DAN Surface Area \mum^2', 'FontSize', 6)

yyaxis left
yticks(500)
xticks(500)

f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date ' Maskednc82_DANVol.eps']);

%% fraction DAN volume of region

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
x = [0:32];
x = x + cumsum([repmat([1, 0],1,12), 1, 1, 1, repmat([1, 0],1,3)]);
ha = setupFigure(1,1, 'AxesWidth', 16, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

fracVol = DANVolume./labelVol;
b = bar(x, fracVol(lidx), 0.5)
set(gca, 'YScale', 'log')
b(1).FaceColor = 'flat';
b(1).CData = cfG;
b(1).EdgeColor = ceG;
ylabel('Fractional DAN Volume', 'FontSize', 6)


yyaxis right
plot(x, DANVolume(lidx), 'color',ceO, 'linewidth',1)
ylabel('DAN Volume \mum^3', 'FontSize', 6)
ylim([1 10^5])
set(gca, 'YScale', 'log')
hAx=gca;             % handles to current figure, axes
hAx.XTick = x;
hAx.XTickLabel= ln;
hAx.XTickLabelRotation = 45;

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' FractionalDANVolume_log.eps']);

%% compile summary table

flybrainStats = zeros(33,11);
flybrainStats(:,1) = labelVol(lidx)'; % domain volume
flybrainStats(:,2) = DANVolume(lidx)'; % DAN volume
flybrainStats(:,3) = DANSurfArea(lidx)'; % DAN Surf Area
flybrainStats(:,4) = ntotalSyn(lidx)'; % # nc82 total
flybrainStats(:,5) = nDANSyn_incPeriph(lidx)' - nDANSyn(lidx)'; % # nc82 DAN vicinity
flybrainStats(:,6) = nDANSyn(lidx)'; % # nc82 DAN associated
flybrainStats(:,7) =  ntotalSyn(lidx)'./labelVol(lidx)'; % # nc82 domain density
flybrainStats(:,8) =  nDANSyn(lidx)'./labelVol(lidx)'; % # DAN ass nc82 domain density
flybrainStats(:,9) =  nDANSyn(lidx)'./ntotalSyn(lidx)' * 100; % percent nc82 DAN ass
flybrainStats(:,10) =  nDANSyn(lidx)'./DANVolume(lidx)'; % # DAN ass nc82 DAN density
flybrainStats(:,11) =  nDANSyn(lidx)'./DANSurfArea(lidx)'; % # DAN ass nc82 DAN surface area

%% plot MB-alpha3 stereotypy
cd('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/FlyStereotypy_MB')
load('MBStereotypy_num_nc82_masked_withMBfromFullBrain.mat')

nc82_left = [lMB, MB_nc82([1,3,4])];
nc82_right = [rMB, MB_nc82([2,5,6])];
t = 'MB-\alpha3';

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO] = GU_colorCodes;
xlabels = {'left', 'right'};
fset = loadFigureSettings('print');
ha = setupFigure(1,1, 'AxesWidth', 2.25, 'AxesHeight', 3.5,'SameAxes', true,...
    'XSpace', [1.5 0.85 0.5], 'YSpace', [2.5 0.5 0.5]);
set(ha, 'FontSize', 6);

cv = [ceR; ceB;];
cf = [cfR; cfB;];
opts = {'LineWidth', 1, 'GroupDistance', 0, 'BarWidth', 0.6, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BorderWidth', 0.5, 'Angle', 0, 'DetectOutliers', false,...
    'FaceColor', cf, 'EdgeColor', cv, 'ErrorbarColor', cv};

title(ha(1), t, fset.lfont{:}, 'FontSize', 6);
ylabel(ha(1), '# nc82', fset.lfont{:}, 'FontSize', 6);
boxplot2({[nc82_left], [nc82_right]}, [], opts{:}, 'Parent', ha(1));
hold on
plot(ha(1), 1, [nc82_left], '.', 'Color', .6*ceR, 'MarkerSize', 5);
plot(ha(1), 2, [nc82_right], '.', 'Color', .6*ceB, 'MarkerSize',5);
rotateXTickLabels(ha(1), 'Angle', 45, 'AdjustFigure', false);

% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date ' number_nc82_left_vs_right.eps']);

f0 = gcf();
print(f0, [date 'number_nc82_left_vs_right.pdf'],'-dpdf', '-r0')

% plot the same figure in the top, but merged
t = 'MB-\alpha3';

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO] = GU_colorCodes;
xlabels = {'left & right'};
fset = loadFigureSettings('print');
ha = setupFigure(1,1, 'AxesWidth', 1.25, 'AxesHeight', 3.5,'SameAxes', true,...
    'XSpace', [1.5 0.85 0.5], 'YSpace', [2.5 0.5 0.5]);
set(ha, 'FontSize', 6);

cv = [ceG; ];
cf = [cfG; ];
opts = {'LineWidth', 1, 'GroupDistance', 0, 'BarWidth', 0.6, 'XTickLabel', xlabels,...
    'AdjustFigure', false, 'BorderWidth', 0.5, 'Angle', 0, 'DetectOutliers', false,...
    'FaceColor', cf, 'EdgeColor', cv, 'ErrorbarColor', cv};

title(ha(1), t, fset.lfont{:}, 'FontSize', 6);
ylabel(ha(1), '# nc82', fset.lfont{:}, 'FontSize', 6);
boxplot2({[nc82_left, nc82_right]}, [], opts{:}, 'Parent', ha(1));
hold on
plot(ha(1), 1, [nc82_left, nc82_right], '.', 'Color', .6*ceG, 'MarkerSize', 5);
rotateXTickLabels(ha(1), 'Angle', 45, 'AdjustFigure', false);

% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date ' number_nc82_combined.eps']);
f0 = gcf();
print(f0, [date 'number_nc82_combined.pdf'],'-dpdf', '-r0')

% plot the same figure in the top, but merged
t = 'MB-\alpha3';

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO] = GU_colorCodes;
xlabels = {'left & right'};
fset = loadFigureSettings('print');
ha = setupFigure(1,1, 'AxesWidth', 3.5, 'AxesHeight', 3.5,'SameAxes', true,...
    'XSpace', [1 0.85 0.5], 'YSpace', [1.5 0.5 0.5]);
set(ha, 'FontSize', 6);

title(ha(1), t, fset.lfont{:}, 'FontSize', 6);
ylabel(ha(1), '# nc82', fset.lfont{:}, 'FontSize', 6);
x = [1,3,7,9,13,15, 19,26];
y = [lMB, rMB, MB_nc82(1:2), MB_nc82(4:5), MB_nc82(3), MB_nc82(6)];
xticks(x)
xticklabels({'Gal4-L-WF', 'Gal4-R-WF','WT-L-F1', 'WT-R-F1', 'WT-L-F3', 'WT-R-F3', 'WT-L-F2', 'WT-R-RBPF1'})
bar(x,y, 'FaceColor', cfO, 'EdgeColor', ceO)
rotateXTickLabels(ha(1), 'Angle', 45, 'AdjustFigure', false);

f0 = gcf();
print(f0, [date 'number_nc82_perAnimal_combined.pdf'],'-dpdf', '-r0')