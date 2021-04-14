
% this script is to check the distance mapping of each of the synapsed
% mapped in the alphalobes 0,1,2

% the xyz corrdinate are in pixels --> 8x8x8 nm
% Gokul Upadhyayula, 2017


% load synapse data
% load('/Users/Gokul/ImAnalysisScripts/GU_Repository/ExMScripts/FIBSEM_SynapseData.mat')
load('/Users/GU/ImAnalysisScripts/GU_Repository/ExMScripts/FIBSEM_SynapseData.mat')
% load('D:\Gokul\Software\GU_Repository\ExMScripts\FIBSEM_SynapseData.mat')

% pixel conversion
px = 8; %nm / pixel
fibsem_synapse.x = fibsem_synapse.x * px;
fibsem_synapse.y = fibsem_synapse.y * px;
fibsem_synapse.z = fibsem_synapse.z * px;


% index synapses in lobe
lobe0idx = fibsem_synapse.loben == 0;
lobe1idx = fibsem_synapse.loben == 1;
lobe2idx = fibsem_synapse.loben == 2;
lobe3idx = fibsem_synapse.loben == 3;

figure,
scatter3(fibsem_synapse.x(lobe0idx), fibsem_synapse.y(lobe0idx), fibsem_synapse.z(lobe0idx), 5, 'filled' ), hold on
scatter3(fibsem_synapse.x(lobe1idx), fibsem_synapse.y(lobe1idx), fibsem_synapse.z(lobe1idx), 5, 'filled' ), hold on
scatter3(fibsem_synapse.x(lobe2idx), fibsem_synapse.y(lobe2idx), fibsem_synapse.z(lobe2idx), 5, 'filled' ), hold on
scatter3(fibsem_synapse.x(lobe3idx), fibsem_synapse.y(lobe3idx), fibsem_synapse.z(lobe3idx), 5, 'filled' )
axis equal
view(3)
%% calulate volume of alpha 2
figure,
scatter3(fibsem_synapse.x(lobe3idx), fibsem_synapse.y(lobe3idx), fibsem_synapse.z(lobe3idx), 5, 'filled' ), hold on
%% alpha shape of alpha lobe 3
ha = setupFigure(2,5, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.25 1.25 1.25], 'YSpace', [1.5 1.25 1]);


rad = [250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000];

for i = 1:numel(rad)
shp = alphaShape(fibsem_synapse.x(lobe3idx), fibsem_synapse.y(lobe3idx), fibsem_synapse.z(lobe3idx),rad(i));
axes(ha(i))
plot(shp), axis equal
view(3)
v(i) = shp.volume * (10^-3)^3; % in um3
d(i) = sum(lobe3idx)/v(i); % density in # per um3
title(['Alpha Radius:' num2str(rad(i))]);
end

ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.75 1.25 1.25], 'YSpace', [1.5 1.25 1]);
 plot(rad, v, 'x-')
 xlabel('AlphaShape Radius')
 ylabel('Volume (\mum^3)')
%  set(gca, 'XScale', 'log')

%% alpha shape of alpha lobe 1,2,3
ha = setupFigure(2,5, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.25 1.25 1.25], 'YSpace', [1.5 1.25 1]);


% rad = [250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000];
rad = [2000, 4000, 8000, 16000];

for i = 1:numel(rad)
shp = alphaShape(fibsem_synapse.x(~lobe0idx), fibsem_synapse.y(~lobe0idx), fibsem_synapse.z(~lobe0idx),rad(i));
axes(ha(i))
plot(shp), axis equal
view(3)
v(i) = shp.volume * (10^-3)^3; % in um3
d(i) = sum(~lobe0idx)/v(i); % density in # per um3
title(['Alpha Radius:' num2str(rad(i))]);
end

ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.75 1.25 1.25], 'YSpace', [1.5 1.25 1]);
 plot(rad, v, 'x-')
 xlabel('AlphaShape Radius')
 ylabel('Volume (\mum^3)')
%  set(gca, 'XScale', 'log')
%%
lobes = {lobe1idx,lobe2idx,lobe3idx};
for i = 1:3
    X = [fibsem_synapse.x(lobes{i})'; fibsem_synapse.y(lobes{i})'; fibsem_synapse.z(lobes{i})']';
    Mdl = KDTreeSearcher(X);
    [~,D] = knnsearch(Mdl, X, 'k',2);
    distances{i} = D(:,2);
end
%%
ha = setupFigure(5,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.25 0.85 1.25], 'YSpace', [1.5 0.85 1]);
t = {'Alpha Lobe 1','Alpha Lobe 2','Alpha Lobe 3'};

i = 1;
axes(ha(i))
GU_stairsHistPlot({distances{i}}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true,'PercentileEnd', 100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
yyaxis left
ylabel('Relative Frequency')
title(t{i})
legend(['n=' num2str(numel(distances{i}))])
xticks([0:250:1000]);

i = 2;
axes(ha(i))
GU_stairsHistPlot({[],distances{i}}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true,'PercentileEnd', 100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
yyaxis left
ylabel('Relative Frequency')
title(t{i})
legend(['n=' num2str(numel(distances{i}))])
xticks([0:250:1000]);

i = 3;
axes(ha(i))
GU_stairsHistPlot({[], [], distances{i}}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true, 'PercentileEnd', 100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
yyaxis left
ylabel('Relative Frequency')
title(t{i})
legend(['n=' num2str(numel(distances{i}))])
xticks([0:250:1000]);

i = 4;
axes(ha(i))
GU_stairsHistPlot(distances, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true, 'PercentileEnd', 100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
yyaxis left
ylabel('Relative Frequency')
xticks([0:250:1000]);

i = 5;
axes(ha(i))
GU_stairsHistPlot({[],[],[],[],vertcat(distances{:})}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true, 'PercentileEnd', 100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
yyaxis left
ylabel('Relative Frequency')
xlabel('Distance to Closest Synapse (nm)')
xticks([0:250:1000]);

f0 = gcf();
print(f0, '-painters','-dpdf', '-loose',[date ' FIBSEM_SynapseDistance.pdf']);

%% DAN neurons
DANidx = {'PPL1-06-A';'PPL1-06-B'};
for j = 1:numel(DANidx)
    idx = fibsem_synapse.KC == DANidx{j};
    allidx = fibsem_synapse.KC == DANidx{1} | fibsem_synapse.KC == DANidx{2};
%     X = [fibsem_synapse.x(~allidx)'; fibsem_synapse.y(~allidx)'; fibsem_synapse.z(~allidx)']'; % non dan
    X = [fibsem_synapse.x'; fibsem_synapse.y'; fibsem_synapse.z']'; % all
    S = [fibsem_synapse.x(idx)'; fibsem_synapse.y(idx)'; fibsem_synapse.z(idx)']'; % non dan
    Mdl = KDTreeSearcher(X);
    [~,D] = knnsearch(Mdl, S, 'k',2);
    distances{j} = D(:,2);
end


%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(1,1, 'AxesWidth',10, 'AxesHeight', 10,'SameAxes', false,...
    'XSpace', [1.25 0.85 1.25], 'YSpace', [1.5 0.85 1]);


i = 1;
axes(ha(i))
title(ha(i), {'Blue:Alpha1, Red:Alpha2, Green:Alpha3'})
scatter3(fibsem_synapse.x(lobe1idx)/1000, fibsem_synapse.y(lobe1idx)/1000, fibsem_synapse.z(lobe1idx)/1000, 1, 'filed', 'MarkerFaceColor', cfB ), hold on
scatter3(fibsem_synapse.x(lobe2idx)/1000, fibsem_synapse.y(lobe2idx)/1000, fibsem_synapse.z(lobe2idx)/1000, 1, 'filed', 'MarkerFaceColor', [0.2 0.9 0.9] ), hold on
scatter3(fibsem_synapse.x(lobe3idx)/1000, fibsem_synapse.y(lobe3idx)/1000, fibsem_synapse.z(lobe3idx)/1000, 1, 'filed', 'MarkerFaceColor',  ceB)
axis equal
ylim([10 60])
zlim([20 100])
xlim([10 60])
zticks(0:20:100)
yticks(0:20:60)
xticks(0:20:60)
xlabel('um')
ylabel('um')
zlabel('um')
view(3)
grid on
set(gca, 'Zdir', 'reverse')
set(gca, 'Xdir', 'reverse')
set(gca, 'Ydir', 'reverse')
zt = get(gca, 'ZTick');
set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))


xt = get(gca, 'XTick');
set(gca, 'XTick',zt, 'XTickLabel',fliplr(xt))
yt = get(gca, 'YTick');
set(gca, 'YTick',zt, 'YTickLabel',fliplr(yt))
% i = 2;
% axes(ha(i))
% scatter3(fibsem_synapse.x(allidx), fibsem_synapse.y(allidx), fibsem_synapse.z(allidx), 1, 'MarkerFaceColor', ceK ), hold on
% scatter3(fibsem_synapse.x(lobe1idx), fibsem_synapse.y(lobe1idx), fibsem_synapse.z(lobe1idx), 0.5, 'filed', 'MarkerFaceColor', cfB ), hold on
% scatter3(fibsem_synapse.x(lobe2idx), fibsem_synapse.y(lobe2idx), fibsem_synapse.z(lobe2idx), 0.5, 'filed', 'MarkerFaceColor', ceB ), hold on
% scatter3(fibsem_synapse.x(lobe3idx), fibsem_synapse.y(lobe3idx), fibsem_synapse.z(lobe3idx), 0.5, 'filed', 'MarkerFaceColor', [0.2 0.9 0.9] )
% axis equal
% ylim([2 6e4])
% zlim([1 11e4])
% xlim([2 5e4])
% xlabel('nm')
% ylabel('nm')
% zlabel('nm')
% view(3)
% grid on
f0 = gcf();
print(f0, '-painters','-dpdf', '-loose',[date ' FIBSEM_DANSynapseLocation.pdf']);

ha = setupFigure(2,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.25 0.85 1.25], 'YSpace', [1.5 0.85 1]);
t = DANidx;

i = 1;
axes(ha(i))
GU_stairsHistPlot({[],[],[],distances{i}}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true, 'PercentileEnd',100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
xticks([0:250:1000]);
yyaxis left
ylabel('Relative Frequency')
legend( ['n = ' num2str(numel(distances{i}))])
title(t{i})

i = 2;
axes(ha(i))
GU_stairsHistPlot({[],[],[],distances{i}}, 'BinSize', 10, 'CreateFigure', false, 'ShowECDF', true, 'PercentileEnd',100);
set(gca,'fontsize',6, 'FontName', 'Helvetica')
xlim([0 1000])
xticks([0:250:1000]);
yyaxis left
ylabel('Relative Frequency')
legend( ['n = ' num2str(numel(distances{i}))])
xlabel('Closest DAN Synapse Distance to All Synapse')
title(t{i})
f0 = gcf();
print(f0, '-painters','-dpdf', '-loose',[date ' FIBSEM_DANSynapseDistance.pdf']);
%% calculate density

lobe = lobe1idx;

AZdia = .2;% diameter of activezone
clusterSizeMin = 0;
sigmas = [1.3, 1.3];
ef = 4;
px = 0.097; %in um
bandwidth = AZdia/2/.097*ef;
V = 1;% in um3 surveyed volume
r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
re = r*ef; % expanded radius in um
bandwidth2 = re/px;

% second round of clustering to group
fprintf('clustering to calculate spot density (diameter -- 1um3)*expansion factor^3...')
tic
ptData(:,1) = fibsem_synapse.x(lobe);
ptData(:,2) = fibsem_synapse.y(lobe);
ptData(:,3) =  fibsem_synapse.z(lobe);
    
ClusCenters = cat(1, clusterInfo.ptClusterCenter);
[clusterInfo2, ~, ~] = MeanShiftClustering( ClusCenters, bandwidth2);
toc

ClusCenters2 = cat(1, clusterInfo2.ptClusterCenter);
cm = jet(numel([clusterInfo2.numPoints]));
cm = flip(cm,2);
[~,idx] = sort([clusterInfo2.numPoints], 'descend');