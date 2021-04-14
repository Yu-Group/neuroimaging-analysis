%% figures for the Expansion paper
%soma0
S0 = load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/2016_01_26_Sample4_Y10stitchedDeconvolved_final/Soma1/20160126_Soma1_holeParameters_rp3.mat');
hvS0 = [S0.params.holeVol];
%
%Axon1
A1 = load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/2016_03_30_Expansion_Sample3_YFP_MitoPosition2_TileStitched_z0p15/Axon1/20160330_Axon1_2ch_holeParameters_rp3.mat');
hvA1 = [A1.params.holeVol];
%
%Dendrite1
D1 = load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/2016_03_30_Expansion_Sample3_YFP_MitoPosition2_TileStitched_z0p15/Dendrite1/20160330_Dendrite1_2ch_holeParameters_rp3.mat');
hvD1 = [D1.params.holeVol];
%
% Soma1
S1 = load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/2016_03_30_Expansion_Sample3_YFP_MitoPosition2_TileStitched_z0p15/Soma1/20160330_Soma1_2ch_holeParameters_rp3.mat');
hvS1 = [S1.params.holeVol];

% Soma2 - 3ch
S2 = load('/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/20161211_ExM_YFP2_threeColor_z0p18/20161211_3ch_holeParameters_rp3.mat');
hvS2 = [S2.params.holeVol];
%%
% load('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\22-Jan-2018parameters_rp3_allholes.mat', 'holeVol')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MitoLysoSegmentations/22-Jan-2018parameters_rp3_allholes.mat', 'holeVol')



% size of all holes from axon, dendrite, soma
ef = 3.95; % expansion factor
% hv = log10([hvS0; hvA1; hvD1; hvS1; hvS2])/ef^3; %%%%%
fprintf('Median-allHoles: %d \n', median(holeVol/ef^3))
fprintf('MAD-allHoles: %d \n', mad(holeVol/ef^3,1))


hv = log10([holeVol]/ef^3);
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
%     'XSpace', [1.5 0.85 0.5], 'YSpace', [1 0.5 0.5]);
% set(ha, 'FontSize', 6);
% axes(ha)
%
% % histogram(hv, 50, 'FaceColor', cfK, 'FaceAlpha', 1, 'EdgeColor',ceK)
% % hold on
% histogram(log10([hvS0; hvS1; hvS2]), 50, 'FaceColor', cfB, 'FaceAlpha', 0.5, 'EdgeColor',ceB)
% % set(gca, 'xScale', 'log')
% xlim([-1.5 3]); ylim([0 350]);grid on
% xticklabels([10^-1,10^0, 10, 10^2, 10^3])
% ylabel('# holes');
% xlabel('volume (µm^3)');
% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',['SomasHoleVolume.eps']);

%
ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 0.85 0.5], 'YSpace', [1 0.5 0.5]);
set(ha, 'FontSize', 6);
axes(ha)

histogram(hv, 40, 'FaceColor', cfB, 'FaceAlpha', 1, 'EdgeColor',ceB)
% set(gca, 'xScale', 'log')
xlim([-5.5 3]);
% ylim([0 10]);
grid on
xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
ylabel('# holes');
xlabel('volume (µm^3)');
legend(['n = ' num2str(numel(holeVol))]);
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'AllHoleVolume.eps']);
%%
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MitoLysoSegmentations/22-Jan-2018parameters_rp3_allholes.mat')
MA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), rp3, 'unif',0);
MA = horzcat(MA{:})*.097/ef;
fprintf('MA Median-allHoles: %d \n', median(MA))
fprintf('MA MAD-allHoles: %d \n', mad(MA,1))

MA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), rp3, 'unif',0);
MA = horzcat(MA{:});
fprintf('AR Median-allHoles: %d \n', median(MA))
fprintf('AR MAD-allHoles: %d \n', mad(MA,1))



%% scatter plot of major and minor axes of the plots
% aspect ratio = major:minor
soma = [S0.rp3 S1.rp3 S2.rp3];
somaAR = (arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), soma))';
% somaAR(isinf(somaAR)) = nan;
somaVol = log10([hvS0; hvS1; hvS2]);

AxonAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), A1.rp3)';
AxonAR = AxonAR(1:26);
DendriteAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), D1.rp3)';
DendriteAR = DendriteAR(1:189);


fprintf('MA Median-allHoles: %d \n', median(MA))
fprintf('MA MAD-allHoles: %d \n', mad(MA,1))

fprintf('AR Median-allHoles: %d \n', median(MA))
fprintf('AR MAD-allHoles: %d \n', mad(MA,1))

%%

ms = 10;

ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
X = [somaVol(~isinf(somaAR))'; somaAR(~isinf(somaAR))'];
densityScatter(X, params);
colormap(jet)
xlim([-1.5 3]); ylim([0 20]);grid on
xticklabels([10^-1,10^0, 10, 10^2])
yticks(0:5:20)
ylabel('Aspect Ratio');
xlabel('expanded volume (µm^3)');
title('Soma')
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',['Soma_AR_vsHoleVolume.eps']);



ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
X = [log10(hvA1([1:21,23:26]))'; AxonAR([1:21,23:26])'];
densityScatter(X, params);
colormap(jet)
xlim([-1.5 3]); ylim([0 20]);grid on
xticklabels([10^-1,10^0, 10, 10^2])
yticks(0:5:20)
ylabel('Aspect Ratio');
xlabel('expanded volume (µm^3)');
title('Axon')
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',['Axon_AR_vsHoleVolume.eps']);



ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
X = [log10(hvD1(~isinf(DendriteAR)))'; DendriteAR(~isinf(DendriteAR))'];
densityScatter(X, params);
colormap(jet)
xlim([-1.5 3]); ylim([0 20]);grid on
xticklabels([10^-1,10^0, 10, 10^2])
yticks(0:5:20)
ylabel('Aspect Ratio');
xlabel('expanded volume (µm^3)');
title('Dendrite')
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',['Dendrite_AR_vsHoleVolume.eps']);

%% sort and select holes with mito and lyso Ab

% mito signal
X = [S1.params.intCh2Int]./[S1.params.holeVol]';
tS1_mito = thresholdOtsu(X);
figure, histogram(X,200), hold on
line([tS1_mito tS1_mito], [0 100])
idxS1_mito = X >tS1_mito;

X = [S2.params.intCh2Int]./[S2.params.holeVol]';
tS2_mito = thresholdOtsu(X);
figure, histogram(X,200), hold on
line([tS2_mito tS2_mito], [0 100])
idxS2_mito = X >tS2_mito;

X = [A1.params.intCh2Int(1:26)]./[A1.params.holeVol]';
tA1_mito = thresholdOtsu(X);
figure, histogram(X,10), hold on
line([tA1_mito tA1_mito], [0 100])
idxA1_mito = X >tA1_mito;

X = [D1.params.intCh2Int(1:189)]./[D1.params.holeVol]';
tD1_mito = thresholdOtsu(X);
figure, histogram(X,50), hold on
line([tD1_mito tD1_mito], [0 100])
idxD1_mito = X >tD1_mito;

% lyso signal

X = [S2.params.intCh3Int]./[S2.params.holeVol]';
tS2_lyso = thresholdOtsu(X);
figure, histogram(X,200), hold on
line([tS2_lyso tS2_lyso], [0 100])
idxS2_lyso = X >tS2_lyso;

%%

soma = [S1.rp3 S2.rp3];
somaAR = (arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), soma))';
somaAR_mito = somaAR([idxS1_mito idxS2_mito]);

somaAR = (arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), S2.rp3))';
somaAR_lyso = somaAR([idxS2_lyso]);

% somaAR(isinf(somaAR)) = nan;
somaVol_mito = log10([hvS1; hvS2]);
somaVol_mito = somaVol_mito([idxS1_mito idxS2_mito]);
somaVol_lyso = log10([hvS2]);
somaVol_lyso = somaVol_lyso(idxS2_lyso);

AxonAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), A1.rp3)';
AxonAR = AxonAR(1:26);
AxonAR = AxonAR(idxA1_mito);
AxonVol = log10(hvA1(idxA1_mito));

DendriteAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), D1.rp3)';
DendriteAR = DendriteAR(1:189);
DendriteAR = DendriteAR(idxD1_mito);
DendriteVol = log10(hvD1(idxD1_mito));

ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
X1 = [[somaVol_mito', AxonVol', DendriteVol']];
X2 = [somaAR_mito', AxonAR', DendriteAR'];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1; X2], params);
colormap(jet)
xlim([-1.5 3]); ylim([0 20]);grid on
xticklabels([10^-1,10^0, 10, 10^2])
yticks(0:5:20)
ylabel('Aspect Ratio');
xlabel('expanded volume (µm^3)');
title('Mitochondria')
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',['Mitochondria_AR_vsHoleVolume.eps']);



ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
X1 = [[somaVol_lyso']];
X2 = [somaAR_lyso'];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1;X2], params);
colormap(jet)
xlim([-1.5 3]); ylim([0 20]);grid on
xticklabels([10^-1,10^0, 10, 10^2])
yticks(0:5:20)
ylabel('Aspect Ratio');
xlabel('expanded volume (µm^3)');
title('Lysosome')
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',['lyso_AR_vsHoleVolume.eps']);

%% mito vs lysosome scatter plot
% load('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\parameters_rp3_lyso_mito_all.mat')
% load('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\22-Jan-2018parameters_rp3_lyso_mito_all.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MitoLysoSegmentations/22-Jan-2018parameters_rp3_lyso_mito_all.mat')
ef = 3.95; % expansion factor


mitoAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), mitorp3);
mitoMA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), mitorp3)*.097/ef;
mitoVol = log10(mitoholeVol/ef^3);
ms = 5;
ha = setupFigure(2,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};

X1 = [mitoVol];
X2 = [mitoMA];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1; X2], params);
colormap(jet)
xlim([-5.5 3]);
ylim([0 10]);
grid on
xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
yticks(0:2:10)
ylabel('Major Axis (µm)');
xlabel('volume (µm^3)');
title('Mitochondria')

fprintf('Major Axis \n')
fprintf('Median-mito: %d \n', median(X2))
fprintf('MAD-mito: %d \n', mad(X2,1))

axes(ha(2))
X1 = [mitoVol];
X2 = [mitoAR];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1; X2], params);
colormap(jet)
xlim([-5.5 3]);
ylim([0 30]);
grid on
xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
yticks(0:5:30)
ylabel('Aspect Ratio');
xlabel('volume (µm^3)');
legend(['n=' num2str(numel(X1))]);
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',['D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\' date 'Mitochondria_AR_vsHoleVolume.eps']);

fprintf('Volume \n')
fprintf('Median-mito: %d \n', 10^median(X1))
fprintf('MAD-mito: %d \n', mad(10.^X1,1))
fprintf('Aspect Ratio \n')
fprintf('Median-mito: %d \n', median(X2))
fprintf('MAD-mito: %d \n', mad(X2,1))
% lyso
clear idx
for j = 1:numel(lysorp3)
    idx(j) = ~isempty([lysorp3(j).MajorAxisLength]);
end

lysoAR = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), lysorp3(idx));
lysoMA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), lysorp3(idx))*.097/ef;

lysoVol = log10(lysoholeVol(idx)/ef^3);
ms = 5;
ha = setupFigure(2,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
axes(ha(1))
a = jet(256);
params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};

X1 = [lysoVol];
X2 = [lysoMA];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1; X2], params);
colormap(jet)
xlim([-5.5 3]);
ylim([0 10]);
grid on
xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
yticks(0:2:10)
ylabel('Major Axis (µm)');
xlabel('volume (µm^3)');
title('Lysosomes')
legend(['n=' num2str(numel(X1))]);
fprintf('Major Axis \n')
fprintf('Median-lyso: %d \n', median(X2))
fprintf('MAD-lyso: %d \n', mad(X2,1))


axes(ha(2))
X1 = [lysoVol];
X2 = [lysoAR];
X1= X1(~isinf(X2));
X2= X2(~isinf(X2));
densityScatter([X1; X2], params);
colormap(jet)
xlim([-5.5 3]);
ylim([0 30]);
grid on
xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
yticks(0:5:30)
ylabel('Major Axis (µm)');
xlabel('volume (µm^3)');

fprintf('volume \n')
fprintf('Median-lyso: %d \n', 10^median(X1))
fprintf('MAD-lyso: %d \n', mad(10.^X1,1))
fprintf('Aspect Ratio \n')
fprintf('Median-lyso: %d \n', median(X2))
fprintf('MAD-lyso: %d \n', mad(X2,1))
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',['D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM_YFP2_threeColor_z0p18\' date 'Lyso_AR_vsHoleVolume.eps']);

%% plots for all soma all dendrites
ef = 3.95; % expansion factor
% load('D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM\YFP2_threecolor\dataLoad_filtered.mat')
load('/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MitoLysoSegmentations/dataLoad_filtered.mat')

Axons = 2:4;
Dend = 5:9;
Soma = [10, 12:15];
cellType = {'Axons', 'Dendrites', 'Soma'};
cellTypeidx = {Axons, Dend, Soma};
for k = 3%1:3
    mitoMA = [];
    mitoAR = [];
    mitoVol = [];
    lysoMA = [];
    lysoAR = [];
    lysoVol = [];
    for i = cellTypeidx{k}
        d = load([data(i).source  'parameters_rp3_lyso_mito_all.mat']);
        alld(i) = d;
        clear idxL idxM
        for j = 1:numel(d.lysorp3)
            idxL(j) = ~isempty([d.lysorp3(j).MajorAxisLength]);
        end
        
        for j = 1:numel(d.mitorp3)
            idxM(j) = ~isempty([d.mitorp3(j).MajorAxisLength]);
        end
        
        mitoMA = [mitoMA arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), d.mitorp3(idxM))*.097/ef];
        mitoAR = [mitoAR arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), d.mitorp3(idxM))];
        mitoVol = [mitoVol log10(d.mitoholeVol(idxM)/ef^3)];
        
        lysoMA = [lysoMA arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), d.lysorp3(idxL))*.097/ef];
        lysoAR = [lysoAR arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength])./min([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), d.lysorp3(idxL))];
        lysoVol = [lysoVol  log10(d.lysoholeVol(idxL)/ef^3)];
    end
    cellType{k}
    ms = 5;
    ha = setupFigure(2,2, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
        'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
    
    a = jet(256);
    params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
    
    axes(ha(1))
    X1 = [mitoVol];
    X2 = [mitoMA];
    X1= X1(~isinf(X2));
    X2= X2(~isinf(X2));
    densityScatter([X1; X2], params);
    colormap(jet)
    xlim([-5.5 3]);
    ylim([0 10]);
    grid on
    xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
    yticks(0:2:10)
    ylabel('Major Axis (µm)');
    xlabel('volume (µm^3)');
    title('Mitochondria')
    legend(['n=' num2str(numel(X1))]);
    
    fprintf('Major Axis \n')
    fprintf('Median-mito: %d \n', median(X2))
    fprintf('MAD-mito: %d \n', mad(X2,1))
    fprintf('Volume \n')
    fprintf('Median-mito: %d \n', 10^median(X1))
    fprintf('MAD-mito: %d \n', mad(10.^X1,1))
    
    axes(ha(2))
    X1 = [lysoVol];
    X2 = [lysoMA];
    X1= X1(~isinf(X2));
    X2= X2(~isinf(X2));
    densityScatter([X1; X2], params);
    colormap(jet)
    xlim([-5.5 3]);
    ylim([0 10]);
    grid on
    xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
    yticks(0:2:10)
    %     ylabel('Major Axis (µm)');
    xlabel('volume (µm^3)');
    title('Lysosomes')
    legend(['n=' num2str(numel(X1))]);
    
    fprintf('Major Axis \n')
    fprintf('Median-lyso: %d \n', median(X2))
    fprintf('MAD-lyso: %d \n', mad(X2,1))
    fprintf('Volume \n')
    fprintf('Median-lyso: %d \n', 10^median(X1))
    fprintf('MAD-lyso: %d \n', mad(10.^X1,1))
    
    axes(ha(3))
    X1 = [mitoVol];
    X2 = [mitoAR];
    X1= X1(~isinf(X2));
    X2= X2(~isinf(X2));
    densityScatter([X1; X2], params);
    colormap(jet)
    xlim([-5.5 3]);
    ylim([0 30]);
    grid on
    xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
    yticks(0:5:30)
    ylabel('Aspect Ratio');
    xlabel('volume (µm^3)');
    %     title('Mitochondria')
    %     legend(['n=' num2str(numel(X1))]);
    fprintf('Aspect Ratio \n')
    fprintf('Median-mito: %d \n', median(X2))
    fprintf('MAD-mito: %d \n', mad(X2,1))
    
    
    axes(ha(4))
    X1 = [lysoVol];
    X2 = [lysoAR];
    X1= X1(~isinf(X2));
    X2= X2(~isinf(X2));
    densityScatter([X1; X2], params);
    colormap(jet)
    xlim([-5.5 3]);
    ylim([0 30]);
    grid on
    xticklabels({'10^{-4}','10^{-2}', '10^0', '10^2'})
    yticks(0:5:30)
    %     ylabel('Major Axis (µm)');
    xlabel('volume (µm^3)');
    %     title('Lysosomes')
    %     legend(['n=' num2str(numel(X1))]);
    fprintf('Aspect Ratio \n')
    fprintf('Median-lyso: %d \n', median(X2))
    fprintf('MAD-lyso: %d \n', mad(X2,1))
    
    f0 = gcf();
%     print(f0, '-painters','-depsc', '-loose',['D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\20161211_ExM\YFP2_threecolor\' cellType{k} 'MitoLyso_AR_vsHoleVolume.eps']);
end